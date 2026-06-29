from __future__ import annotations

import csv
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


BASE_URL = "https://glass.hku.hk/archive/NDVI/MODIS/250M/"
PRODUCT_PREFIX = "GLASS13D01.V10"
PRODUCTION_SUFFIX_BY_YEAR = {
    2019: "2023064",
    2020: "2023064",
    2021: "2023075",
    2022: "2024177",
    2023: "2025012",
    2024: "2026057",
}

YEARS = range(2019, 2025)
START_DOY = 33
END_DOY = 209
DOY_STEP = 8

TILES = ["h18v03", "h18v04"]

OUTPUT_DIR = Path("data/glass_ndvi_modis_250m")
CHECKLIST_PATH = OUTPUT_DIR / "glass_ndvi_download_checklist.csv"

OVERWRITE = False

REQUEST_TIMEOUT_SECONDS = 120
REQUEST_RETRIES = 3
REQUEST_RETRY_SLEEP_SECONDS = 5
USER_AGENT = "Mozilla/5.0 GLASS-NDVI-downloader"


@dataclass(frozen=True)
class Target:
    year: int
    doy: int
    tile: str

    @property
    def date(self):
        return datetime.strptime(f"{self.year}{self.doy:03d}", "%Y%j").date()

    @property
    def day_dir_url(self) -> str:
        return f"{BASE_URL}{self.year}/{self.doy:03d}/"

    @property
    def filename(self) -> str:
        production_suffix = PRODUCTION_SUFFIX_BY_YEAR[self.year]
        return f"{PRODUCT_PREFIX}.A{self.year}{self.doy:03d}.{self.tile}.{production_suffix}.hdf"

    @property
    def file_url(self) -> str:
        return f"{self.day_dir_url}{self.filename}"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def target_days() -> list[Target]:
    return [
        Target(year=year, doy=doy, tile=tile.lower())
        for year in YEARS
        for doy in range(START_DOY, END_DOY + 1, DOY_STEP)
        for tile in TILES
    ]


def request_bytes(url: str) -> bytes:
    last_error: Exception | None = None

    for attempt in range(1, REQUEST_RETRIES + 1):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                return response.read()
        except (HTTPError, URLError, TimeoutError) as error:
            last_error = error
            if isinstance(error, HTTPError) and error.code == 404:
                raise
            if attempt < REQUEST_RETRIES:
                time.sleep(REQUEST_RETRY_SLEEP_SECONDS)

    assert last_error is not None
    raise last_error


def local_path_for(target: Target) -> Path:
    return OUTPUT_DIR / target.tile / str(target.year) / target.filename


def file_is_complete(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def download_file(url: str, output_path: Path) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_suffix(output_path.suffix + ".part")

    content = request_bytes(url)
    temp_path.write_bytes(content)
    temp_path.replace(output_path)

    return output_path.stat().st_size


def load_checklist() -> dict[tuple[str, int, int], dict[str, str]]:
    if not CHECKLIST_PATH.exists():
        return {}

    with CHECKLIST_PATH.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return {
            (
                row["tile"],
                int(row["year"]),
                int(row["doy"]),
            ): row
            for row in reader
        }


def save_checklist(rows: dict[tuple[str, int, int], dict[str, str]]) -> None:
    CHECKLIST_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tile",
        "year",
        "doy",
        "date",
        "status",
        "remote_url",
        "output_path",
        "file_size_bytes",
        "error",
        "updated_at",
    ]

    with CHECKLIST_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow(rows[key])


def checklist_row(
    target: Target,
    status: str,
    output_path: Path | None = None,
    file_size_bytes: int | None = None,
    error: str | None = None,
) -> dict[str, str]:
    return {
        "tile": target.tile,
        "year": str(target.year),
        "doy": f"{target.doy:03d}",
        "date": target.date.isoformat(),
        "status": status,
        "remote_url": target.file_url,
        "output_path": str(output_path) if output_path else "",
        "file_size_bytes": str(file_size_bytes) if file_size_bytes is not None else "",
        "error": error or "",
        "updated_at": utc_now(),
    }


def should_skip(row: dict[str, str] | None) -> bool:
    if OVERWRITE or row is None:
        return False

    if row.get("status") != "completed":
        return False

    output_path = row.get("output_path")
    return bool(output_path) and file_is_complete(Path(output_path))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    checklist = load_checklist()
    targets = target_days()

    print(
        f"Downloading GLASS NDVI targets: {len(targets)} files, "
        f"years {min(YEARS)}-{max(YEARS)}, DOY {START_DOY:03d}-{END_DOY:03d} step {DOY_STEP}, "
        f"tiles {', '.join(TILES)}"
    )
    print(f"Output directory: {OUTPUT_DIR}")

    for index, target in enumerate(targets, start=1):
        print(
            f"[{index}/{len(targets)}] {target.date.isoformat()} "
            f"DOY {target.doy:03d} {target.tile}"
        )

        key = (target.tile, target.year, target.doy)
        existing_row = checklist.get(key)
        if should_skip(existing_row):
            continue

        output_path = local_path_for(target)
        if file_is_complete(output_path) and not OVERWRITE:
            checklist[key] = checklist_row(
                target,
                status="completed",
                output_path=output_path,
                file_size_bytes=output_path.stat().st_size,
            )
            save_checklist(checklist)
            continue

        try:
            file_size = download_file(target.file_url, output_path)
            checklist[key] = checklist_row(
                target,
                status="completed",
                output_path=output_path,
                file_size_bytes=file_size,
            )
        except Exception as error:
            checklist[key] = checklist_row(
                target,
                status="download_error",
                output_path=output_path,
                error=repr(error),
            )

        save_checklist(checklist)

    print(f"Done. Checklist written to: {CHECKLIST_PATH}")


if __name__ == "__main__":
    main()
