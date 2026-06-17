from __future__ import annotations

import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import urllib3
from pygeodes import Config, Geodes


urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

CONFIG_FILE = Path("config.json")
SIF_CSV = Path("data/ba_sif_mgrs_crop_composition.csv")

COLLECTION_ID = "THEIA_REFLECTANCE_SENTINEL2_L3A"
TILE = "32UPU"
GRID_CODE = f"T{TILE}"

BASE_OUT_DIR = Path("data/geodes_wasp_zips") / TILE
CHECKLIST_PATH = BASE_OUT_DIR / f"geodes_{TILE}_download_checklist.csv"

RECHECK_MISSING = False
CURL_RETRIES = "20"
CURL_RETRY_DELAY_SECONDS = "10"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def monthly_token(year: int, month: int) -> str:
    return f"{year}{month:02d}15"


def monthly_datetime(year: int, month: int) -> str:
    return f"{year}-{month:02d}-15T00:00:00Z"


def year_month_parts(year_month: str) -> tuple[int, int]:
    year, month = year_month.split("-")
    return int(year), int(month)


def item_identifier(item) -> str:
    return item.find("identifier") or ""


def configure_geodes() -> tuple[Geodes, Config]:
    geodes = Geodes()
    conf = Config.from_file(str(CONFIG_FILE))
    conf.download_dir = str(BASE_OUT_DIR.resolve())
    conf.checksum_error = False
    conf.use_async_requests = False
    geodes.set_conf(conf)

    if not conf.api_key:
        raise RuntimeError(f"No api_key found in {CONFIG_FILE}")

    return geodes, conf


def target_months_from_sif() -> list[str]:
    sif = pd.read_csv(SIF_CSV)

    if "mgrs_tile" in sif.columns:
        tile_values = (
            sif["mgrs_tile"]
            .astype(str)
            .str.strip()
            .str.upper()
            .str.replace("^T", "", regex=True)
        )
        sif = sif.loc[tile_values.eq(TILE)].copy()

    if sif.empty:
        available_tiles = pd.read_csv(SIF_CSV, usecols=["mgrs_tile"])
        tile_counts = (
            available_tiles["mgrs_tile"]
            .astype(str)
            .str.strip()
            .str.upper()
            .str.replace("^T", "", regex=True)
            .value_counts()
            .head(20)
        )
        raise RuntimeError(
            f"No rows found for tile {TILE} in {SIF_CSV}. "
            f"Top available tiles:\n{tile_counts.to_string()}"
        )

    sif["Delta_Date"] = pd.to_datetime(sif["Delta_Date"], errors="coerce")
    invalid_dates = sif["Delta_Date"].isna().sum()
    if invalid_dates:
        print(f"Warning: {invalid_dates} Delta_Date values could not be parsed and will be ignored.")

    months = sorted(sif["Delta_Date"].dt.strftime("%Y-%m").dropna().unique())

    if not months:
        raise RuntimeError(f"No Delta_Date month-year values found for {TILE} in {SIF_CSV}")

    return months


def load_checklist(target_months: list[str]) -> pd.DataFrame:
    base_rows = pd.DataFrame(
        {
            "year_month": target_months,
            "target_date": [monthly_datetime(*year_month_parts(ym)) for ym in target_months],
            "tile": TILE,
            "grid_code": GRID_CODE,
            "collection_id": COLLECTION_ID,
        }
    )

    if CHECKLIST_PATH.exists():
        existing = pd.read_csv(CHECKLIST_PATH)
        checklist = base_rows.merge(
            existing.drop(columns=["target_date", "tile", "grid_code", "collection_id"], errors="ignore"),
            on="year_month",
            how="left",
        )
    else:
        checklist = base_rows

    defaults = {
        "item_found": False,
        "item_id": None,
        "identifier": None,
        "archive_title": None,
        "download_url": None,
        "status": "pending",
        "output_path": None,
        "file_size_bytes": None,
        "error": None,
        "started_at": None,
        "completed_at": None,
    }

    for column, default in defaults.items():
        if column not in checklist.columns:
            checklist[column] = default

    checklist["status"] = checklist["status"].fillna("pending")
    return checklist


def save_checklist(checklist: pd.DataFrame) -> None:
    CHECKLIST_PATH.parent.mkdir(parents=True, exist_ok=True)
    checklist.to_csv(CHECKLIST_PATH, index=False)


def completed_file_is_valid(row: pd.Series) -> bool:
    output_path = row.get("output_path")
    if not isinstance(output_path, str) or not output_path:
        return False

    path = Path(output_path)
    if not path.exists() or path.stat().st_size <= 0:
        return False

    recorded_size = row.get("file_size_bytes")
    if pd.isna(recorded_size):
        return True

    try:
        return path.stat().st_size == int(recorded_size)
    except (TypeError, ValueError):
        return True


def should_skip(row: pd.Series) -> bool:
    status = row.get("status")

    if status == "completed" and completed_file_is_valid(row):
        return True

    if status == "missing" and not RECHECK_MISSING:
        return True

    return False


def matching_items(items, token: str):
    expected_bits = [token, f"_T{TILE}_"]
    return [
        item
        for item in items
        if all(bit in item_identifier(item) for bit in expected_bits)
    ]


def find_monthly_tile_item(geodes: Geodes, year: int, month: int):
    date_time = monthly_datetime(year, month)
    token = monthly_token(year, month)

    query = {
        "grid:code": {"eq": GRID_CODE},
        "start_datetime": {"eq": date_time},
    }

    try:
        items = geodes.search_items(
            query=query,
            get_all=True,
            return_df=False,
            quiet=True,
            collections=[COLLECTION_ID],
        )
    except Exception as exc:
        print(f"  exact start_datetime query failed: {exc}")
        items = []

    matches = matching_items(items, token)
    if matches:
        return matches[0]

    print("  falling back to grid-only query and local identifier filtering")
    items = geodes.search_items(
        query={"grid:code": {"eq": GRID_CODE}},
        get_all=True,
        return_df=False,
        quiet=True,
        collections=[COLLECTION_ID],
    )
    matches = matching_items(items, token)

    if not matches:
        return None

    if len(matches) > 1:
        print("  multiple matching items found; using the first")
        for item in matches:
            print(f"    {item.id} | {item_identifier(item)}")

    return matches[0]


def output_path_for_item(item, year: int) -> Path:
    archive_title = item.data_asset.title or f"{item_identifier(item)}.zip"
    return (BASE_OUT_DIR / str(year) / archive_title).resolve()


def curl_path() -> str:
    found = shutil.which("curl.exe") or shutil.which("curl")
    if not found:
        raise RuntimeError("curl was not found on PATH")
    return found


def curl_download(download_url: str, output_path: Path, api_key: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        curl_path(),
        "-k",
        "-L",
        "-C",
        "-",
        "--retry",
        CURL_RETRIES,
        "--retry-delay",
        CURL_RETRY_DELAY_SECONDS,
        "--retry-all-errors",
        "-H",
        f"X-API-Key: {api_key}",
        "-X",
        "GET",
        download_url,
        "--output",
        str(output_path),
    ]

    subprocess.run(cmd, check=True)


def update_row(checklist: pd.DataFrame, idx: int, **values) -> None:
    for key, value in values.items():
        checklist.loc[idx, key] = value


def download_one_month(geodes: Geodes, conf: Config, checklist: pd.DataFrame, idx: int) -> None:
    row = checklist.loc[idx]
    year_month = row["year_month"]
    year, month = year_month_parts(year_month)

    print(f"\n{year_month} / {TILE}")

    if should_skip(row):
        print(f"  skipping: {row.get('status')} -> {row.get('output_path')}")
        return

    update_row(
        checklist,
        idx,
        status="running",
        error=None,
        started_at=utc_now(),
        completed_at=None,
    )
    save_checklist(checklist)

    try:
        item = find_monthly_tile_item(geodes, year, month)

        if item is None:
            print("  missing on GEODES")
            update_row(
                checklist,
                idx,
                item_found=False,
                status="missing",
                error=None,
                completed_at=utc_now(),
            )
            return

        identifier = item_identifier(item)
        outfile = output_path_for_item(item, year)
        download_url = item.data_asset.href

        update_row(
            checklist,
            idx,
            item_found=True,
            item_id=item.id,
            identifier=identifier,
            archive_title=item.data_asset.title,
            download_url=download_url,
            output_path=str(outfile),
        )
        save_checklist(checklist)

        existing_size = outfile.stat().st_size if outfile.exists() else 0
        print(f"  archive: {item.data_asset.title}")
        print(f"  output: {outfile}")
        if existing_size > 0:
            print(f"  resuming from existing file size: {existing_size} bytes")

        curl_download(download_url, outfile, conf.api_key)

        file_size = outfile.stat().st_size if outfile.exists() else None
        print(f"  completed: {outfile} ({file_size} bytes)")

        update_row(
            checklist,
            idx,
            status="completed",
            output_path=str(outfile),
            file_size_bytes=file_size,
            error=None,
            completed_at=utc_now(),
        )

    except subprocess.CalledProcessError as exc:
        print(f"  curl failed with exit code {exc.returncode}")
        update_row(
            checklist,
            idx,
            status="failed",
            error=f"curl exit code {exc.returncode}",
            completed_at=utc_now(),
        )
    except Exception as exc:
        print(f"  failed: {exc}")
        update_row(
            checklist,
            idx,
            status="failed",
            error=repr(exc),
            completed_at=utc_now(),
        )
    finally:
        save_checklist(checklist)


def main() -> None:
    BASE_OUT_DIR.mkdir(parents=True, exist_ok=True)
    target_months = target_months_from_sif()
    checklist = load_checklist(target_months)
    save_checklist(checklist)

    geodes, conf = configure_geodes()

    print(f"Downloading GEODES archives for {TILE} with curl resume")
    print(f"Months: {len(target_months)}")
    print(f"Checklist: {CHECKLIST_PATH}")
    print(f"Output base: {BASE_OUT_DIR.resolve()}")

    for idx in checklist.index:
        download_one_month(geodes, conf, checklist, idx)

    final = pd.read_csv(CHECKLIST_PATH)
    print("\nFinal status counts")
    print(final["status"].value_counts(dropna=False).to_string())
    print(f"\nWrote checklist: {CHECKLIST_PATH}")


if __name__ == "__main__":
    main()
