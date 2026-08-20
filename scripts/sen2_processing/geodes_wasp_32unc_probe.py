import os
import re
import shutil
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests
import urllib3
from pygeodes import Config, Geodes


urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

SIF_CSV = Path("data/ns_sif_mgrs_crop_composition.csv")
OUT_DIR = Path("data/geodes_wasp_probe")
TILE = "32UNC"
GEODES_GRID_CODE = f"T{TILE}"

# First modeling/download focus.
BANDS = ["B8", "B4", "B11"]
ASSET_KEYS = [f"FRC_{band}" for band in BANDS]

# Keep the first download small in scope: one exact monthly product, one band.
DO_ONE_DOWNLOAD_TEST = True
DOWNLOAD_TEST_BAND = "B8"

# If you already know the GEODES collection id, put it here.
# Otherwise the script searches collections for "wasp" and chooses a likely
# Sentinel-2 L3A WASP collection.
COLLECTION_ID_OVERRIDE = "THEIA_REFLECTANCE_SENTINEL2_L3A"

# Optional API key. You can also set GEODES_API_KEY in your environment.
GEODES_API_KEY = os.environ.get("GEODES_API_KEY")
CONFIG_FILE = Path("config.json")


def month_token(year_month: str) -> str:
    return year_month.replace("-", "") + "15"


def month_window(year_month: str) -> tuple[str, str]:
    start = pd.Timestamp(f"{year_month}-01")
    end = start + pd.DateOffset(months=1) - pd.Timedelta(seconds=1)
    return (
        start.strftime("%Y-%m-%dT00:00:00Z"),
        end.strftime("%Y-%m-%dT23:59:59Z"),
    )


def item_tile(item_id: str) -> str | None:
    match = re.search(r"_T([0-9]{2}[A-Z]{3})_", item_id)
    return match.group(1) if match else None


def item_text(item) -> str:
    parts = [item.id]
    for key in [
        "dataset",
        "identifier",
        "product",
        "product:type",
        "datetime",
        "start_datetime",
        "end_datetime",
        "grid:code",
    ]:
        value = item.find(key)
        if value is not None:
            parts.append(str(value))

    for asset_key, asset in item.assets.items():
        parts.extend([asset_key, asset.title or "", asset.href or ""])

    return " ".join(parts)


def load_target_months() -> list[str]:
    sif = pd.read_csv(SIF_CSV, parse_dates=["Delta_Date"])

    if "mgrs_tile" in sif.columns:
        sif = sif.loc[sif["mgrs_tile"].astype(str).eq(TILE)].copy()

    months = sorted(sif["Delta_Date"].dt.strftime("%Y-%m").unique())

    if not months:
        raise RuntimeError(f"No SIF months found for tile {TILE} in {SIF_CSV}")

    return months


def choose_wasp_collection(geodes: Geodes) -> str:
    if COLLECTION_ID_OVERRIDE:
        return COLLECTION_ID_OVERRIDE

    collections = geodes.search_collections(return_df=False, quiet=True)

    collection_rows = []
    for collection in collections:
        title = collection.title or ""
        description = collection.description or ""
        text = f"{collection.id} {title} {description}".lower()

        positive_terms = [
            "sentinel",
            "sentinel-2",
            "sentinel2",
            "s2",
            "l3a",
            "wasp",
            "reflectance",
            "theia",
        ]
        negative_terms = ["venus", "landsat", "pleiades", "spot", "sar", "radar"]

        score = sum(3 if term in {"sentinel", "sentinel-2", "sentinel2", "s2"} else 1
                    for term in positive_terms if term in text)
        score -= sum(5 for term in negative_terms if term in text)

        collection_rows.append(
            {
                "collection_id": collection.id,
                "title": title,
                "description": description,
                "score": score,
            }
        )

    collections_df = pd.DataFrame(collection_rows).sort_values(
        ["score", "collection_id"],
        ascending=[False, True],
    )
    collections_df.to_csv(OUT_DIR / "geodes_wasp_collection_candidates.csv", index=False)

    if collections_df.empty:
        raise RuntimeError("GEODES returned no collections.")

    plausible = collections_df.loc[
        collections_df["collection_id"].str.contains("SENTINEL|S2|THEIA", case=False, na=False)
        | collections_df["title"].str.contains("SENTINEL|S2|THEIA", case=False, na=False)
    ].copy()
    plausible = plausible.loc[
        ~plausible["collection_id"].str.contains("VENUS", case=False, na=False)
        & ~plausible["title"].str.contains("VENUS", case=False, na=False)
    ]
    plausible = plausible.sort_values(["score", "collection_id"], ascending=[False, True])

    print("GEODES Sentinel/L3A collection candidates:")
    print(plausible.head(20)[["collection_id", "title", "score"]].to_string(index=False))

    if plausible.empty or plausible.iloc[0]["score"] <= 0:
        print("\nTop 20 collections overall:")
        print(collections_df.head(20)[["collection_id", "title", "score"]].to_string(index=False))
        raise RuntimeError(
            "Could not confidently choose a Sentinel-2 L3A/WASP collection. "
            "Open data/geodes_wasp_probe/geodes_wasp_collection_candidates.csv, "
            "find the correct collection id, then set COLLECTION_ID_OVERRIDE."
        )

    chosen = plausible.iloc[0]["collection_id"]
    print(f"\nUsing collection: {chosen}\n")
    return chosen


def search_exact_month_items(
    geodes: Geodes,
    collection_id: str,
    year_month: str,
) -> list:
    start, end = month_window(year_month)

    query = {
        "datetime": {"gte": start, "lte": end},
        "grid:code": {"eq": GEODES_GRID_CODE},
    }

    try:
        items = geodes.search_items(
            query=query,
            get_all=True,
            return_df=False,
            quiet=True,
            collections=[collection_id],
        )
    except Exception as exc:
        print(f"Query with grid:code failed for {year_month}: {exc}")
        print("Retrying date-only query and filtering locally.")
        items = geodes.search_items(
            query={"datetime": {"gte": start, "lte": end}},
            get_all=True,
            return_df=False,
            quiet=True,
            collections=[collection_id],
        )

    expected_date = month_token(year_month)

    return [
        item
        for item in items
        if item.find("grid:code") == GEODES_GRID_CODE
        and expected_date in item_text(item)
    ]


def diagnostic_probe(geodes: Geodes, collection_id: str, target_months: list[str]):
    print(f"\nDiagnostic probe for {collection_id} / {GEODES_GRID_CODE}")

    probe_rows = []

    probe_queries = [
        {
            "label": "grid_only",
            "query": {"grid:code": {"eq": GEODES_GRID_CODE}},
        },
        {
            "label": "first_month_date_only",
            "query": {
                "datetime": {
                    "gte": month_window(target_months[0])[0],
                    "lte": month_window(target_months[0])[1],
                }
            },
        },
        {
            "label": "first_month_grid_and_date",
            "query": {
                "datetime": {
                    "gte": month_window(target_months[0])[0],
                    "lte": month_window(target_months[0])[1],
                },
                "grid:code": {"eq": GEODES_GRID_CODE},
            },
        },
    ]

    for probe in probe_queries:
        try:
            items = geodes.search_items(
                query=probe["query"],
                get_all=False,
                page=1,
                return_df=False,
                quiet=True,
                collections=[collection_id],
            )
        except Exception as exc:
            print(f"  {probe['label']}: query failed: {exc}")
            continue

        print(f"  {probe['label']}: returned {len(items)} item(s)")

        for item in items[:10]:
            probe_rows.append(
                {
                    "probe": probe["label"],
                    "item_id": item.id,
                    "dataset": item.find("dataset"),
                    "identifier": item.find("identifier"),
                    "datetime": item.find("datetime"),
                    "start_datetime": item.find("start_datetime"),
                    "end_datetime": item.find("end_datetime"),
                    "grid_code": item.find("grid:code"),
                    "asset_keys": ";".join(item.assets.keys()),
                    "text_excerpt": item_text(item)[:500],
                }
            )

    if probe_rows:
        probe_path = OUT_DIR / f"geodes_{TILE}_diagnostic_probe.csv"
        pd.DataFrame(probe_rows).to_csv(probe_path, index=False)
        print(f"  wrote diagnostic probe: {probe_path}\n")


def iter_asset_records(item) -> Iterable[dict]:
    for key, asset in item.assets.items():
        yield {
            "asset_key": key,
            "asset_title": asset.title,
            "asset_href": asset.href,
            "asset_roles": "|".join(asset.roles or []),
            "asset_type": asset.media_type,
        }


def asset_matches_band(asset_key: str, asset_title: str, band: str) -> bool:
    needles = [
        f"FRC_{band}".upper(),
        f"_{band}.".upper(),
        f"_{band}_".upper(),
        f"{band}.TIF".upper(),
    ]
    haystack = f"{asset_key} {asset_title}".upper()
    return any(needle in haystack for needle in needles)


def item_band_availability(item, bands: list[str]) -> dict:
    asset_records = list(iter_asset_records(item))
    out = {}

    for band in bands:
        matching_assets = [
            record
            for record in asset_records
            if asset_matches_band(record["asset_key"], record["asset_title"], band)
        ]
        out[f"{band}_asset_count"] = len(matching_assets)
        out[f"{band}_asset_keys"] = ";".join(record["asset_key"] for record in matching_assets)
        out[f"{band}_asset_hrefs"] = ";".join(record["asset_href"] for record in matching_assets)

    return out


def safe_list_archive_files(geodes: Geodes, item) -> list[str]:
    try:
        return geodes.list_item_files(item)
    except Exception as exc:
        print(f"Could not list archive files for {item.id}: {exc}")
        return []


def archive_file_matches_band(filename: str, band: str) -> bool:
    upper = filename.upper()
    return (
        f"FRC_{band}".upper() in upper
        or f"_{band}.TIF" in upper
        or f"_{band}_" in upper
    )


def add_archive_band_availability(geodes: Geodes, rows: list[dict], items_by_id: dict):
    for row in rows:
        item = items_by_id.get(row["item_id"])
        if item is None:
            continue

        files = safe_list_archive_files(geodes, item)
        row["archive_file_count"] = len(files)

        for band in BANDS:
            matches = [filename for filename in files if archive_file_matches_band(filename, band)]
            row[f"{band}_archive_file_count"] = len(matches)
            row[f"{band}_archive_files"] = ";".join(matches)


def download_one_band_file(geodes: Geodes, availability: pd.DataFrame, items_by_id: dict):
    count_col = f"{DOWNLOAD_TEST_BAND}_archive_file_count"
    if count_col not in availability.columns:
        print(
            f"No {count_col} column found because no downloadable archive files "
            "were discovered; skipping one-file download test."
        )
        return

    candidates = availability.loc[
        availability[count_col].fillna(0).astype(int).gt(0)
    ].copy()

    if candidates.empty:
        print(
            f"No archive file found for {DOWNLOAD_TEST_BAND}; "
            "skipping one-file download test."
        )
        return

    row = candidates.sort_values(["year_month", "item_id"]).iloc[0]
    item = items_by_id[row["item_id"]]
    filename = row[f"{DOWNLOAD_TEST_BAND}_archive_files"].split(";")[0]

    print("\nOne-file download test")
    print(f"  item: {item.id}")
    print(f"  band: {DOWNLOAD_TEST_BAND}")
    print(f"  archive member: {filename}")

    before_files = set(OUT_DIR.rglob("*"))
    geodes.download_item_files(item, filenames=[filename])
    after_files = set(OUT_DIR.rglob("*"))

    new_files = sorted(after_files - before_files)

    if not new_files:
        print("Download command finished, but no new file was detected under OUT_DIR.")
        return

    downloaded = new_files[0]
    final_path = OUT_DIR / f"{item.id}_{DOWNLOAD_TEST_BAND}{Path(downloaded).suffix}"
    shutil.move(str(downloaded), final_path)
    print(f"  saved as: {final_path}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    geodes = Geodes()

    if CONFIG_FILE.exists():
        conf = Config.from_file(str(CONFIG_FILE))
    else:
        conf = Config(api_key=GEODES_API_KEY)

    conf.download_dir = str(OUT_DIR.resolve())
    conf.checksum_error = False
    conf.use_async_requests = False
    geodes.set_conf(conf)

    target_months = load_target_months()
    collection_id = choose_wasp_collection(geodes)
    diagnostic_probe(geodes, collection_id, target_months)

    rows = []
    items_by_id = {}

    for year_month in target_months:
        print(f"Querying {year_month} / {TILE}")
        items = search_exact_month_items(geodes, collection_id, year_month)

        if not items:
            rows.append(
                {
                    "year_month": year_month,
                    "mgrs_tile": TILE,
                    "item_id": None,
                    "collection_id": collection_id,
                    "item_found": False,
                }
            )
            continue

        for item in items:
            items_by_id[item.id] = item
            row = {
                "year_month": year_month,
                "mgrs_tile": TILE,
                "grid_code": item.find("grid:code"),
                "item_id": item.id,
                "collection_id": collection_id,
                "item_found": True,
                "datetime": item.datetime.isoformat() if item.datetime else None,
                "geodes_datetime": item.find("datetime"),
                "start_datetime": item.find("start_datetime"),
                "end_datetime": item.find("end_datetime"),
                "dataset": item.find("dataset"),
                "identifier": item.find("identifier"),
                "asset_keys": ";".join(item.assets.keys()),
            }
            row.update(item_band_availability(item, BANDS))
            rows.append(row)

    availability = pd.DataFrame(rows)

    # GEODES often provides products as archives; this step checks whether
    # B8/B4/B11 are present inside the archive even if they are not direct assets.
    add_archive_band_availability(geodes, rows, items_by_id)
    availability = pd.DataFrame(rows)

    availability_path = OUT_DIR / f"geodes_{TILE}_wasp_band_availability.csv"
    availability.to_csv(availability_path, index=False)

    print("\nAvailability summary")
    summary_cols = [
        "year_month",
        "item_found",
        "item_id",
        "B8_archive_file_count",
        "B4_archive_file_count",
        "B11_archive_file_count",
        "B8_asset_count",
        "B4_asset_count",
        "B11_asset_count",
    ]
    present_summary_cols = [col for col in summary_cols if col in availability.columns]
    print(availability[present_summary_cols].to_string(index=False))
    print(f"\nWrote: {availability_path}")

    if DO_ONE_DOWNLOAD_TEST:
        download_one_band_file(geodes, availability, items_by_id)


if __name__ == "__main__":
    main()
