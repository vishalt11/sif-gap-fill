from pathlib import Path
from zipfile import ZipFile


TILE_ROOT = Path(__file__).resolve().parent / "data" / "geodes_wasp_zips" / "32UMC"


def is_safe_zip_member(target_dir: Path, member_name: str) -> bool:
    target_dir = target_dir.resolve()
    member_path = (target_dir / member_name).resolve()
    return member_path == target_dir or target_dir in member_path.parents


def unzip_then_delete(zip_path: Path) -> None:
    year_dir = zip_path.parent

    print(f"Extracting: {zip_path}")
    with ZipFile(zip_path) as archive:
        unsafe_members = [
            name for name in archive.namelist()
            if not is_safe_zip_member(year_dir, name)
        ]
        if unsafe_members:
            raise RuntimeError(f"Unsafe paths in {zip_path}: {unsafe_members[:5]}")

        archive.extractall(year_dir)

    zip_path.unlink()
    print(f"Deleted zip: {zip_path}")


def main() -> None:
    if not TILE_ROOT.is_dir():
        raise FileNotFoundError(f"Tile folder does not exist: {TILE_ROOT}")

    zip_count = 0
    for year_dir in sorted(path for path in TILE_ROOT.iterdir() if path.is_dir()):
        print(f"\nYear: {year_dir.name}")
        year_zips = sorted(year_dir.glob("*.zip"))

        if not year_zips:
            print("  No zip files found.")
            continue

        for zip_path in year_zips:
            unzip_then_delete(zip_path)
            zip_count += 1

    print(f"\nDone. Extracted and deleted {zip_count} zip file(s).")


if __name__ == "__main__":
    main()
