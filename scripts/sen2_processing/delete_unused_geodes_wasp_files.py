"""Remove Sentinel-2 files that are not needed for the selected indices.

The script keeps the six bands used by NDMI, NDVI, EVI, NIRv, and NDRE:
B2, B4, B5, B8, B8A, and B11. It marks all other FRC band TIFFs for
deletion. It also marks WGT_R1 and WGT_R2 TIFFs inside MASKS directories.

The default mode is a dry run. Pass --delete to perform the deletions.
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


DEFAULT_ROOT = Path(__file__).resolve().parent / "data" / "geodes_wasp_zips"

# GEODES filenames omit the leading zero used in Sentinel-2 band notation:
# B02 -> B2, B04 -> B4, B05 -> B5, and B08/B08A -> B8/B8A.
KEEP_BANDS = frozenset({"B2", "B4", "B5", "B8", "B8A", "B11"})

BAND_FILE_PATTERN = re.compile(
    r"_FRC_(B(?:8A|\d{1,2}))\.tif$",
    flags=re.IGNORECASE,
)
WGT_MASK_PATTERN = re.compile(
    r"_(WGT_R[12])\.tif$",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class DeletionCandidate:
    path: Path
    reason: str
    size_bytes: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Delete unused GEODES Sentinel-2 band TIFFs and WGT_R1/WGT_R2 "
            "mask TIFFs. The default is a dry run."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help=f"Dataset root to scan (default: {DEFAULT_ROOT})",
    )
    parser.add_argument(
        "--delete",
        action="store_true",
        help="Actually delete the listed files. Without this flag, only preview them.",
    )
    parser.add_argument(
        "--list-all",
        action="store_true",
        help="List every candidate instead of only the first 20.",
    )
    return parser.parse_args()


def deletion_reason(path: Path) -> str | None:
    """Return why a TIFF should be deleted, or None when it should be kept."""
    band_match = BAND_FILE_PATTERN.search(path.name)
    if band_match:
        band = band_match.group(1).upper()
        if band not in KEEP_BANDS:
            return f"unused band {band}"
        return None

    if path.parent.name.casefold() == "masks":
        mask_match = WGT_MASK_PATTERN.search(path.name)
        if mask_match:
            return f"mask {mask_match.group(1).upper()}"

    return None


def find_candidates(root: Path) -> list[DeletionCandidate]:
    candidates: list[DeletionCandidate] = []

    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.casefold() != ".tif":
            continue

        reason = deletion_reason(path)
        if reason is not None:
            candidates.append(
                DeletionCandidate(
                    path=path,
                    reason=reason,
                    size_bytes=path.stat().st_size,
                )
            )

    return sorted(candidates, key=lambda item: str(item.path).casefold())


def format_bytes(size_bytes: int) -> str:
    size = float(size_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.2f} {unit}"
        size /= 1024

    raise AssertionError("Unreachable")


def print_summary(candidates: list[DeletionCandidate]) -> None:
    counts = Counter(candidate.reason for candidate in candidates)

    print("Deletion candidates by type:")
    for reason, count in sorted(counts.items()):
        print(f"  {reason}: {count}")

    total_size = sum(candidate.size_bytes for candidate in candidates)
    print(f"Total: {len(candidates)} file(s), {format_bytes(total_size)}")


def print_candidates(
    candidates: list[DeletionCandidate],
    root: Path,
    list_all: bool,
) -> None:
    visible = candidates if list_all else candidates[:20]

    print("\nFiles:")
    for candidate in visible:
        print(f"  [{candidate.reason}] {candidate.path.relative_to(root)}")

    hidden_count = len(candidates) - len(visible)
    if hidden_count:
        print(f"  ... {hidden_count} more (use --list-all to show every file)")


def delete_candidates(candidates: list[DeletionCandidate]) -> int:
    deleted_count = 0
    failures: list[tuple[Path, OSError]] = []

    for candidate in candidates:
        try:
            candidate.path.unlink()
            deleted_count += 1
        except OSError as error:
            failures.append((candidate.path, error))

    print(f"\nDeleted {deleted_count} of {len(candidates)} candidate file(s).")

    if failures:
        print("Failed deletions:")
        for path, error in failures:
            print(f"  {path}: {error}")
        return 1

    return 0


def main() -> int:
    args = parse_args()
    root = args.root.expanduser().resolve()

    if not root.is_dir():
        raise FileNotFoundError(f"Dataset root does not exist or is not a directory: {root}")

    print(f"Dataset root: {root}")
    print(f"Mode: {'DELETE' if args.delete else 'DRY RUN'}")
    print(f"Bands retained: {', '.join(sorted(KEEP_BANDS))}\n")

    candidates = find_candidates(root)
    if not candidates:
        print("No matching files found. Nothing to delete.")
        return 0

    print_summary(candidates)
    print_candidates(candidates, root, args.list_all)

    if not args.delete:
        print(
            f"\nDry run only: no files were deleted. "
            f"Run {Path(__file__).name} --delete to apply these deletions."
        )
        return 0

    return delete_candidates(candidates)


if __name__ == "__main__":
    raise SystemExit(main())
