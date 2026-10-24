#!/usr/bin/env python3
"""Regenerate SHA-256 manifest after intentional edits to this release folder."""

import csv
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
excluded = {"manifest.csv", "validation/current_check.json"}
rows = []
for path in sorted(p for p in root.rglob("*") if p.is_file()):
    relative = path.relative_to(root).as_posix()
    if relative in excluded or "__pycache__" in path.parts:
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    rows.append((relative, path.stat().st_size, digest.hexdigest()))

with (root / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["path", "bytes", "sha256"])
    writer.writerows(rows)
print(f"Manifest updated: {len(rows)} files")
