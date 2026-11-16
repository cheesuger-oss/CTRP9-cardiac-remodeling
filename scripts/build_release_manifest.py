#!/usr/bin/env python3
"""Create or verify the public-release SHA-256 manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "release_manifest.json"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def inventory() -> list[dict[str, object]]:
    rows = []
    for path in sorted(ROOT.rglob("*")):
        if (
            not path.is_file()
            or path == MANIFEST
            or ".git" in path.parts
            or "__pycache__" in path.parts
            or path.suffix == ".pyc"
        ):
            continue
        rows.append({
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": digest(path),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    current = inventory()
    if args.verify:
        stored = json.loads(MANIFEST.read_text(encoding="utf-8"))
        ok = stored["files"] == current
        print(f"RESULT: {'PASS' if ok else 'FAIL'} ({len(current)} files)")
        return 0 if ok else 1
    payload = {"release": "v0.1.0", "files": current}
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"WROTE {MANIFEST} ({len(current)} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
