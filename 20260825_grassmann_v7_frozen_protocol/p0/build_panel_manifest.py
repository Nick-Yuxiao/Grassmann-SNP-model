from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def line_count(path: Path) -> int | None:
    if path.suffix.lower() not in {".txt", ".tsv", ".csv", ".json", ".yaml", ".yml"}:
        return None
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create an auditable T01 data-branch manifest without copying data.")
    parser.add_argument("--branch", choices=["A", "B"], required=True)
    parser.add_argument("--input", action="append", type=Path, required=True, help="Repeat for every panel/sample/variant/config input.")
    parser.add_argument("--approval-reference", required=True, help="Non-secret approval record ID for branch A, or synthetic-design record ID for branch B.")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    resolved = [path.resolve() for path in args.input]
    missing = [str(path) for path in resolved if not path.is_file()]
    if missing:
        raise SystemExit("Missing T01 inputs: " + ", ".join(missing))
    if args.branch == "B" and not any("synt" in path.name.lower() for path in resolved):
        raise SystemExit("Branch B requires an explicitly named synthetic-generator/config input.")

    records = []
    for path in resolved:
        stat = path.stat()
        records.append({
            "logical_name": path.name,
            "bytes": stat.st_size,
            "sha256": sha256(path),
            "line_count_if_text": line_count(path),
        })
    total_bytes = sum(item["bytes"] for item in records)
    payload = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "branch": args.branch,
        "claim_limit": (
            "approved individual-level biobank scope" if args.branch == "A"
            else "1KGP+HGDP plus explicitly synthetic augmentation; no real-biobank claim"
        ),
        "approval_or_design_reference": args.approval_reference,
        "inputs": records,
        "total_bytes": total_bytes,
        "status": "SIGNED_INPUTS_PRESENT",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    panel_path = args.output_dir / "PANEL_MANIFEST.json"
    decision_path = args.output_dir / "DATA_BRANCH_DECISION.md"
    panel_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    decision_path.write_text(
        "# T01 data branch decision\n\n"
        f"- Branch: **{args.branch}**\n"
        f"- Approval/design reference: `{args.approval_reference}`\n"
        f"- Claim limit: {payload['claim_limit']}\n"
        f"- Inputs: {len(records)} files, {total_bytes} bytes\n"
        f"- PANEL_MANIFEST SHA-256: `{sha256(panel_path)}`\n\n"
        "This record freezes input identities, not scientific outcomes. Any input change requires a new manifest.\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PASS", "branch": args.branch, "inputs": len(records), "total_bytes": total_bytes}, indent=2))


if __name__ == "__main__":
    main()
