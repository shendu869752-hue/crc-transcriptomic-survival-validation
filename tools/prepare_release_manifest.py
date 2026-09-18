from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reanalysis_v3" / "manifests" / "release_manifest_sha256.csv"
REQUIRED = [
    ROOT / "README.md",
    ROOT / "LICENSE",
    ROOT / "CITATION.cff",
    ROOT / ".zenodo.json",
    ROOT / "run_v3_local.ps1",
    ROOT / "run_tests.ps1",
    ROOT / "reanalysis_v3" / "scripts" / "01_core_reanalysis_v3.R",
    ROOT / "reanalysis_v3" / "results" / "run_status.csv",
    ROOT / "reanalysis_v3" / "results" / "code_manifest.csv",
    ROOT / "reanalysis_v3" / "results" / "input_manifest.csv",
]
TEXT_EXTENSIONS = {".R", ".ps1", ".py", ".md", ".json", ".cff", ".txt", ".csv", ".gitignore"}
SECRET_PATTERNS = {
    "GitHub token": re.compile(r"gh[opsu]_[A-Za-z0-9]{20,}"),
    "generic private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "password assignment": re.compile(r"(?i)password\s*[=:]\s*['\"][^'\"]+['\"]"),
    "API key assignment": re.compile(r"(?i)api[_-]?key\s*[=:]\s*['\"][^'\"]+['\"]"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


errors: list[str] = []
for path in REQUIRED:
    if not path.is_file():
        errors.append(f"Missing required file: {path.relative_to(ROOT)}")

with (ROOT / ".zenodo.json").open("r", encoding="utf-8") as handle:
    zenodo = json.load(handle)
if zenodo.get("version") != "1.0.0":
    errors.append(".zenodo.json version must be 1.0.0")
if zenodo.get("upload_type") != "software":
    errors.append(".zenodo.json upload_type must be software")

rows = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or ".git" in path.parts or path == MANIFEST:
        continue
    relative = path.relative_to(ROOT).as_posix()
    if path.stat().st_size >= 100_000_000:
        errors.append(f"File exceeds GitHub 100 MB limit: {relative}")
    if path.suffix in TEXT_EXTENSIONS or path.name in {"CITATION.cff", ".gitignore"}:
        content = path.read_text(encoding="utf-8", errors="replace")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(content):
                errors.append(f"Potential {label} in {relative}")
    rows.append((relative, path.stat().st_size, sha256(path)))

MANIFEST.parent.mkdir(parents=True, exist_ok=True)
with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["relative_path", "bytes", "sha256"])
    writer.writerows(rows)

print(f"FILES={len(rows)}")
print(f"BYTES={sum(row[1] for row in rows)}")
print(f"MANIFEST={MANIFEST.relative_to(ROOT).as_posix()}")
if errors:
    print("RELEASE_AUDIT_FAIL")
    for error in errors:
        print(error)
    raise SystemExit(1)
print("RELEASE_AUDIT_PASS")
