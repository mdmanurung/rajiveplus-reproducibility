"""Independently authored processed-input audit. No fits or clinical inference."""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone

import requests
import yaml

PIN = "12cae2bd385919cb2ffb48abb07813be9c440fc3"
REPOSITORY = "https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy"


def check(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def json_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def load_config(root, name):
    return yaml.safe_load((Path(root) / "config" / f"{name}.yaml").read_text())


def covid_config(root):
    cfg = load_config(root, "datasets")["covid"]
    check(cfg["commit"] == PIN and cfg["repository"] == REPOSITORY, "Unexpected source pin")
    check(cfg["orientation"] == "features_by_samples", "Orientation must be features_by_samples")
    check(cfg["raw_counts"] is False and cfg["track"] == "processed_cell_type_means_and_bulk_omics",
          "Processed expression is NOT raw counts; count-aware analysis forbidden")
    check(len(cfg["files"]) == 7, "Exactly six matrices and one metadata file required")
    names = [x["name"] for x in cfg["files"]]
    check(len(set(names)) == 7 and all(Path(n).name == n for n in names), "Invalid/duplicate filenames")
    check(sum(x["kind"] == "expression" for x in cfg["files"]) == 4, "Four expression blocks required")
    for item in cfg["files"]:
        check(len(item["sha256"]) == 64 and all(c in "0123456789abcdef" for c in item["sha256"]),
              "Missing/invalid SHA256")
    return cfg


def verified_files(root, cfg):
    outputs = []
    for item in cfg["files"]:
        path = Path(root) / cfg["cache_dir"] / item["name"]
        check(path.is_file() and path.stat().st_size > 0, f"Missing/nonregular input: {item['name']}")
        digest = sha256(path)
        check(digest == item["sha256"], f"Checksum mismatch: {item['name']}")
        outputs.append({"name": item["name"], "sha256": digest, "bytes": path.stat().st_size})
    return outputs


def run_acquire(root, output):
    start = time.monotonic()
    cfg = covid_config(root)
    cache = Path(root) / cfg["cache_dir"]
    cache.mkdir(parents=True, exist_ok=True)
    for item in cfg["files"]:
        path = cache / item["name"]
        if path.exists():
            check(path.is_file() and sha256(path) == item["sha256"],
                  f"Checksum mismatch in cache; remove corrupt file explicitly: {item['name']}")
            continue
        url = (f"https://raw.githubusercontent.com/ByronSyun/DIVAS_COVID19_CaseStudy/"
               f"{cfg['commit']}/{cfg['path']}/{item['name']}")
        tmp = path.with_suffix(path.suffix + ".download")
        try:
            with requests.get(url, stream=True, timeout=(20, 120)) as response:
                response.raise_for_status()
                with tmp.open("wb") as f:
                    for chunk in response.iter_content(1024 * 1024):
                        f.write(chunk)
            check(sha256(tmp) == item["sha256"], f"Downloaded checksum mismatch: {item['name']}")
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
    json_write(output, {"status": "PASS", "source": REPOSITORY, "commit": PIN,
                        "source_path": cfg["path"], "files": verified_files(root, cfg),
                        "wall_time_s": time.monotonic() - start,
                        "data_license": cfg["data_license"], "redistributed": False})


def unique_nonempty(values, label):
    check(all(v and v == v.strip() for v in values), f"Empty/whitespace {label}")
    check(len(values) == len(set(values)), f"Duplicate {label}")


def metadata(path, cfg):
    with Path(path).open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        unique_nonempty(header, "metadata columns")
        rows = list(reader)
    check(all(len(x) == len(header) for x in rows), "Metadata row width")
    data = [dict(zip(header, row)) for row in rows]
    required = {"sample_id", "clinical_id", "timepoint", "wos", "divas_sample_id"}
    check(required <= set(header), "Missing required metadata columns")
    check(len(data) == cfg["expected_samples"], "Metadata sample row count")
    samples = [x["sample_id"] for x in data]
    unique_nonempty(samples, "metadata sample IDs")
    unique_nonempty([x["divas_sample_id"] for x in data], "DIVAS sample IDs")
    check(samples == [x["divas_sample_id"] for x in data], "Metadata DIVAS/sample ID disagreement")
    for row in data:
        check(row["clinical_id"] and row["clinical_id"] == row["clinical_id"].strip(), "Missing patient ID")
        check(row["timepoint"] in cfg["expected_timepoints"], "Unknown timepoint")
        check(row["sample_id"].endswith("_" + row["timepoint"]), "Sample/timepoint mismatch")
        try:
            value = float(row["wos"])
        except ValueError as exc:
            raise ValueError("Non-numeric WOS source recoding") from exc
        check(math.isfinite(value) and value in {0, 1, 1.5, 2, 3, 4, 5, 6, 7},
              "Invalid WOS source recoding")
    pairs = [(x["clinical_id"], x["timepoint"]) for x in data]
    check(len(pairs) == len(set(pairs)), "Duplicate patient/timepoint")
    groups = {}
    for row in data:
        groups.setdefault(row["clinical_id"], []).append(row)
    check(len(groups) == cfg["expected_patients"], "Patient count")
    check(sum(len(rows) == 2 for rows in groups.values()) == cfg["expected_paired_patients"],
          "Paired patient count")
    check(all(len(rows) in (1, 2) for rows in groups.values()), "Unexpected patient multiplicity")
    check({tp: sum(r["timepoint"] == tp for r in data) for tp in cfg["expected_timepoints"]}
          == cfg["expected_timepoints"], "Timepoint count")
    return data


def matrix(path, samples, expected_rows):
    features = set()
    noninteger = 0
    with Path(path).open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        unique_nonempty(header[1:], "matrix sample IDs")
        check(header[1:] == samples, "Exact sample order/orientation mismatch versus metadata")
        for row in reader:
            check(len(row) == len(samples) + 1, "Matrix row width")
            check(row[0] and row[0] == row[0].strip(), "Empty/whitespace feature ID")
            check(row[0] not in features, "Duplicate feature ID")
            features.add(row[0])
            try:
                values = [float(x) for x in row[1:]]
            except ValueError as exc:
                raise ValueError("Non-numeric/missing matrix value") from exc
            check(all(math.isfinite(x) for x in values), "Nonfinite matrix value")
            noninteger += sum(not x.is_integer() for x in values)
    check(len(features) == expected_rows, "Feature row count")
    return {"features": len(features), "samples": len(samples), "noninteger_values": noninteger,
            "duplicate_features": 0, "duplicate_samples": 0, "missing_values": 0,
            "nonfinite_values": 0, "exact_metadata_sample_order": True}


def manifests(root, rows):
    rank = {"T1": 1, "T2": 2}
    earliest = {}
    for row in rows:
        patient = row["clinical_id"]
        if patient not in earliest or rank[row["timepoint"]] < rank[earliest[patient]["timepoint"]]:
            earliest[patient] = row
    primary = [r for r in rows if r is earliest[r["clinical_id"]]]
    check(len(primary) == len(earliest), "Ambiguous earliest selection")
    directory = Path(root) / "data/manifests"
    directory.mkdir(parents=True, exist_ok=True)
    for arm, selected in [("all", rows), ("primary", primary), ("sensitivity", rows)]:
        path = directory / f"covid_{arm}.csv"
        with path.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["source_row_1based", "sample_id", "clinical_id", "timepoint",
                             "wos_source_recoding", "resampling_group"])
            positions = {r["sample_id"]: i + 1 for i, r in enumerate(rows)}
            for row in selected:
                writer.writerow([positions[row["sample_id"]], row["sample_id"], row["clinical_id"],
                                 row["timepoint"], row["wos"], row["clinical_id"]])
    return {"all": len(rows), "primary": len(primary), "sensitivity": len(rows),
            "split_unit": "clinical_id", "primary_rule": "earliest_available_T1_then_T2",
            "primary_order": "source_order_retained_after_selection"}


def validate_manifests(root, rows):
    """Reject hand-edited selection/grouping indexes, not merely hash them."""
    rank = {"T1": 1, "T2": 2}
    earliest = {}
    for row in rows:
        patient = row["clinical_id"]
        if patient not in earliest or rank[row["timepoint"]] < rank[earliest[patient]["timepoint"]]:
            earliest[patient] = row
    columns = ["source_row_1based", "sample_id", "clinical_id", "timepoint",
               "wos_source_recoding", "resampling_group"]
    for arm in ("all", "primary", "sensitivity"):
        expected = []
        for i, row in enumerate(rows, 1):
            if arm == "primary" and row is not earliest[row["clinical_id"]]:
                continue
            expected.append([str(i), row["sample_id"], row["clinical_id"], row["timepoint"],
                             row["wos"], row["clinical_id"]])
        with (Path(root) / f"data/manifests/covid_{arm}.csv").open(newline="") as f:
            actual = list(csv.reader(f))
        check(actual == [columns] + expected,
              f"Invalid {arm} manifest: source/order/earliest/grouping mismatch")


def run_audit(root, output):
    start = time.monotonic()
    cfg = covid_config(root)
    files = verified_files(root, cfg)
    cache = Path(root) / cfg["cache_dir"]
    rows = metadata(cache / cfg["metadata_file"], cfg)
    samples = [r["sample_id"] for r in rows]
    summaries = {}
    for item in cfg["files"]:
        if item["kind"] != "metadata":
            summaries[item["name"]] = matrix(cache / item["name"], samples, item["rows"])
            if item["kind"] == "expression":
                check(summaries[item["name"]]["noninteger_values"] > 0,
                      "Expected processed noninteger expression; raw/count track mismatch")
    arms = manifests(root, rows)
    validate_manifests(root, rows)
    json_write(output, {"status": "PASS", "source_commit": PIN, "files": files,
                        "sample_count": len(rows), "patient_count": cfg["expected_patients"],
                        "paired_patient_count": cfg["expected_paired_patients"],
                        "timepoints": cfg["expected_timepoints"], "blocks": summaries,
                        "sample_order_sha256": hashlib.sha256(json.dumps(samples).encode()).hexdigest(),
                        "manifests": arms, "track": cfg["track"], "raw_counts": False,
                        "wos_1_5_meaning": cfg["wos_1_5_meaning"],
                        "selection_limitation": cfg["selection_limitation"],
                        "wall_time_s": time.monotonic() - start,
                        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024})


def git(root, *args):
    proc = subprocess.run(["git", "-C", str(Path(root).resolve()), *args],
                          text=True, capture_output=True)
    return proc.stdout.strip() if proc.returncode == 0 else "UNAVAILABLE_SOURCE_ARCHIVE"


def identity_hashes(root, paths):
    return {str(p.relative_to(root)): sha256(p) for p in sorted(paths)}


def run_receipt(root, output):
    root = Path(root)
    audit_path = root / "results/manifests/audit.json"
    audit = json.loads(audit_path.read_text())
    check(audit["status"] == "PASS", "Input audit has not passed")
    cfg = covid_config(root)
    verified_files(root, cfg)
    rows = metadata(root / cfg["cache_dir"] / cfg["metadata_file"], cfg)
    validate_manifests(root, rows)
    source_paths = [root / "scripts/stage0.py", root / "workflow/Snakefile", root / "requirements.txt"]
    manifest_paths = sorted((root / "data/manifests").glob("covid_*.csv"))
    check(len(manifest_paths) == 3, "Expected three local manifests")
    json_write(output, {
        "status": "PASS", "scope": "processed_acquisition_audit_provenance_only",
        "analysis_set": "pilot", "scientific_readiness": "NOT_READY", "fits": "NOT_RUN",
        "inference": "NOT_RUN", "seed": "NOT_APPLICABLE_NO_RANDOM_OPERATIONS",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(), "source_commit": PIN,
        "paper_git_commit": git(root, "rev-parse", "HEAD"),
        "paper_git_branch": git(root, "branch", "--show-current"),
        "paper_tracked_changes": git(root, "status", "--porcelain", "--untracked-files=no"),
        "config_sha256": identity_hashes(root, list((root / "config").glob("*.yaml"))),
        "source_sha256": identity_hashes(root, source_paths),
        "input_sha256": {x["name"]: x["sha256"] for x in audit["files"]},
        "local_manifest_sha256": identity_hashes(root, manifest_paths),
        "audit_sha256": sha256(audit_path),
        "acquisition_sha256": sha256(root / "results/manifests/acquisition.json"),
        "environment": {"python": sys.version, "platform": platform.platform(),
                        "executor": "local", "slurm": "UNAVAILABLE_NOT_USED",
                        "versions": {name: importlib.metadata.version(name)
                                     for name in ("snakemake", "pytest", "requests", "PyYAML")},
                        "thread_environment": {k: os.environ.get(k, "UNSET") for k in
                                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}},
        "exit_code": 0, "audit_wall_time_s": audit["wall_time_s"],
        "audit_peak_rss_mb": audit["peak_rss_mb"],
        "limitations": [cfg["selection_limitation"], cfg["wos_1_5_meaning"],
                        "No new raw-count preprocessing; no fit, calibration, or clinical inference.",
                        "Data terms not inferred from upstream code license."]})


def registry(root):
    cfg = load_config(root, "figures")
    check(cfg["main"] == [f"Fig{i}" for i in range(1, 6)], "Exactly five main figures required")
    check(cfg["supplementary"] == [f"S{i}" for i in range(1, 11)], "Exactly ten standard supplements required")
    check(cfg["conditional_supplementary"] == ["S11", "S12"], "Conditional supplements must be S11/S12")
    check(cfg["tables"] == [f"TableS{i}" for i in range(1, 9)], "Exactly eight supplementary tables required")
    return cfg


def run_readiness(root, output):
    cfg = registry(root)
    records = []
    for group in ("main", "supplementary", "conditional_supplementary", "tables"):
        for identifier in cfg[group]:
            prereq = cfg["prerequisite_pattern"].format(id=identifier)
            records.append({"id": identifier, "group": group, "status": "NOT_READY",
                            "analysis_set": "pilot", "scientific_prerequisite": prereq,
                            "prerequisite_present": (Path(root) / prereq).is_file(),
                            "renderer": "NOT_IMPLEMENTED",
                            "conditional_inclusion": "NOT_APPROVED" if group == "conditional_supplementary" else "NOT_APPLICABLE"})
    json_write(output, {"scientific_readiness": "NOT_READY", "fits": "NOT_RUN",
                        "figures_and_tables": records})


def require_figure(root, identifier):
    cfg = registry(root)
    figures = cfg["main"] + cfg["supplementary"] + cfg["conditional_supplementary"]
    check(identifier in figures, f"Unregistered scientific figure: {identifier}")
    prerequisite = Path(root) / cfg["prerequisite_pattern"].format(id=identifier)
    check(prerequisite.is_file(),
          f"Missing scientific prerequisites for {identifier}: {prerequisite}; validated fits/metrics/design NOT_RUN")
    raise RuntimeError(f"{identifier}: scientific renderer NOT_IMPLEMENTED; no substitute/fake PDF is generated")
