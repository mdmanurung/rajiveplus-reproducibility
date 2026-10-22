"""Small synthetic plumbing fixtures only; never used as scientific results."""
import csv
from pathlib import Path
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import stage0 as s


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        csv.writer(f).writerows(rows)
    return path


@pytest.fixture
def fixture_metadata(tmp_path):
    header = ["sample_id", "clinical_id", "timepoint", "wos", "divas_sample_id"]
    rows = [["toy_A_T2", "toy_A", "T2", "1.5", "toy_A_T2"],
            ["toy_B_T2", "toy_B", "T2", "3", "toy_B_T2"],
            ["toy_A_T1", "toy_A", "T1", "1", "toy_A_T1"]]
    cfg = {"expected_samples": 3, "expected_patients": 2, "expected_paired_patients": 1,
           "expected_timepoints": {"T1": 1, "T2": 2}}
    return header, rows, cfg


def test_primary_is_earliest_not_first_and_sensitivity_grouped(tmp_path, fixture_metadata):
    header, rows, cfg = fixture_metadata
    data = s.metadata(write_csv(tmp_path / "meta.csv", [header] + rows), cfg)
    result = s.manifests(tmp_path, data)
    assert result["all"] == 3 and result["primary"] == 2 and result["sensitivity"] == 3
    with (tmp_path / "data/manifests/covid_primary.csv").open() as f:
        primary = list(csv.DictReader(f))
    assert [r["sample_id"] for r in primary] == ["toy_B_T2", "toy_A_T1"]
    with (tmp_path / "data/manifests/covid_sensitivity.csv").open() as f:
        sensitivity = list(csv.DictReader(f))
    assert [r["sample_id"] for r in sensitivity] == [r[0] for r in rows]
    assert all(r["clinical_id"] == r["resampling_group"] for r in sensitivity)
    assert sensitivity[0]["wos_source_recoding"] == "1.5"
    s.validate_manifests(tmp_path, data)


@pytest.mark.parametrize("mutation", ["drop", "group", "order", "duplicate"])
def test_modified_manifests_rejected(tmp_path, fixture_metadata, mutation):
    header, rows, cfg = fixture_metadata
    data = s.metadata(write_csv(tmp_path / "meta.csv", [header] + rows), cfg)
    s.manifests(tmp_path, data)
    path = tmp_path / "data/manifests/covid_sensitivity.csv"
    with path.open() as f:
        contents = list(csv.reader(f))
    if mutation == "drop":
        contents.pop()
    elif mutation == "group":
        contents[1][-1] = "wrong_patient"
    elif mutation == "order":
        contents[1], contents[2] = contents[2], contents[1]
    elif mutation == "duplicate":
        contents[2] = contents[1]
    write_csv(path, contents)
    with pytest.raises(ValueError, match="Invalid sensitivity manifest"):
        s.validate_manifests(tmp_path, data)


@pytest.mark.parametrize("mutation,message", [
    ("duplicate_sample", "Duplicate metadata sample"),
    ("divas_mismatch", "DIVAS/sample"),
    ("duplicate_pair", "Duplicate patient/timepoint"),
    ("unknown_time", "Unknown timepoint"),
    ("nonfinite_wos", "Invalid WOS"),
    ("missing_patient", "Missing patient"),
    ("patient_count", "Patient count"),
    ("row_count", "sample row count"),
    ("duplicate_header", "Duplicate metadata columns"),
])
def test_bad_metadata_fails(tmp_path, fixture_metadata, mutation, message):
    header, rows, cfg = fixture_metadata
    if mutation == "duplicate_sample":
        rows[1][0] = rows[0][0]
    elif mutation == "divas_mismatch":
        rows[0][4] = "wrong"
    elif mutation == "duplicate_pair":
        rows[1][1] = "toy_A"
    elif mutation == "unknown_time":
        rows[0][2] = "T3"
    elif mutation == "nonfinite_wos":
        rows[0][3] = "NaN"
    elif mutation == "missing_patient":
        rows[0][1] = ""
    elif mutation == "patient_count":
        cfg["expected_patients"] = 99
    elif mutation == "row_count":
        rows.pop()
    elif mutation == "duplicate_header":
        header[-1] = "sample_id"
    with pytest.raises(ValueError, match=message):
        s.metadata(write_csv(tmp_path / "meta.csv", [header] + rows), cfg)


@pytest.mark.parametrize("data,features,message", [
    ([["", "b", "a"], ["f", 1, 2]], 1, "sample order"),
    ([["", "a", "a"], ["f", 1, 2]], 1, "Duplicate matrix sample"),
    ([["", "a", "b"], ["f", 1, 2], ["f", 3, 4]], 2, "Duplicate feature"),
    ([["", "a", "b"], ["", 1, 2]], 1, "feature ID"),
    ([["", "a", "b"], ["f", "nan", 2]], 1, "Nonfinite"),
    ([["", "a", "b"], ["f", "inf", 2]], 1, "Nonfinite"),
    ([["", "a", "b"], ["f", "", 2]], 1, "Non-numeric"),
    ([["", "a", "b"], ["f", "oops", 2]], 1, "Non-numeric"),
    ([["", "a", "b"], ["f", 1]], 1, "row width"),
    ([["", "a", "b"], ["f", 1, 2]], 2, "row count"),
    ([["", "f"], ["a", 1], ["b", 2]], 2, "sample order"),
])
def test_bad_matrix_fails(tmp_path, data, features, message):
    with pytest.raises(ValueError, match=message):
        s.matrix(write_csv(tmp_path / "matrix.csv", data), ["a", "b"], features)


def test_valid_matrix(tmp_path):
    path = write_csv(tmp_path / "matrix.csv", [["", "a", "b"], ["f", 0.5, 2]])
    assert s.matrix(path, ["a", "b"], 1)["noninteger_values"] == 1


def test_corrupt_checksum_rejected(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("bad")
    cfg = {"cache_dir": ".", "files": [{"name": p.name, "sha256": "0" * 64}]}
    with pytest.raises(ValueError, match="Checksum mismatch"):
        s.verified_files(tmp_path, cfg)


def test_source_and_semantics_are_fail_closed(tmp_path):
    (tmp_path / "config").mkdir()
    config = yaml.safe_load((ROOT / "config/datasets.yaml").read_text())
    config["covid"]["raw_counts"] = True
    (tmp_path / "config/datasets.yaml").write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="NOT raw counts"):
        s.covid_config(tmp_path)


def test_figure_registry_and_missing_prerequisites(tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "config/figures.yaml").write_bytes((ROOT / "config/figures.yaml").read_bytes())
    s.run_readiness(tmp_path, tmp_path / "readiness.json")
    assert s.registry(tmp_path)["main"] == ["Fig1", "Fig2", "Fig3", "Fig4", "Fig5"]
    with pytest.raises(ValueError, match="Missing scientific prerequisites"):
        s.require_figure(tmp_path, "Fig1")
    assert not (tmp_path / "results/figures/Fig1.pdf").exists()
    p = tmp_path / "results/scientific/Fig1.json"
    p.parent.mkdir(parents=True)
    p.write_text("{}")
    with pytest.raises(RuntimeError, match="NOT_IMPLEMENTED"):
        s.require_figure(tmp_path, "Fig1")


def test_future_gates_not_weakened():
    gates = s.load_config(ROOT, "simulations")["future_gates"]
    assert gates["rank_null"] == {"datasets": 200, "type_I_max": 0.10}
    assert gates["rank_one"] == {"datasets": 100, "exact_recovery_min": 0.80}
    assert gates["rank_two"] == {"datasets": 100, "exact_recovery_min": 0.80, "over_selection_max": 0.10}
