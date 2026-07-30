#!/usr/bin/env python3
"""Build the linked-image include for the snow and firn QA report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import string
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import quote


# The accepted first-pass scope is deliberately fixed even though inventory
# counts are discovered at runtime.
FAMILIES = ("promice", "esm_snowmip", "laugh_tests")
FIRN_FAMILIES = ("retmip", "imau", "sumup", "research_site")
FAMILY_TITLES = {
    "promice": "PROMICE",
    "esm_snowmip": "ESM-SnowMIP",
    "laugh_tests": "Laugh/Colbeck tests",
}
FIRN_FAMILY_TITLES = {
    "retmip": "RetMIP",
    "imau": "IMAU",
    "sumup": "SUMup",
    "research_site": "Research sites",
}
# These retired plot groups may remain in a case directory from an older run.
# They are outside the seasonal report contract even before a clean overwrite.
RETIRED_FIGURE_STEMS = frozenset(
    ("measurement_geometry", "measurement_geometry_2")
)
# The setup contract defines these eight channels as the complete native model-
# forcing interface.  Coverage below is evidence only; the manifest remains the
# sole authority for the advisory forcing-ready decision and complete windows.
PROMICE_REQUIRED_MET_CHANNELS = (
    "tair",
    "swd",
    "lwd",
    "albedo",
    "wspd",
    "rh",
    "psfc",
    "ppt",
)
# Coverage is reported only for model-development legs in the accepted seasonal
# scope. Evaluation payloads duplicate native observations and would otherwise
# inflate the tables without adding a distinct forcing/evaluation source.
COVERAGE_LEGS = {
    "promice": (
        ("promice", "met"),
        ("promice", "userdata"),
        ("mar3.11", "met"),
        ("mar3.11", "userdata"),
        ("merra2", "met"),
        ("merra2", "userdata"),
        ("racmo2.3p3", "met"),
        ("racmo2.3p3", "userdata"),
    ),
    "esm_snowmip": (("esm_snowmip", "met"),),
}
SOURCE_TITLES = {
    "promice": "PROMICE",
    "mar3.11": "MAR 3.11",
    "merra2": "MERRA-2",
    "racmo2.3p3": "RACMO 2.3p3",
    "esm_snowmip": "ESM-SnowMIP",
}
NATIVE_SUPPORT_FIGURE_STEMS = frozenset(
    ("radiation_fluxes", "quality_flags", "subsurface_temperature_string")
)
STRICT_DAILY_FIGURE_STEMS = frozenset(
    (
        "energy_balance",
        "turbulent_fluxes",
        "temperature_state",
        "mass_balance_terms",
        "mass_balance_components",
        "mar_combined_mass_diagnostics",
        "surface_height_depth",
        "meteorological_diagnostics",
    )
)


def load_qa(
    qa_path: Path,
    accepted_families: tuple[str, ...] = FAMILIES,
) -> tuple[dict, datetime, dict[str, set[str]]]:
    """Load passing QA and derive its exact nonempty case sets."""
    # Fail before figure discovery when the authoritative QA result is absent
    # or empty, because a gallery alone is not acceptance evidence.
    if not qa_path.is_file() or qa_path.stat().st_size == 0:
        raise ValueError(f"missing or empty artifact-QA JSON: {qa_path}")

    # Convert parse and schema failures into one concise command-line error.
    try:
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
        schema_version = qa["schema_version"]
        families = qa["families"]
        artifacts = qa["artifacts"]
        channels = qa["channels"]
        findings = qa["findings"]
        summary = qa["summary"]
        generated_at = datetime.fromisoformat(
            str(qa["generated_at"]).replace("Z", "+00:00")
        )
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid artifact-QA JSON: {qa_path}") from error

    # Validate the small schema surface consumed below so malformed reports
    # fail as QA errors instead of leaking Python attribute or key exceptions.
    required_counts = (
        "family_count",
        "case_count",
        "artifact_count",
        "channel_count",
        "finding_count",
        "placeholder_count",
        "error_count",
        "blocker_count",
        "warning_count",
    )
    if schema_version != "1.0":
        raise ValueError("artifact-QA schema_version must be 1.0")
    if not isinstance(families, list) or not isinstance(summary, dict):
        raise ValueError(f"invalid artifact-QA JSON: {qa_path}")
    if not isinstance(artifacts, list) or not isinstance(channels, list):
        raise ValueError(f"invalid artifact-QA JSON: {qa_path}")
    if not isinstance(findings, list):
        raise ValueError(f"invalid artifact-QA JSON: {qa_path}")
    if any(type(summary.get(key)) is not int or summary[key] < 0 for key in required_counts):
        raise ValueError(f"invalid artifact-QA JSON: {qa_path}")

    # Require the exact accepted family set; duplicate or unrelated families
    # would make the report look broader than the QA run actually was.
    if any(not isinstance(family, str) or not family for family in families):
        raise ValueError("artifact-QA family entries must be nonempty strings")
    if len(families) != len(accepted_families) or set(families) != set(
        accepted_families
    ):
        raise ValueError(
            f"artifact-QA families must be exactly {', '.join(accepted_families)}"
        )

    # Reconcile finding records with every severity summary so a stale or
    # hand-edited summary cannot authorize figures its records would reject.
    severities = ("error", "blocker", "warning", "placeholder")
    severity_counts = {severity: 0 for severity in severities}
    for finding in findings:
        if not isinstance(finding, dict):
            raise ValueError("artifact-QA findings contain an invalid severity")
        severity = finding.get("severity")
        if not isinstance(severity, str) or severity not in severity_counts:
            raise ValueError("artifact-QA findings contain an invalid severity")
        severity_counts[severity] += 1
    for severity, count in severity_counts.items():
        if summary[f"{severity}_count"] != count:
            raise ValueError("artifact-QA finding severity counts do not match summary")
    if summary["finding_count"] != len(findings):
        raise ValueError("artifact-QA finding count does not match its records")

    # Artifact case identifiers define the complete report inventory. Manifest
    # rows carry an empty case id and therefore do not create a figure case.
    expected_cases = {family: set() for family in accepted_families}
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ValueError("artifact-QA artifacts must be objects")
        family = artifact.get("dataset_family")
        case = artifact.get("case_id")
        if not isinstance(family, str) or family not in expected_cases:
            raise ValueError("artifact-QA artifacts contain an invalid family or case id")
        if not isinstance(case, str):
            raise ValueError("artifact-QA artifacts contain an invalid family or case id")
        if case:
            expected_cases[family].add(case)
    if any(not cases for cases in expected_cases.values()):
        raise ValueError("artifact QA must declare at least one nonempty case per family")
    if summary["artifact_count"] != len(artifacts):
        raise ValueError("artifact-QA artifact count does not match its records")
    if summary["channel_count"] != len(channels):
        raise ValueError("artifact-QA channel count does not match its records")
    if summary["case_count"] != sum(len(cases) for cases in expected_cases.values()):
        raise ValueError("artifact-QA case count does not match its artifact cases")

    # The first-class audit defines warnings as documented, nonblocking
    # diagnostics. Authorization still fails closed on its passed flag or any
    # error/blocker, while warning records remain reconciled and visible below.
    blocking = ("error_count", "blocker_count")
    if summary.get("passed") is not True or any(summary.get(key) != 0 for key in blocking):
        raise ValueError("artifact QA did not pass with zero errors and blockers")
    if summary.get("family_count") != len(accepted_families):
        raise ValueError("artifact-QA summary family count does not match its family set")
    if generated_at.tzinfo is None or generated_at.utcoffset() is None:
        raise ValueError("artifact-QA generated_at must include a timezone")

    return qa, generated_at, expected_cases


def load_promice_readiness(
    qa: dict,
    expected_cases: dict[str, set[str]],
) -> dict[str, dict]:
    """Load canonical native-PROMICE readiness metadata from its audited manifest."""
    # Bind the manifest to the audit's evaluation root and manifest artifact so
    # readiness cannot silently come from a different staging tree.
    evaluation_root = qa.get("evaluation_data_root")
    if not isinstance(evaluation_root, str) or not evaluation_root:
        raise ValueError("artifact QA must declare evaluation_data_root")
    manifest_path = (Path(evaluation_root) / "promice" / "manifest.json").resolve()
    audited_manifests = [
        artifact
        for artifact in qa["artifacts"]
        if isinstance(artifact, dict)
        and artifact.get("dataset_family") == "promice"
        and artifact.get("kind") == "manifest"
        and isinstance(artifact.get("path"), str)
        and Path(artifact.get("path", "")).resolve() == manifest_path
    ]
    if len(audited_manifests) != 1:
        raise ValueError("PROMICE readiness manifest is not the audited manifest")
    if not manifest_path.is_file() or manifest_path.stat().st_size == 0:
        raise ValueError(f"missing or empty PROMICE manifest: {manifest_path}")

    # Every artifact, including this manifest, passed the shared byte-identity
    # gate before readiness extraction. Read the verified manifest only to copy
    # its advisory decision; do not maintain a second hashing path here.
    try:
        manifest_bytes = manifest_path.read_bytes()
    except OSError as error:
        raise ValueError(f"invalid PROMICE manifest: {manifest_path}") from error
    try:
        manifest = json.loads(manifest_bytes)
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ValueError(f"invalid PROMICE manifest: {manifest_path}") from error
    if manifest.get("dataset_family") != "promice" or not isinstance(
        manifest.get("cases"), list
    ):
        raise ValueError(f"invalid PROMICE manifest: {manifest_path}")

    # Require exact case agreement and the canonical scalar/window fields.  The
    # report displays these values verbatim and never recomputes readiness.
    readiness: dict[str, dict] = {}
    for case in manifest["cases"]:
        if not isinstance(case, dict) or not isinstance(case.get("case_id"), str):
            raise ValueError(f"invalid PROMICE manifest: {manifest_path}")
        case_id = case["case_id"]
        if not case_id or case_id in readiness:
            raise ValueError(f"invalid PROMICE manifest: {manifest_path}")
        colocation = case.get("colocation")
        native = colocation.get("promice") if isinstance(colocation, dict) else None
        if not isinstance(native, dict):
            raise ValueError(f"PROMICE case {case_id} lacks native readiness metadata")
        ready = native.get("forcing_ready")
        reason = native.get("forcing_ready_reason")
        windows = native.get("forcing_complete_windows")
        if type(ready) is not bool or not isinstance(reason, str) or not isinstance(windows, list):
            raise ValueError(f"PROMICE case {case_id} has invalid readiness metadata")
        for window in windows:
            if not isinstance(window, dict):
                raise ValueError(f"PROMICE case {case_id} has invalid complete windows")
            if not all(isinstance(window.get(key), str) for key in ("start_time", "end_time")):
                raise ValueError(f"PROMICE case {case_id} has invalid complete windows")
            count = window.get("sample_count")
            if type(count) not in (int, float) or count < 0 or int(count) != count:
                raise ValueError(f"PROMICE case {case_id} has invalid complete windows")
        readiness[case_id] = {
            "forcing_ready": ready,
            "forcing_ready_reason": reason,
            "forcing_complete_windows": windows,
        }
    if set(readiness) != expected_cases["promice"]:
        raise ValueError("PROMICE manifest cases do not match artifact QA")
    return readiness


def record_count(record: dict, name: str) -> int:
    """Return one reconciled nonnegative integer count from a QA record."""
    # MATLAB JSON may encode integral counts as either JSON integers or floats.
    value = record.get(name)
    if type(value) not in (int, float) or value < 0 or int(value) != value:
        raise ValueError(f"artifact-QA record has invalid {name}")
    return int(value)


def file_identity(path: Path) -> tuple[int, str]:
    """Return one file's byte count and streaming SHA-256 digest."""
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as artifact_file:
        while chunk := artifact_file.read(8 * 1024 * 1024):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def verify_audited_artifacts(qa: dict) -> None:
    """Require every accepted artifact to match the bytes audited by QA."""
    # A passed audit must not contain a non-existing artifact. Validate every
    # identity before readiness or coverage extraction so no quantitative row
    # can outlive the file that supplied it.
    identities: dict[Path, tuple[int, str]] = {}
    for artifact in qa["artifacts"]:
        if artifact.get("exists") is not True:
            raise ValueError("accepted artifact QA contains a non-existing artifact")
        raw_path = artifact.get("path")
        if not isinstance(raw_path, str) or not raw_path:
            raise ValueError("accepted artifact QA contains an invalid artifact path")
        path = Path(raw_path).resolve()
        expected_sha256 = artifact.get("artifact_sha256")
        if (
            not isinstance(expected_sha256, str)
            or len(expected_sha256) != 64
            or any(character not in string.hexdigits for character in expected_sha256)
        ):
            raise ValueError(f"audited artifact lacks a valid artifact_sha256: {path}")
        try:
            expected_size = record_count(artifact, "artifact_size_bytes")
        except ValueError as error:
            raise ValueError(
                f"audited artifact has invalid artifact_size_bytes: {path}"
            ) from error
        if not path.is_file():
            raise ValueError(f"missing audited artifact: {path}")
        try:
            if path not in identities:
                identities[path] = file_identity(path)
            actual_size, actual_sha256 = identities[path]
        except OSError as error:
            raise ValueError(f"cannot read audited artifact: {path}") from error
        if expected_size != actual_size:
            raise ValueError(
                f"audited artifact artifact_size_bytes does not match disk: {path}"
            )
        if expected_sha256.lower() != actual_sha256:
            raise ValueError(
                f"audited artifact artifact_sha256 does not match disk: {path}"
            )


def coverage_evidence(
    qa: dict,
    expected_cases: dict[str, set[str]],
) -> dict[str, dict[str, list[dict]]]:
    """Collect reconciled per-variable QA coverage for the requested source legs."""
    # Index selected channel rows by their exact artifact identity. This keeps
    # extraction linear and prevents similarly named variables in another file
    # or table from being attributed to the wrong source leg.
    allowed = {family: set(legs) for family, legs in COVERAGE_LEGS.items()}
    channel_rows: dict[tuple[str, str, str, str, str], list[dict]] = defaultdict(list)
    for row in qa["channels"]:
        if not isinstance(row, dict):
            raise ValueError("artifact-QA channels must be objects")
        family = row.get("dataset_family")
        source = row.get("source")
        kind = row.get("kind")
        if family not in allowed or (source, kind) not in allowed[family]:
            continue
        case_id = row.get("case_id")
        path = row.get("path")
        if (
            not isinstance(case_id, str)
            or case_id not in expected_cases[family]
            or not isinstance(path, str)
            or not path
        ):
            raise ValueError("artifact-QA coverage channel has invalid case or path")
        channel_rows[(family, case_id, source, kind, path)].append(row)

    # Each selected source/kind leg must resolve to one artifact, because two
    # candidate files would make a source-level coverage row ambiguous.
    evidence = {
        family: {case_id: [] for case_id in expected_cases[family]}
        for family in COVERAGE_LEGS
    }
    seen_legs: set[tuple[str, str, str, str]] = set()
    for artifact in qa["artifacts"]:
        family = artifact.get("dataset_family")
        source = artifact.get("source")
        kind = artifact.get("kind")
        if family not in allowed or (source, kind) not in allowed[family]:
            continue
        case_id = artifact.get("case_id")
        path = artifact.get("path")
        if (
            not isinstance(case_id, str)
            or case_id not in expected_cases[family]
            or not isinstance(path, str)
            or not path
        ):
            raise ValueError("artifact-QA coverage artifact has invalid case or path")
        leg_key = (family, case_id, source, kind)
        if leg_key in seen_legs:
            raise ValueError(
                f"artifact-QA has duplicate {family}/{case_id}/{source}/{kind} artifacts"
            )
        seen_legs.add(leg_key)
        artifact_key = (*leg_key, path)
        rows = channel_rows.pop(artifact_key, [])
        if not rows:
            raise ValueError(
                f"artifact-QA {family}/{case_id}/{source}/{kind} has no channel records"
            )

        # Reconcile every displayed denominator and reject duplicate variable
        # names instead of silently double-counting one artifact table.
        channels: list[dict] = []
        names: set[str] = set()
        for row in rows:
            name = row.get("channel")
            if not isinstance(name, str) or not name or name in names:
                raise ValueError(
                    f"artifact-QA {family}/{case_id}/{source}/{kind} has invalid channels"
                )
            names.add(name)
            sample = record_count(row, "sample_count")
            finite = record_count(row, "finite_count")
            missing = record_count(row, "missing_count")
            missing_runs = record_count(row, "missing_run_count")
            longest_missing_run = record_count(row, "longest_missing_run_samples")
            if finite + missing != sample:
                raise ValueError(
                    f"artifact-QA {family}/{case_id}/{source}/{kind}/{name} counts do not reconcile"
                )
            if (
                (missing == 0 and (missing_runs != 0 or longest_missing_run != 0))
                or (
                    missing > 0
                    and not (
                        1 <= missing_runs <= missing
                        and 1 <= longest_missing_run <= missing
                    )
                )
            ):
                raise ValueError(
                    f"artifact-QA {family}/{case_id}/{source}/{kind}/{name} "
                    "missing-run counts do not reconcile"
                )
            unit = row.get("unit", "")
            if not isinstance(unit, str):
                raise ValueError(
                    f"artifact-QA {family}/{case_id}/{source}/{kind}/{name} has invalid unit"
                )
            channels.append(
                {
                    "name": name,
                    "unit": unit,
                    "sample": sample,
                    "finite": finite,
                    "missing": missing,
                    "missing_runs": missing_runs,
                    "longest_missing_run": longest_missing_run,
                }
            )
        channels.sort(key=lambda row: row["name"])
        evidence[family][case_id].append(
            {"source": source, "kind": kind, "artifact": artifact, "channels": channels}
        )

    # Selected channel rows without an exact artifact are not trustworthy
    # evidence and must not disappear merely because the report is artifact-led.
    if channel_rows:
        raise ValueError("artifact-QA coverage channels do not match audited artifacts")

    # Native PROMICE met/userdata and ESM-SnowMIP met are mandatory. RCM legs
    # remain optional per case and are reported whenever the manifest audit
    # includes them.
    for case_id, rows in evidence["promice"].items():
        pairs = {(row["source"], row["kind"]) for row in rows}
        if not {("promice", "met"), ("promice", "userdata")} <= pairs:
            raise ValueError(f"PROMICE case {case_id} lacks native coverage artifacts")
    for case_id, rows in evidence["esm_snowmip"].items():
        pairs = {(row["source"], row["kind"]) for row in rows}
        if pairs != {("esm_snowmip", "met")}:
            raise ValueError(f"ESM-SnowMIP case {case_id} must have one audited met artifact")

    # Apply the declared source/kind order so tables are stable across QA JSON
    # record order and regeneration.
    for family, cases in evidence.items():
        order = {leg: index for index, leg in enumerate(COVERAGE_LEGS[family])}
        for rows in cases.values():
            rows.sort(key=lambda row: order[(row["source"], row["kind"])])
    return evidence


def one_coverage_artifact(rows: list[dict], source: str, kind: str) -> dict:
    """Return one already-reconciled coverage artifact by source and kind."""
    matches = [row for row in rows if row["source"] == source and row["kind"] == kind]
    if len(matches) != 1:
        raise ValueError(f"coverage evidence must contain one {source} {kind} artifact")
    return matches[0]


def promice_evidence(
    qa: dict,
    readiness: dict[str, dict],
    coverage: dict[str, dict[str, list[dict]]],
) -> dict[str, dict]:
    """Pair manifest readiness with transparent native met/userdata QA counts."""
    evidence: dict[str, dict] = {}
    for case_id in sorted(readiness):
        coverage_rows = coverage["promice"][case_id]
        met = one_coverage_artifact(coverage_rows, "promice", "met")
        userdata = one_coverage_artifact(coverage_rows, "promice", "userdata")
        met_rows = met["channels"]
        by_name: dict[str, dict] = {}
        for name in PROMICE_REQUIRED_MET_CHANNELS:
            matches = [row for row in met_rows if row["name"] == name]
            if len(matches) != 1:
                raise ValueError(
                    f"PROMICE case {case_id} must have one audited met channel {name}"
                )
            row = matches[0]
            by_name[name] = {
                "sample": row["sample"],
                "finite": row["finite"],
                "missing": row["missing"],
            }

        # Native userdata is evaluation-support evidence, not a second forcing
        # readiness calculation.  Summarize every audited channel transparently.
        userdata_rows = userdata["channels"]
        userdata_sample = sum(row["sample"] for row in userdata_rows)
        userdata_finite = sum(row["finite"] for row in userdata_rows)
        userdata_missing = sum(row["missing"] for row in userdata_rows)

        placeholders = {
            finding.get("channel", "")
            for finding in qa["findings"]
            if finding.get("severity") == "placeholder"
            and finding.get("dataset_family") == "promice"
            and finding.get("case_id") == case_id
            and finding.get("source") == "promice"
            and finding.get("kind") == "met"
        }
        evidence[case_id] = {
            **readiness[case_id],
            "met_artifact": met["artifact"],
            "userdata_artifact": userdata["artifact"],
            "met_channels": by_name,
            "met_placeholders": placeholders,
            "coverage_artifacts": coverage_rows,
            "userdata_channel_count": len(userdata_rows),
            "userdata_finite_channel_count": sum(
                row["finite"] > 0 for row in userdata_rows
            ),
            "userdata_sample": userdata_sample,
            "userdata_finite": userdata_finite,
            "userdata_missing": userdata_missing,
        }
    return evidence


def collect_figures(
    figures_root: Path,
    generated_at: datetime,
    expected_cases: dict[str, set[str]],
    ignored_roots: frozenset[str] = frozenset(),
    retired_stems: frozenset[str] = RETIRED_FIGURE_STEMS,
) -> dict[str, list[Path]]:
    """Return deterministic current PNG inventories for the accepted families."""
    # plotVerificationArtifacts owns the canonical
    # <preview>/figures/<family>/<case> tree. Reject a missing or relocated
    # figures subroot instead of accepting the obsolete <preview>/<family> form.
    if not figures_root.is_dir():
        raise ValueError(f"missing figures root: {figures_root}")

    # Inspect every PNG under this scope so unrelated families cannot enter the
    # report merely by adding another directory. Explicitly ignored roots are
    # validated by their own independent QA call below.
    figures = {family: [] for family in expected_cases}
    extras: set[str] = set()
    zero_byte: list[str] = []
    stale: list[str] = []
    png_cases = {family: set() for family in expected_cases}
    case_problems: list[str] = []
    cutoff = generated_at.timestamp()

    # Compare immediate case directories before reading images so empty or
    # unexpected directories cannot disappear from the report inventory.
    for family in expected_cases:
        family_root = figures_root / family
        actual_dirs = (
            {entry.name for entry in family_root.iterdir() if entry.is_dir()}
            if family_root.is_dir()
            else set()
        )
        missing_dirs = sorted(expected_cases[family] - actual_dirs)
        unexpected_dirs = sorted(actual_dirs - expected_cases[family])
        if missing_dirs:
            case_problems.append(f"{family} missing case directories: {', '.join(missing_dirs)}")
        if unexpected_dirs:
            case_problems.append(
                f"{family} unexpected case directories: {', '.join(unexpected_dirs)}"
            )

    for image_path in sorted(figures_root.rglob("*"), key=lambda path: path.as_posix()):
        if not image_path.is_file() or image_path.suffix.lower() != ".png":
            continue
        relative = image_path.relative_to(figures_root)
        family = relative.parts[0]
        if family not in figures:
            if family not in ignored_roots:
                extras.add(family)
            continue
        if image_path.stem in retired_stems:
            continue
        case = relative.parts[1] if len(relative.parts) > 2 else ""
        png_cases[family].add(case)
        stat = image_path.stat()
        if stat.st_size == 0:
            zero_byte.append(relative.as_posix())
        if stat.st_mtime < cutoff:
            stale.append(relative.as_posix())
        figures[family].append(image_path)

    # Report all inventory defects together so one run gives a complete and
    # deterministic correction list without weakening any gate.
    problems: list[str] = []
    problems.extend(case_problems)
    for family in expected_cases:
        missing_cases = sorted(expected_cases[family] - png_cases[family])
        unexpected_cases = sorted(png_cases[family] - expected_cases[family])
        if missing_cases:
            problems.append(f"{family} missing PNG cases: {', '.join(missing_cases)}")
        if unexpected_cases:
            problems.append(f"{family} unexpected PNG cases: {', '.join(unexpected_cases)}")
    if extras:
        problems.append(f"extra PNG families: {', '.join(sorted(extras))}")
    if zero_byte:
        problems.append(f"zero-byte PNGs: {', '.join(zero_byte)}")
    if stale:
        problems.append(f"PNGs older than artifact QA: {', '.join(stale)}")
    if problems:
        raise ValueError("; ".join(problems))

    return figures


def validate_firn_figure_index(
    index_path: Path,
    figures_root: Path,
    figures: dict[str, list[Path]],
    expected_cases: dict[str, set[str]],
) -> None:
    """Require the MATLAB figure ledger to equal the nested firn PNG tree."""
    # The post-promotion driver round-trips this CSV before sealing completion;
    # the report independently binds its rows to the exact canonical PNG files.
    if not index_path.is_file() or index_path.stat().st_size == 0:
        raise ValueError(f"missing or empty firn figure index: {index_path}")
    with index_path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"dataset_family", "case_id", "figure_file"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("firn figure index lacks required columns")
        rows = list(reader)
    if not rows:
        raise ValueError("firn figure index contains no rows")

    indexed: list[Path] = []
    indexed_cases: dict[str, set[str]] = {family: set() for family in expected_cases}
    root = figures_root.resolve()
    for row in rows:
        family = row["dataset_family"]
        case = row["case_id"]
        filename = Path(row["figure_file"])
        if (
            family not in expected_cases
            or case not in expected_cases[family]
            or not filename.is_absolute()
        ):
            raise ValueError("firn figure index contains an invalid family, case, or path")
        resolved = filename.resolve()
        if not resolved.is_relative_to(root):
            raise ValueError("firn figure index contains a path outside figures/firn")
        relative = resolved.relative_to(root)
        if len(relative.parts) < 3 or relative.parts[:2] != (family, case):
            raise ValueError("firn figure index path disagrees with its family or case")
        indexed.append(resolved)
        indexed_cases[family].add(case)

    tree = sorted(path.resolve() for paths in figures.values() for path in paths)
    if len(indexed) != len(set(indexed)) or sorted(indexed) != tree:
        raise ValueError("firn figure index does not exactly match the firn PNG tree")
    if indexed_cases != expected_cases:
        raise ValueError("firn figure index case map does not match firn artifact QA")


def display_name(value: str) -> str:
    """Convert a stable path token into a compact display label."""
    # Replace control characters before the label reaches Markdown, HTML, or
    # browser accessibility text.
    printable = "".join(character if character.isprintable() else " " for character in value)
    return printable.replace("_", " ").replace("-", " ").title()


def markdown_heading(value: str) -> str:
    """Escape a case token for literal use in a Markdown heading."""
    # Escape all ASCII punctuation so Pandoc cannot reinterpret a case token as
    # links, emphasis, math, superscript, strikeout, or raw markup.
    return "".join(
        f"\\{character}" if character in string.punctuation else character
        for character in display_name(value)
    )


def markdown_cell(value: object) -> str:
    """Escape one compact scalar for literal use in a Markdown table cell."""
    # Keep source-provided reasons on one line and prevent their punctuation
    # from adding table columns.
    return str(value).replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def coverage_anchor(family: str, case: str) -> str:
    """Return one safe stable fragment identifier for a coverage table."""
    # A short digest keeps arbitrary source case tokens out of raw HTML IDs.
    token = hashlib.sha256(f"{family}/{case}".encode("utf-8")).hexdigest()[:12]
    return f"coverage-{family}-{token}"


def coverage_text(finite: int, sample: int) -> str:
    """Format a finite-value percentage without making a readiness claim."""
    # An empty audited channel has no meaningful denominator.
    return "—" if sample == 0 else f"{100 * finite / sample:.2f}%"


def cadence_text(artifact: dict) -> str:
    """Format the audited artifact cadence in the smallest useful unit."""
    seconds = artifact.get("cadence_seconds")
    if type(seconds) not in (int, float) or seconds <= 0:
        return "—"
    if seconds % 3600 == 0:
        return f"{seconds / 3600:g} h"
    if seconds % 60 == 0:
        return f"{seconds / 60:g} min"
    return f"{seconds:g} s"


def complete_windows_text(windows: list[dict]) -> str:
    """Render canonical manifest windows without deriving or filtering them."""
    if not windows:
        return "None recorded"
    return "<br>".join(
        f"{markdown_cell(window['start_time'])} to "
        f"{markdown_cell(window['end_time'])} "
        f"({record_count(window, 'sample_count'):,} samples)"
        for window in windows
    )


def coverage_case_lines(title: str, artifacts: list[dict]) -> list[str]:
    """Render one shared per-source, per-variable artifact-QA coverage table."""
    # Keep absent-value volume, contiguous-run count, and maximum run length in
    # separate columns; none is a substitute for either of the others.
    lines = [
        f"#### {title}",
        "",
        (
            "| source | artifact kind | status | cadence | variable | unit | "
            "finite / audited samples | finite coverage | missing samples | "
            "contiguous missing runs | longest missing run (samples) |"
        ),
        "|---|---|---|---:|---|---|---:|---:|---:|---:|---:|",
    ]
    for entry in artifacts:
        artifact = entry["artifact"]
        source = SOURCE_TITLES.get(entry["source"], display_name(entry["source"]))
        for row in entry["channels"]:
            lines.append(
                "| "
                + " | ".join(
                    (
                        markdown_cell(source),
                        markdown_cell(entry["kind"]),
                        markdown_cell(artifact.get("status", "")),
                        cadence_text(artifact),
                        f"`{markdown_cell(row['name'])}`",
                        markdown_cell(row["unit"] or "—"),
                        f"{row['finite']:,} / {row['sample']:,}",
                        coverage_text(row["finite"], row["sample"]),
                        f"{row['missing']:,}",
                        f"{row['missing_runs']:,}",
                        f"{row['longest_missing_run']:,}",
                    )
                )
                + " |"
            )
    return lines + [""]


def promice_summary_lines(evidence: dict[str, dict]) -> list[str]:
    """Build the report-wide scan table for native PROMICE readiness."""
    ready_count = sum(case["forcing_ready"] for case in evidence.values())
    window_count = sum(
        bool(case["forcing_complete_windows"]) for case in evidence.values()
    )
    lines = [
        "## PROMICE native forcing readiness",
        "",
        (
            "The **forcing ready** value and complete windows below are copied "
            "from each case manifest. Coverage percentages come from the first-class "
            "artifact audit and do not override that advisory manifest decision."
        ),
        "",
        (
            f"**Ready sites: {ready_count} of {len(evidence)}. Sites with at least one "
            f"recorded complete window: {window_count} of {len(evidence)}.**"
        ),
        "",
        (
            "| case | forcing ready | complete windows | met status | minimum required "
            "coverage | required channels with missing samples/placeholders | manifest reason |"
        ),
        "|---|---|---|---|---:|---|---|",
    ]
    for case_id, case in evidence.items():
        channels = case["met_channels"]
        percentages = [
            100 * row["finite"] / row["sample"] if row["sample"] else 0
            for row in channels.values()
        ]
        missing_channels = []
        for name, row in channels.items():
            if row["missing"] or name in case["met_placeholders"]:
                suffix = " placeholder" if name in case["met_placeholders"] else ""
                missing_channels.append(f"{name}: {row['missing']:,}{suffix}")
        lines.append(
            "| "
            + " | ".join(
                (
                    markdown_cell(case_id),
                    "yes" if case["forcing_ready"] else "no",
                    complete_windows_text(case["forcing_complete_windows"]),
                    markdown_cell(case["met_artifact"].get("status", "")),
                    f"{min(percentages):.2f}%",
                    markdown_cell(", ".join(missing_channels) if missing_channels else "none"),
                    markdown_cell(case["forcing_ready_reason"] or "—"),
                )
            )
            + " |"
        )
    return lines + [""]


def promice_case_lines(case_id: str, case: dict) -> list[str]:
    """Build one case's manifest decision and link to source-leg coverage."""
    met = case["met_artifact"]
    userdata = case["userdata_artifact"]
    lines = [
        "#### Quantitative readiness evidence",
        "",
        (
            f"Manifest forcing ready: **{'yes' if case['forcing_ready'] else 'no'}**. "
            f"Complete windows: {complete_windows_text(case['forcing_complete_windows'])}. "
            f"Reason: {markdown_cell(case['forcing_ready_reason'] or '—')}."
        ),
        "",
        (
            f"Native PROMICE met artifact: status `{markdown_cell(met.get('status', ''))}`, "
            f"cadence {cadence_text(met)}, span "
            f"`{markdown_cell(met.get('time_start', ''))}` to "
            f"`{markdown_cell(met.get('time_end', ''))}`."
        ),
        "",
        (
            f'<p><a href="#{coverage_anchor("promice", case_id)}">Open audited '
            "per-source channel coverage in the final appendix.</a></p>"
        ),
        "",
    ]
    lines.extend(
        [
            (
                "Native PROMICE userdata is evaluation-support evidence, not a forcing-"
                "readiness decision. "
                f"Status `{markdown_cell(userdata.get('status', ''))}`, cadence "
                f"{cadence_text(userdata)}, span "
                f"`{markdown_cell(userdata.get('time_start', ''))}` to "
                f"`{markdown_cell(userdata.get('time_end', ''))}`; "
                f"{case['userdata_finite_channel_count']:,} of "
                f"{case['userdata_channel_count']:,} audited channels contain finite data; "
                f"{case['userdata_finite']:,} of {case['userdata_sample']:,} audited "
                f"channel values are finite "
                f"({coverage_text(case['userdata_finite'], case['userdata_sample'])}), "
                f"with {case['userdata_missing']:,} missing samples."
            ),
            "",
        ]
    )
    return lines


def figure_interpretation(stem: str) -> str:
    """Label readiness views separately from completeness-gated diagnostics."""
    # Split suffixes are renderer inventory details, not different aggregation
    # policies, so resolve them back to the canonical group stem.
    canonical = stem.removesuffix("_2")
    if canonical == "met_forcing":
        return (
            "Readiness view: staged model-met support. Use with the manifest "
            "decision and linked required-variable coverage table in the final appendix."
        )
    if canonical in NATIVE_SUPPORT_FIGURE_STEMS:
        return (
            "Native-support view: actual source timestamps with explicit gap breaks; "
            "use this to inspect evaluation/source support."
        )
    if canonical == "energy_balance":
        return (
            "Daily diagnostic: dense subdaily series require complete native UTC days. "
            "When MAR userdata lacks a radiation balance, the renderer derives only "
            "the missing display term in memory: swn = swd * (1 - albedo), lwn uses "
            "icemodel.surface.net_longwave_radiation(tsfc, lwd) with the same surface "
            "emissivity and Kelvin temperature convention as icemodel.processmet, and "
            "netr = swn + lwn. Existing native terms and the staged artifact remain "
            "unchanged."
        )
    if canonical == "surface_albedo":
        return (
            "Daily diagnostic: radiometer observations require a complete native "
            "UTC timestamp grid, source solar screening, and at least six hours of "
            "valid support before computing the ratio of reflected to incoming daily "
            "shortwave energy. Native MAR/MERRA-2 modeled surface-state albedo uses an "
            "exact-grid arithmetic finite-sample daily mean and retains valid "
            "polar-night coverage. RACMO/model-ratio albedo uses positive SWD as the "
            "energy weight without the radiometer six-hour gate; an all-zero-SWD "
            "polar-night day remains undefined. Native daily products without "
            "collocated SWD retain their finite-sample values."
        )
    if canonical in STRICT_DAILY_FIGURE_STEMS:
        return (
            "Daily diagnostic: dense subdaily series require complete native UTC days. "
            "A missing daily point does not by itself prove a source or forcing gap."
        )
    return ""


def build_report(qa_path: Path, preview_root: Path, output_path: Path) -> int:
    """Validate inputs and write a deterministic linked-image Markdown include."""
    # Resolve all paths once and keep the authoritative QA and generated output
    # inside their dedicated preview subtrees.
    preview_root = preview_root.resolve()
    qa_path = qa_path.resolve()
    output_path = output_path.resolve()
    if not qa_path.is_relative_to(preview_root / "qa"):
        raise ValueError("artifact-QA JSON must be beneath <preview-root>/qa")
    if not output_path.is_relative_to(preview_root / "report"):
        raise ValueError("generated include must be beneath <preview-root>/report")

    qa, generated_at, expected_cases = load_qa(qa_path)
    verify_audited_artifacts(qa)
    readiness = load_promice_readiness(qa, expected_cases)
    coverage = coverage_evidence(qa, expected_cases)
    evidence = promice_evidence(qa, readiness, coverage)
    figures_root = preview_root / "figures"
    figures = collect_figures(
        figures_root,
        generated_at,
        expected_cases,
        ignored_roots=frozenset(("firn",)),
    )

    # The canonical firn proof owns a separate exact QA scope and nested figure
    # root. Reuse the same schema and inventory checks before adding its appendix.
    firn_qa_path = preview_root / "qa/firn/artifact_qa.json"
    firn_qa, firn_generated_at, firn_expected_cases = load_qa(
        firn_qa_path, FIRN_FAMILIES
    )
    verify_audited_artifacts(firn_qa)
    firn_figures_root = figures_root / "firn"
    firn_figures = collect_figures(
        firn_figures_root,
        firn_generated_at,
        firn_expected_cases,
        retired_stems=frozenset(),
    )
    firn_index_path = preview_root / "qa/firn_figure_index.csv"
    validate_firn_figure_index(
        firn_index_path, firn_figures_root, firn_figures, firn_expected_cases
    )
    if firn_index_path.stat().st_mtime_ns < firn_qa_path.stat().st_mtime_ns:
        raise ValueError("firn figure index is older than firn artifact-QA JSON")

    # Link both audit representations from the rendered output directory. The
    # concise Markdown is for review; JSON remains the quantitative authority.
    qa_summary_path = qa_path.with_suffix(".md")
    if not qa_summary_path.is_file() or qa_summary_path.stat().st_size == 0:
        raise ValueError(f"missing or empty artifact-QA Markdown: {qa_summary_path}")
    if qa_summary_path.stat().st_mtime_ns < qa_path.stat().st_mtime_ns:
        raise ValueError("artifact-QA Markdown is older than artifact-QA JSON")
    firn_qa_summary_path = firn_qa_path.with_suffix(".md")
    if not firn_qa_summary_path.is_file() or firn_qa_summary_path.stat().st_size == 0:
        raise ValueError(
            f"missing or empty firn artifact-QA Markdown: {firn_qa_summary_path}"
        )
    if firn_qa_summary_path.stat().st_mtime_ns < firn_qa_path.stat().st_mtime_ns:
        raise ValueError("firn artifact-QA Markdown is older than its JSON")
    qa_relative = Path(os.path.relpath(qa_path, output_path.parent)).as_posix()
    qa_href = quote(qa_relative, safe="/._-")
    qa_summary_relative = Path(
        os.path.relpath(qa_summary_path, output_path.parent)
    ).as_posix()
    qa_summary_href = quote(qa_summary_relative, safe="/._-")
    firn_qa_relative = Path(
        os.path.relpath(firn_qa_path, output_path.parent)
    ).as_posix()
    firn_qa_href = quote(firn_qa_relative, safe="/._-")
    firn_qa_summary_relative = Path(
        os.path.relpath(firn_qa_summary_path, output_path.parent)
    ).as_posix()
    firn_qa_summary_href = quote(firn_qa_summary_relative, safe="/._-")
    firn_index_relative = Path(
        os.path.relpath(firn_index_path, output_path.parent)
    ).as_posix()
    firn_index_href = quote(firn_index_relative, safe="/._-")
    summary = qa["summary"]
    lines = [
        "<!-- Generated by icemodel.verification.report; do not edit. -->",
        "",
        "## Artifact QA status",
        "",
        (
            f"Artifact QA passed at `{qa['generated_at']}` for "
            f"{summary['case_count']} cases and {summary['artifact_count']} artifacts. "
            f"It recorded {summary['placeholder_count']} informational "
            f"{'placeholder' if summary['placeholder_count'] == 1 else 'placeholders'}, "
            f"{summary['warning_count']} documented "
            f"{'warning' if summary['warning_count'] == 1 else 'warnings'}, and "
            "zero errors or blockers."
        ),
        "",
        (
            f'<p><a href="{qa_summary_href}">Open the human-readable artifact-QA '
            f'summary</a> or <a href="{qa_href}">the machine-readable artifact-QA '
            "result</a>.</p>"
        ),
        "",
        "## How to read the figures",
        "",
        (
            "Use `met_forcing.png` and the native-support figures to judge staged "
            "forcing/evaluation support, together with the quantitative audit tables. "
            "The manifest's **forcing ready** field is the authoritative advisory "
            "decision; coverage percentages do not replace it."
        ),
        "",
        (
            "Figures labeled **Daily diagnostic** are visualization reductions. Dense "
            "subdaily series appear only for complete native UTC days, so a blank day in "
            "those panels is not by itself evidence that the staged met or userdata is "
            "missing. `surface_albedo.png` is source-aware: radiometer observations "
            "require the complete timestamp grid, at least six hours of valid support, "
            "and source solar-elevation screening before weighting finite albedo by "
            "positive downwelling shortwave. That result is the ratio of reflected to "
            "incoming daily shortwave energy. Native MAR/MERRA-2 modeled surface-state "
            "albedo retains valid polar-night coverage through an exact-grid arithmetic "
            "finite-sample daily mean. RACMO/model-ratio albedo instead uses positive "
            "SWD as the energy weight without the radiometer six-hour gate; an "
            "all-zero-SWD polar-night day remains undefined. Native daily products "
            "without collocated SWD retain their finite-sample values."
        ),
        "",
        (
            "Coverage tables keep three distinct audit metrics: `missing_count` is the "
            "number of individual missing samples, `missing_run_count` is the number of "
            "contiguous missing runs, and `longest_missing_run_samples` is the maximum "
            "run length in samples."
        ),
        "",
    ]
    lines.extend(promice_summary_lines(evidence))

    # Group each family by its first directory component so the Quarto table of
    # contents provides direct navigation without a second inventory format.
    total = 0
    for family in FAMILIES:
        family_root = figures_root / family
        grouped: dict[str, list[Path]] = defaultdict(list)
        for image_path in figures[family]:
            relative = image_path.relative_to(family_root)
            case = relative.parts[0]
            grouped[case].append(image_path)

        lines.extend(
            [
                f"## {FAMILY_TITLES[family]}",
                "",
                (
                    f"{len(figures[family])} "
                    f"{'figure' if len(figures[family]) == 1 else 'figures'} across "
                    f"{len(grouped)} {'case' if len(grouped) == 1 else 'cases'}."
                ),
                "",
            ]
        )
        for case in sorted(grouped):
            lines.extend([f"### {markdown_heading(case)}", ""])
            if family == "promice":
                lines.extend(promice_case_lines(case, evidence[case]))
            elif family == "esm_snowmip":
                lines.extend(
                    [
                        (
                            f'<p><a href="#{coverage_anchor(family, case)}">Open '
                            "audited staged-met channel coverage in the final "
                            "appendix.</a></p>"
                        ),
                        "",
                    ]
                )
            for image_path in grouped[case]:
                relative = Path(os.path.relpath(image_path, output_path.parent)).as_posix()
                href = quote(relative, safe="/._-")
                label = image_path.relative_to(preview_root).as_posix()
                alt = html.escape(display_name(image_path.stem), quote=True)
                caption = html.escape(label)
                interpretation = figure_interpretation(image_path.stem)
                interpretation_html = (
                    f'<br><span class="figure-interpretation">{html.escape(interpretation)}</span>'
                    if interpretation
                    else ""
                )
                lines.extend(
                    [
                        "<figure>",
                        f'  <a href="{href}"><img src="{href}" loading="lazy" alt="{alt}"></a>',
                        f"  <figcaption><code>{caption}</code>{interpretation_html}</figcaption>",
                        "</figure>",
                        "",
                    ]
                )
                total += 1

    # Append the canonical firn-family set without introducing another report
    # generator. Its compact QA and MATLAB figure ledger remain direct links.
    firn_summary = firn_qa["summary"]
    lines.extend(
        [
            "## Firn artifact QA status",
            "",
            (
                f"Firn artifact QA passed at `{firn_qa['generated_at']}` for "
                f"{firn_summary['case_count']} cases and "
                f"{firn_summary['artifact_count']} artifacts. It recorded "
                f"{firn_summary['placeholder_count']} informational placeholders, "
                f"{firn_summary['warning_count']} documented warnings, and zero "
                "errors or blockers."
            ),
            "",
            (
                f'<p><a href="{firn_qa_summary_href}">Open the firn artifact-QA '
                f'summary</a>, <a href="{firn_qa_href}">the machine-readable '
                f'result</a>, or <a href="{firn_index_href}">the exact figure '
                "ledger</a>.</p>"
            ),
            "",
        ]
    )
    for family in FIRN_FAMILIES:
        family_root = firn_figures_root / family
        grouped: dict[str, list[Path]] = defaultdict(list)
        for image_path in firn_figures[family]:
            case = image_path.relative_to(family_root).parts[0]
            grouped[case].append(image_path)
        lines.extend(
            [
                f"## {FIRN_FAMILY_TITLES[family]}",
                "",
                (
                    f"{len(firn_figures[family])} "
                    f"{'figure' if len(firn_figures[family]) == 1 else 'figures'} "
                    f"across {len(grouped)} "
                    f"{'case' if len(grouped) == 1 else 'cases'}."
                ),
                "",
            ]
        )
        for case in sorted(grouped):
            lines.extend([f"### {markdown_heading(case)}", ""])
            for image_path in grouped[case]:
                relative = Path(
                    os.path.relpath(image_path, output_path.parent)
                ).as_posix()
                href = quote(relative, safe="/._-")
                label = image_path.relative_to(preview_root).as_posix()
                alt = html.escape(display_name(image_path.stem), quote=True)
                interpretation = figure_interpretation(image_path.stem)
                interpretation_html = (
                    f'<br><span class="figure-interpretation">{html.escape(interpretation)}</span>'
                    if interpretation
                    else ""
                )
                lines.extend(
                    [
                        "<figure>",
                        f'  <a href="{href}"><img src="{href}" loading="lazy" alt="{alt}"></a>',
                        f"  <figcaption><code>{html.escape(label)}</code>{interpretation_html}</figcaption>",
                        "</figure>",
                        "",
                    ]
                )
                total += 1

    # Keep large quantitative tables after every figure so the case-by-case
    # visual scan remains uninterrupted while retaining direct case links.
    lines.extend(
        [
            "## Audited per-source channel coverage appendix",
            "",
            (
                "Detailed channel tables are collected here to keep the figure "
                "sections compact. Each contributing case links directly to its table."
            ),
            "",
        ]
    )
    for family in COVERAGE_LEGS:
        for case in sorted(coverage[family]):
            lines.extend(
                [
                    f'<span id="{coverage_anchor(family, case)}"></span>',
                    "",
                ]
            )
            lines.extend(
                coverage_case_lines(
                    f"{FAMILY_TITLES[family]} — {markdown_heading(case)}",
                    coverage[family][case],
                )
            )

    # Avoid touching the generated include when its deterministic content is
    # unchanged, which makes no-op regeneration visible through file mtimes.
    content = "\n".join(lines).rstrip() + "\n"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not output_path.exists() or output_path.read_text(encoding="utf-8") != content:
        output_path.write_text(content, encoding="utf-8")

    return total


def main(argv: list[str] | None = None) -> int:
    """Parse paths, build the include, and print a compact inventory summary."""
    # Defaults are anchored to the tracked script rather than the caller's
    # working directory so the documented command works from any directory.
    repo_root = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa", type=Path, default=repo_root / "data/preview/qa/artifact_qa.json")
    parser.add_argument("--preview-root", type=Path, default=repo_root / "data/preview")
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root / "data/preview/report/figures.generated.md",
    )
    args = parser.parse_args(argv)

    # Resolve caller-supplied relative paths before computing report links.
    qa_path = args.qa.resolve()
    preview_root = args.preview_root.resolve()
    output_path = args.output.resolve()
    total = build_report(qa_path, preview_root, output_path)
    print(f"Generated {output_path} with {total} linked figures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
