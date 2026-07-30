"""Focused standard-library tests for the combined snow/firn report generator."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import shutil
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path


# Load the report implementation from the MATLAB verification namespace without
# treating its plus-prefixed package directories as Python packages.
REPO_ROOT = Path(__file__).resolve().parents[2]
REPORT_FILE = (
    REPO_ROOT
    / "icemodel/+icemodel/+verification/+report/build_snow_artifact_qa.py"
)
SPEC = importlib.util.spec_from_file_location("build_snow_artifact_qa", REPORT_FILE)
assert SPEC is not None and SPEC.loader is not None
report_builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report_builder)


class BuildSnowArtifactQaTest(unittest.TestCase):
    """Prove the report gate and relative-link output with tiny fixtures."""

    def make_fixture(
        self,
        cases: dict[str, tuple[str, ...]] | None = None,
    ) -> tuple[Path, Path, Path, float]:
        """Create current case PNGs and a reconciled passing QA JSON."""
        # Keep each fixture independent so failure mutations cannot leak between
        # subtests that exercise different validation branches.
        root = Path(tempfile.mkdtemp(dir=self.temp_root.name))
        preview = root / "data/preview"
        qa_path = preview / "qa/artifact_qa.json"
        output = preview / "report/figures.generated.md"
        qa_time = 1_700_000_000.0
        generated_at = datetime.fromtimestamp(qa_time, timezone.utc).isoformat().replace(
            "+00:00", "Z"
        )
        case_map = cases or {family: ("case_a",) for family in report_builder.FAMILIES}
        evaluation_root = root / "data/eval"
        input_root = root / "data/input"
        manifest_path = evaluation_root / "promice/manifest.json"

        # Keep one false/no-window native case and make any additional PROMICE
        # cases true with an explicit window, covering both canonical states.
        manifest_cases = []
        artifacts = [
            {
                "dataset_family": "promice",
                "case_id": "",
                "source": "",
                "kind": "manifest",
                "path": str(manifest_path),
            }
        ]
        channels = []
        findings = []
        for index, case in enumerate(case_map["promice"]):
            ready = index > 0
            windows = (
                [
                    {
                        "start_time": "2020-01-01 00:00:00",
                        "end_time": "2020-01-01 23:45:00",
                        "sample_count": 96,
                    }
                ]
                if ready
                else []
            )
            reason = "" if ready else "required met placeholder/gap channel(s): ppt"
            manifest_cases.append(
                {
                    "case_id": case,
                    "colocation": {
                        "promice": {
                            "forcing_ready": ready,
                            "forcing_ready_reason": reason,
                            "forcing_complete_windows": windows,
                        }
                    },
                }
            )
            native_paths = {}
            for kind, cadence, samples, status in (
                ("met", 900, 96, "passed" if ready else "placeholder"),
                ("userdata", 3600, 24, "warning"),
            ):
                artifact_path = input_root / kind / f"{case}.mat"
                native_paths[kind] = artifact_path
                artifacts.append(
                    {
                        "dataset_family": "promice",
                        "case_id": case,
                        "source": "promice",
                        "kind": kind,
                        "path": str(artifact_path),
                        "status": status,
                        "time_start": "2020-01-01 00:00:00",
                        "time_end": "2020-01-01 23:45:00",
                        "cadence_seconds": cadence,
                        "n_samples": samples,
                        "n_channels": len(report_builder.PROMICE_REQUIRED_MET_CHANNELS),
                    }
                )
            for name in report_builder.PROMICE_REQUIRED_MET_CHANNELS:
                missing = 96 if name == "ppt" and not ready else 0
                channels.append(
                    {
                        "dataset_family": "promice",
                        "case_id": case,
                        "source": "promice",
                        "kind": "met",
                        "path": str(native_paths["met"]),
                        "channel": name,
                        "sample_count": 96,
                        "finite_count": 96 - missing,
                        "missing_count": missing,
                        "missing_run_count": 1 if missing else 0,
                        "longest_missing_run_samples": missing,
                    }
                )
            channels.append(
                {
                    "dataset_family": "promice",
                    "case_id": case,
                    "source": "promice",
                    "kind": "userdata",
                    "path": str(native_paths["userdata"]),
                    "channel": "swd",
                    "sample_count": 24,
                    "finite_count": 20,
                    "missing_count": 4,
                    "missing_run_count": 2,
                    "longest_missing_run_samples": 3,
                }
            )

            # Include every requested PROMICE model-development source leg so
            # the shared coverage table is exercised on optional RCM artifacts
            # as well as native met/userdata.
            for source, kind, cadence, samples, channel, missing in (
                ("mar3.11", "met", 900, 96, "swd", 0),
                ("mar3.11", "userdata", 3600, 24, "swn", 1),
                ("merra2", "userdata", 3600, 24, "swn", 0),
                ("racmo2.3p3", "userdata", 3600, 24, "swn", 2),
            ):
                artifact_path = input_root / kind / source / f"{case}.mat"
                artifacts.append(
                    {
                        "dataset_family": "promice",
                        "case_id": case,
                        "source": source,
                        "kind": kind,
                        "path": str(artifact_path),
                        "status": "passed",
                        "time_start": "2020-01-01 00:00:00",
                        "time_end": "2020-01-01 23:45:00",
                        "cadence_seconds": cadence,
                        "n_samples": samples,
                        "n_channels": 1,
                    }
                )
                channels.append(
                    {
                        "dataset_family": "promice",
                        "case_id": case,
                        "source": source,
                        "kind": kind,
                        "path": str(artifact_path),
                        "channel": channel,
                        "unit": "W m-2",
                        "sample_count": samples,
                        "finite_count": samples - missing,
                        "missing_count": missing,
                        "missing_run_count": 1 if missing else 0,
                        "longest_missing_run_samples": missing,
                    }
                )
            if not ready:
                findings.append(
                    {
                        "severity": "placeholder",
                        "dataset_family": "promice",
                        "case_id": case,
                        "source": "promice",
                        "kind": "met",
                        "channel": "ppt",
                    }
                )
        manifest_path.parent.mkdir(parents=True)
        manifest_bytes = json.dumps(
            {"dataset_family": "promice", "cases": manifest_cases}
        ).encode("utf-8")
        manifest_path.write_bytes(manifest_bytes)
        # ESM-SnowMIP staged met is model-forcing evidence and therefore gets
        # the same per-variable channel treatment as the PROMICE source legs.
        for case in case_map["esm_snowmip"]:
            artifact_path = input_root / "met/esm_snowmip" / f"{case}.mat"
            artifacts.append(
                {
                    "dataset_family": "esm_snowmip",
                    "case_id": case,
                    "source": "esm_snowmip",
                    "kind": "met",
                    "path": str(artifact_path),
                    "status": "passed",
                    "time_start": "2020-01-01 00:00:00",
                    "time_end": "2020-01-01 23:45:00",
                    "cadence_seconds": 900,
                    "n_samples": 96,
                    "n_channels": 3,
                }
            )
            for name, missing, unit in (
                ("swd", 0, "W m-2"),
                ("lwd", 3, "W m-2"),
                ("tair", 0, "degC"),
            ):
                channels.append(
                    {
                        "dataset_family": "esm_snowmip",
                        "case_id": case,
                        "source": "esm_snowmip",
                        "kind": "met",
                        "path": str(artifact_path),
                        "channel": name,
                        "unit": unit,
                        "sample_count": 96,
                        "finite_count": 96 - missing,
                        "missing_count": missing,
                        "missing_run_count": 1 if missing else 0,
                        "longest_missing_run_samples": missing,
                    }
                )

        # Laugh/Colbeck remains a figure-only family in this readiness report.
        artifacts.extend(
            {
                "dataset_family": "laugh_tests",
                "case_id": case,
                "source": "laugh_tests",
                "kind": "evaluation",
                "path": str(evaluation_root / "laugh_tests" / case / "evaluation.mat"),
                "status": "passed",
            }
            for case in case_map["laugh_tests"]
        )

        # Every accepted artifact carries the exact identity emitted by the
        # first-class audit. Tiny opaque files are sufficient because this test
        # owns report trust, not MAT/JSON payload parsing beyond the manifest.
        for index, artifact in enumerate(artifacts):
            artifact_path = Path(artifact["path"])
            artifact_path.parent.mkdir(parents=True, exist_ok=True)
            if not artifact_path.is_file():
                artifact_path.write_bytes(f"artifact-{index}".encode("utf-8"))
            artifact_bytes = artifact_path.read_bytes()
            artifact.update(
                {
                    "exists": True,
                    "artifact_sha256": hashlib.sha256(artifact_bytes).hexdigest(),
                    "artifact_size_bytes": len(artifact_bytes),
                }
            )
        qa = {
            "schema_version": "1.0",
            "generated_at": generated_at,
            "evaluation_data_root": str(evaluation_root),
            "input_data_root": str(input_root),
            "families": list(report_builder.FAMILIES),
            "artifacts": artifacts,
            "channels": channels,
            "findings": findings,
            "summary": {
                "family_count": len(report_builder.FAMILIES),
                "case_count": sum(len(family_cases) for family_cases in case_map.values()),
                "artifact_count": len(artifacts),
                "channel_count": len(channels),
                "finding_count": len(findings),
                "placeholder_count": len(findings),
                "error_count": 0,
                "blocker_count": 0,
                "warning_count": 0,
                "passed": True,
            },
        }
        qa_path.parent.mkdir(parents=True)
        qa_path.write_text(json.dumps(qa), encoding="utf-8")
        qa_markdown = qa_path.with_suffix(".md")
        qa_markdown.write_text("# Artifact QA\n", encoding="utf-8")
        os.utime(qa_path, (qa_time, qa_time))
        os.utime(qa_markdown, (qa_time, qa_time))

        # PNG validity is owned by plotting; this gate only requires a nonempty
        # current file for every QA-derived case.
        for family, family_cases in case_map.items():
            for case in family_cases:
                image = preview / "figures" / family / case / "plot_one.png"
                image.parent.mkdir(parents=True)
                image.write_bytes(b"\x89PNG\r\n\x1a\n")
                os.utime(image, (qa_time + 1, qa_time + 1))

        # Add the separately audited nested firn tree and the exact MATLAB
        # figure ledger consumed by the combined report.
        firn_qa_path = preview / "qa/firn/artifact_qa.json"
        firn_artifacts = []
        firn_rows = []
        for family in report_builder.FIRN_FAMILIES:
            case = f"{family}_case"
            artifact_path = evaluation_root / family / case / "observations.mat"
            artifact_path.parent.mkdir(parents=True)
            artifact_path.write_bytes(f"{family}-artifact".encode("utf-8"))
            artifact_bytes = artifact_path.read_bytes()
            firn_artifacts.append(
                {
                    "dataset_family": family,
                    "case_id": case,
                    "source": family,
                    "kind": "evaluation",
                    "path": str(artifact_path),
                    "exists": True,
                    "artifact_sha256": hashlib.sha256(artifact_bytes).hexdigest(),
                    "artifact_size_bytes": len(artifact_bytes),
                }
            )
            # Include one interpreted firn figure so the fixture proves that
            # appendix captions carry the same semantics as seasonal figures.
            stem = "surface_albedo" if family == "retmip" else "profile"
            image = preview / "figures/firn" / family / case / f"{stem}.png"
            image.parent.mkdir(parents=True)
            image.write_bytes(b"\x89PNG\r\n\x1a\n")
            os.utime(image, (qa_time + 1, qa_time + 1))
            firn_rows.append(
                {
                    "dataset_family": family,
                    "case_id": case,
                    "figure_file": str(image.resolve()),
                }
            )
        firn_qa = {
            "schema_version": "1.0",
            "generated_at": generated_at,
            "evaluation_data_root": str(evaluation_root),
            "input_data_root": str(input_root),
            "families": list(report_builder.FIRN_FAMILIES),
            "artifacts": firn_artifacts,
            "channels": [],
            "findings": [],
            "summary": {
                "family_count": len(report_builder.FIRN_FAMILIES),
                "case_count": len(report_builder.FIRN_FAMILIES),
                "artifact_count": len(firn_artifacts),
                "channel_count": 0,
                "finding_count": 0,
                "placeholder_count": 0,
                "error_count": 0,
                "blocker_count": 0,
                "warning_count": 0,
                "passed": True,
            },
        }
        firn_qa_path.parent.mkdir(parents=True)
        firn_qa_path.write_text(json.dumps(firn_qa), encoding="utf-8")
        firn_qa_path.with_suffix(".md").write_text("# Firn artifact QA\n", encoding="utf-8")
        os.utime(firn_qa_path, (qa_time, qa_time))
        os.utime(firn_qa_path.with_suffix(".md"), (qa_time, qa_time))
        firn_index = preview / "qa/firn_figure_index.csv"
        with firn_index.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=("dataset_family", "case_id", "figure_file")
            )
            writer.writeheader()
            writer.writerows(firn_rows)
        os.utime(firn_index, (qa_time + 2, qa_time + 2))
        return preview, qa_path, output, qa_time

    def setUp(self) -> None:
        """Allocate one parent temporary directory per test method."""
        self.temp_root = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_root.cleanup)

    def test_build_is_deterministic_and_uses_relative_lazy_links(self) -> None:
        """Generate dynamic inventory text without embedding or copying PNGs."""
        preview, qa_path, output, _ = self.make_fixture()

        # Exercise the public CLI path and capture only its compact status line.
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            status = report_builder.main(
                ["--qa", str(qa_path), "--preview-root", str(preview), "--output", str(output)]
            )
        self.assertEqual(status, 0)
        self.assertIn("with 7 linked figures", stdout.getvalue())
        content = output.read_text(encoding="utf-8")
        self.assertIn(
            'src="../figures/promice/case_a/plot_one.png" loading="lazy"', content
        )
        self.assertIn(
            'src="../figures/firn/retmip/retmip_case/surface_albedo.png" loading="lazy"',
            content,
        )
        self.assertIn("3 cases and 9 artifacts", content)
        self.assertIn("1 informational placeholder, 0 documented warnings", content)
        self.assertIn('href="../qa/artifact_qa.md"', content)
        self.assertIn('href="../qa/artifact_qa.json"', content)
        self.assertIn("## PROMICE native forcing readiness", content)
        self.assertIn("Ready sites: 0 of 1", content)
        self.assertIn("complete window: 0 of 1", content)
        self.assertIn("Manifest forcing ready: **no**", content)
        self.assertIn("Complete windows: None recorded", content)
        self.assertIn(
            "Open audited per-source channel coverage in the final appendix", content
        )
        self.assertIn(
            "Open audited staged-met channel coverage in the final appendix", content
        )
        appendix_start = content.index(
            "## Audited per-source channel coverage appendix"
        )
        self.assertGreater(appendix_start, content.rfind("</figure>"))
        self.assertIn(
            "| PROMICE | met | placeholder | 15 min | `ppt` | — | "
            "0 / 96 | 0.00% | 96 | 1 | 96 |",
            content[appendix_start:],
        )
        self.assertNotIn("| PROMICE | met |", content[:appendix_start])
        self.assertIn("missing samples | contiguous missing runs", content)
        self.assertNotIn("required channels with gaps", content)
        self.assertIn("Native PROMICE userdata is evaluation-support evidence", content)
        self.assertEqual(content.count("1 figure across 1 case."), 7)
        self.assertIn("## Firn artifact QA status", content)
        self.assertIn('href="../qa/firn_figure_index.csv"', content)
        self.assertNotIn("data:image", content)

        # An identical second build must preserve the include mtime and bytes.
        first_mtime = output.stat().st_mtime_ns
        first_content = content
        self.assertEqual(report_builder.build_report(qa_path, preview, output), 7)
        self.assertEqual(output.stat().st_mtime_ns, first_mtime)
        self.assertEqual(output.read_text(encoding="utf-8"), first_content)
        self.assertEqual(len(list((preview / "report").rglob("*.png"))), 0)

    def test_matlab_figure_exports_use_exportgraphics(self) -> None:
        """Keep every MATLAB figure writer on the documented export API."""
        prohibited = re.compile(r"\b(?:saveas|savefig|print)\s*\(")
        offenders = []
        for path in sorted(REPO_ROOT.rglob("*.m")):
            for number, line in enumerate(
                path.read_text(encoding="utf-8").splitlines(), start=1
            ):
                if prohibited.search(line):
                    offenders.append(f"{path.relative_to(REPO_ROOT)}:{number}")
        self.assertEqual(offenders, [])
        style = (REPO_ROOT / "STYLE.local.md").read_text(encoding="utf-8")
        self.assertIn("Use `exportgraphics` for every saved MATLAB figure", style)

    def test_retired_geometry_figures_are_outside_report_inventory(self) -> None:
        """Do not link or count stale geometry plots retired by the renderer."""
        preview, qa_path, output, qa_time = self.make_fixture()

        # Retain both historical group names on disk and make them older than
        # QA; exclusion must happen before stale-file validation and counting.
        case_dir = preview / "figures/promice/case_a"
        retired = [
            case_dir / f"{stem}.png"
            for stem in report_builder.RETIRED_FIGURE_STEMS
        ]
        for image in retired:
            image.write_bytes(b"\x89PNG\r\n\x1a\n")
            os.utime(image, (qa_time - 1, qa_time - 1))

        self.assertEqual(report_builder.build_report(qa_path, preview, output), 7)
        content = output.read_text(encoding="utf-8")
        self.assertTrue(all(image.is_file() for image in retired))
        self.assertNotIn("measurement_geometry", content)

    def test_requires_current_exact_firn_figure_index(self) -> None:
        """Reject missing, escaped, or stale rows in the canonical firn ledger."""
        for mode, message in (
            ("missing", "does not exactly match"),
            ("outside", "outside figures/firn"),
            ("stale", "older than firn artifact-QA"),
        ):
            with self.subTest(mode=mode):
                preview, qa_path, output, qa_time = self.make_fixture()
                index = preview / "qa/firn_figure_index.csv"
                if mode == "stale":
                    os.utime(index, (qa_time - 1, qa_time - 1))
                else:
                    with index.open(encoding="utf-8", newline="") as stream:
                        rows = list(csv.DictReader(stream))
                    if mode == "missing":
                        rows.pop()
                    else:
                        rows[0]["figure_file"] = str(
                            (preview.parent / "outside.png").resolve()
                        )
                    with index.open("w", encoding="utf-8", newline="") as stream:
                        writer = csv.DictWriter(
                            stream,
                            fieldnames=("dataset_family", "case_id", "figure_file"),
                        )
                        writer.writeheader()
                        writer.writerows(rows)
                with self.assertRaisesRegex(ValueError, message):
                    report_builder.build_report(qa_path, preview, output)

    def test_escapes_case_heading_markup_and_controls(self) -> None:
        """Render a QA-authorized metacharacter case as literal Markdown text."""
        special_case = "case_[x](y)#$~^"
        case_map = {
            "promice": (special_case,),
            "esm_snowmip": ("case_a",),
            "laugh_tests": ("case_a",),
        }
        preview, qa_path, output, _ = self.make_fixture(case_map)

        # The heading remains literal while the matching QA-approved path is
        # percent-encoded in its relative image link.
        report_builder.build_report(qa_path, preview, output)
        content = output.read_text(encoding="utf-8")
        self.assertIn("### Case \\[X\\]\\(Y\\)\\#\\$\\~\\^", content)
        self.assertIn(
            "../figures/promice/case_%5Bx%5D%28y%29%23%24~%5E/plot_one.png",
            content,
        )
        anchor = report_builder.coverage_anchor("promice", special_case)
        self.assertIn(f'href="#{anchor}"', content)
        self.assertIn(f'id="{anchor}"', content)
        control_heading = report_builder.markdown_heading("case_[x]\n#next")
        self.assertNotIn("\n", control_heading)
        self.assertIn("\\#Next", control_heading)

    def test_accepts_documented_warnings_and_reports_count(self) -> None:
        """Render reconciled warnings while retaining the passed audit gate."""
        preview, qa_path, output, _ = self.make_fixture()
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
        qa["findings"].extend([{"severity": "warning"}, {"severity": "warning"}])
        qa["summary"]["finding_count"] = len(qa["findings"])
        qa["summary"]["warning_count"] = 2
        qa_path.write_text(json.dumps(qa), encoding="utf-8")
        qa_mtime = qa_path.stat().st_mtime_ns
        os.utime(qa_path.with_suffix(".md"), ns=(qa_mtime, qa_mtime))

        self.assertEqual(report_builder.build_report(qa_path, preview, output), 7)
        content = output.read_text(encoding="utf-8")
        self.assertIn("2 documented warnings, and zero errors or blockers", content)

    def test_reports_manifest_readiness_coverage_and_figure_semantics(self) -> None:
        """Report every requested source leg without confusing samples and runs."""
        case_map = {
            "promice": ("case_a", "case_b"),
            "esm_snowmip": ("case_a",),
            "laugh_tests": ("case_a",),
        }
        preview, qa_path, output, qa_time = self.make_fixture(case_map)

        # Add one figure from each interpretation class; the report must label
        # daily reductions without presenting them as forcing-readiness views.
        case_dir = preview / "figures/promice/case_a"
        for stem in ("met_forcing", "radiation_fluxes", "energy_balance", "surface_albedo"):
            image = case_dir / f"{stem}.png"
            image.write_bytes(b"\x89PNG\r\n\x1a\n")
            os.utime(image, (qa_time + 1, qa_time + 1))

        report_builder.build_report(qa_path, preview, output)
        content = output.read_text(encoding="utf-8")
        self.assertIn("| case_a | no | None recorded | placeholder | 0.00%", content)
        self.assertIn(
            "| case_b | yes | 2020-01-01 00:00:00 to 2020-01-01 23:45:00 (96 samples)",
            content,
        )
        self.assertIn(
            "| MAR 3.11 | met | passed | 15 min | `swd` | W m-2 | "
            "96 / 96 | 100.00% | 0 | 0 | 0 |",
            content,
        )
        self.assertIn(
            "| MAR 3.11 | userdata | passed | 1 h | `swn` | W m-2 | "
            "23 / 24 | 95.83% | 1 | 1 | 1 |",
            content,
        )
        self.assertIn(
            "| MERRA-2 | userdata | passed | 1 h | `swn` | W m-2 | "
            "24 / 24 | 100.00% | 0 | 0 | 0 |",
            content,
        )
        self.assertIn(
            "| RACMO 2.3p3 | userdata | passed | 1 h | `swn` | W m-2 | "
            "22 / 24 | 91.67% | 2 | 1 | 2 |",
            content,
        )
        self.assertIn(
            "| ESM-SnowMIP | met | passed | 15 min | `lwd` | W m-2 | "
            "93 / 96 | 96.88% | 3 | 1 | 3 |",
            content,
        )
        self.assertEqual(content.count('href="#coverage-'), 3)
        self.assertEqual(content.count('<span id="coverage-'), 3)
        self.assertIn("`missing_run_count` is the number of contiguous missing runs", content)
        self.assertIn("Readiness view: staged model-met support", content)
        self.assertIn("Native-support view: actual source timestamps", content)
        self.assertIn("Daily diagnostic: dense subdaily series require complete", content)
        self.assertIn("Daily diagnostic: radiometer observations require a complete", content)
        self.assertIn("ratio of reflected to incoming daily shortwave energy", content)
        self.assertIn("Native MAR/MERRA-2 modeled surface-state albedo", content)
        self.assertIn("valid polar-night coverage", content)
        self.assertIn("RACMO/model-ratio albedo uses positive SWD", content)
        self.assertIn("without the radiometer six-hour gate", content)
        self.assertIn("all-zero-SWD polar-night day remains undefined", content)
        self.assertIn("When MAR userdata lacks a radiation balance", content)
        self.assertIn("icemodel.surface.net_longwave_radiation(tsfc, lwd)", content)
        self.assertIn("same surface emissivity and Kelvin temperature", content)
        self.assertIn("staged artifact remain unchanged", content)
        self.assertNotIn("averages available finite albedo", content)
        self.assertNotIn("Readiness", report_builder.figure_interpretation("energy_balance"))
        self.assertNotIn("Readiness", report_builder.figure_interpretation("surface_albedo"))

    def test_binds_report_to_every_audited_artifact_identity(self) -> None:
        """Fail closed when any audited file identity is absent, stale, or invalid."""
        def remove_file(row: dict) -> None:
            """Delete the selected non-tabulated artifact from disk."""
            Path(row["path"]).unlink()

        def replace_file_bytes(row: dict) -> None:
            """Change file content without changing its audited byte count."""
            path = Path(row["path"])
            content = path.read_bytes()
            path.write_bytes(bytes([content[0] ^ 0xFF]) + content[1:])

        mutations = (
            ("missing sha", lambda row: row.pop("artifact_sha256"), "lacks a valid"),
            (
                "invalid sha",
                lambda row: row.__setitem__("artifact_sha256", "z" * 64),
                "lacks a valid",
            ),
            (
                "sha mismatch",
                lambda row: row.__setitem__("artifact_sha256", "0" * 64),
                "artifact_sha256 does not match",
            ),
            (
                "size mismatch",
                lambda row: row.__setitem__(
                    "artifact_size_bytes", row["artifact_size_bytes"] + 1
                ),
                "artifact_size_bytes does not match",
            ),
            (
                "invalid size",
                lambda row: row.__setitem__("artifact_size_bytes", "invalid"),
                "invalid artifact_size_bytes",
            ),
            (
                "same-size disk drift",
                replace_file_bytes,
                "artifact_sha256 does not match",
            ),
            ("missing file", remove_file, "missing audited artifact"),
            (
                "missing exists flag",
                lambda row: row.pop("exists"),
                "non-existing artifact",
            ),
        )
        for name, mutate, message in mutations:
            with self.subTest(name=name):
                preview, qa_path, output, _ = self.make_fixture()
                qa = json.loads(qa_path.read_text(encoding="utf-8"))
                # The Laugh evaluation never contributes readiness or coverage
                # tables, proving the gate covers every artifact, not only rows
                # otherwise consumed by report generation.
                mutate(qa["artifacts"][-1])
                qa_path.write_text(json.dumps(qa), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    report_builder.build_report(qa_path, preview, output)

    def test_rejects_unreconciled_quantitative_coverage(self) -> None:
        """Reject stale denominators and legacy QA without contiguous-run evidence."""
        mutations = (
            (
                "sample denominator",
                lambda row: row.__setitem__("missing_count", row["missing_count"] + 1),
                "counts do not reconcile",
            ),
            (
                "legacy missing run fields",
                lambda row: row.pop("missing_run_count"),
                "invalid missing_run_count",
            ),
            (
                "run semantics",
                lambda row: row.__setitem__("longest_missing_run_samples", 0),
                "missing-run counts do not reconcile",
            ),
        )
        for name, mutate, message in mutations:
            with self.subTest(name=name):
                preview, qa_path, output, _ = self.make_fixture()
                qa = json.loads(qa_path.read_text(encoding="utf-8"))
                channel = next(
                    row
                    for row in qa["channels"]
                    if row["dataset_family"] == "esm_snowmip" and row["channel"] == "lwd"
                )
                mutate(channel)
                qa_path.write_text(json.dumps(qa), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    report_builder.build_report(qa_path, preview, output)

    def test_rejects_failed_qa_and_unsafe_figure_inventories(self) -> None:
        """Reject failed QA, blocking findings, and unsafe report inputs."""
        # Each mutator produces one requested failure class while retaining an
        # otherwise reconciled fixture, so messages remain diagnostic.
        def failed(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Clear the first-class passed gate without count drift."""
            qa = json.loads(qa_path.read_text(encoding="utf-8"))
            qa["summary"]["passed"] = False
            qa_path.write_text(json.dumps(qa), encoding="utf-8")

        def error(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Replace the placeholder with one reconciled error."""
            qa = json.loads(qa_path.read_text(encoding="utf-8"))
            qa["findings"][0]["severity"] = "error"
            qa["summary"]["error_count"] = 1
            qa["summary"]["placeholder_count"] = 0
            qa_path.write_text(json.dumps(qa), encoding="utf-8")

        def blocker(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Replace the placeholder with one reconciled blocker."""
            qa = json.loads(qa_path.read_text(encoding="utf-8"))
            qa["findings"][0]["severity"] = "blocker"
            qa["summary"]["blocker_count"] = 1
            qa["summary"]["placeholder_count"] = 0
            qa_path.write_text(json.dumps(qa), encoding="utf-8")

        def wrong_family(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Replace one accepted QA family with an unrelated family."""
            qa = json.loads(qa_path.read_text(encoding="utf-8"))
            qa["families"][-1] = "retmip"
            qa_path.write_text(json.dumps(qa), encoding="utf-8")

        def missing(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Remove one case's only figure while retaining its directory."""
            (preview / "figures/promice/case_a/plot_one.png").unlink()

        def zero_byte(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Truncate one current figure."""
            (preview / "figures/promice/case_a/plot_one.png").write_bytes(b"")

        def stale(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Move one figure before the QA generation time."""
            image = preview / "figures/promice/case_a/plot_one.png"
            os.utime(image, (qa_time - 1, qa_time - 1))

        def extra_family(preview: Path, qa_path: Path, qa_time: float) -> None:
            """Add a figure beneath an out-of-scope family."""
            image = preview / "figures/retmip/case_a/plot_one.png"
            image.parent.mkdir(parents=True)
            image.write_bytes(b"png")

        cases = (
            ("failed", failed, "did not pass"),
            ("error", error, "did not pass"),
            ("blocker", blocker, "did not pass"),
            ("wrong family", wrong_family, "must be exactly"),
            ("missing PNG", missing, "missing PNG cases"),
            ("zero byte", zero_byte, "zero-byte PNGs"),
            ("stale", stale, "older than artifact QA"),
            ("extra family", extra_family, "extra PNG families"),
        )
        for name, mutate, message in cases:
            with self.subTest(name=name):
                preview, qa_path, output, qa_time = self.make_fixture()
                mutate(preview, qa_path, qa_time)
                with self.assertRaisesRegex(ValueError, message):
                    report_builder.build_report(qa_path, preview, output)

    def test_requires_exact_qa_case_directory_coverage(self) -> None:
        """Reject both a deleted authorized case and an unexpected case."""
        # A partial deletion must fail even when the family still has another
        # valid figure, preventing family-level checks from masking the case.
        two_promice = {
            "promice": ("case_a", "case_b"),
            "esm_snowmip": ("case_a",),
            "laugh_tests": ("case_a",),
        }
        preview, qa_path, output, _ = self.make_fixture(two_promice)
        shutil.rmtree(preview / "figures/promice/case_b")
        with self.assertRaisesRegex(ValueError, "missing case directories: case_b"):
            report_builder.build_report(qa_path, preview, output)

        # An added case directory and PNG cannot expand report scope beyond the
        # exact cases derived from the same QA artifact records.
        preview, qa_path, output, _ = self.make_fixture()
        image = preview / "figures/promice/case_extra/plot_one.png"
        image.parent.mkdir(parents=True)
        image.write_bytes(b"png")
        with self.assertRaisesRegex(ValueError, "unexpected case directories: case_extra"):
            report_builder.build_report(qa_path, preview, output)

    def test_reconciles_qa_schema_records_counts_and_time(self) -> None:
        """Reject malformed or internally inconsistent first-class QA data."""
        # Scalar and collection replacements exercise every consumed schema
        # boundary without creating a fixture framework or external dependency.
        patches = (
            ("schema", ("schema_version",), "2.0", "schema_version must be 1.0"),
            ("findings type", ("findings",), {}, "invalid artifact-QA JSON"),
            (
                "family entry type",
                ("families",),
                ["promice", "esm_snowmip", 7],
                "family entries must be nonempty strings",
            ),
            ("invalid time", ("generated_at",), "not-a-time", "invalid artifact-QA JSON"),
            (
                "naive time",
                ("generated_at",),
                "2023-11-14T22:13:20",
                "generated_at must include a timezone",
            ),
            ("family count", ("summary", "family_count"), 2, "family count does not match"),
            (
                "severity count",
                ("summary", "placeholder_count"),
                0,
                "severity counts do not match",
            ),
            ("finding count", ("summary", "finding_count"), 0, "finding count does not match"),
            ("artifact count", ("summary", "artifact_count"), 2, "artifact count does not match"),
            ("channel count", ("summary", "channel_count"), 0, "channel count does not match"),
            ("case count", ("summary", "case_count"), 2, "case count does not match"),
        )
        for name, keys, value, message in patches:
            with self.subTest(name=name):
                preview, qa_path, output, _ = self.make_fixture()
                qa = json.loads(qa_path.read_text(encoding="utf-8"))
                target = qa
                for key in keys[:-1]:
                    target = target[key]
                target[keys[-1]] = value
                qa_path.write_text(json.dumps(qa), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    report_builder.build_report(qa_path, preview, output)

        # Invalid record shapes and severities have dedicated diagnostics, and
        # removing a family's records cannot yield an empty expected case set.
        record_cases = (
            "artifact object",
            "artifact family",
            "artifact family type",
            "finding severity",
            "finding severity type",
            "empty cases",
        )
        for name in record_cases:
            with self.subTest(name=name):
                preview, qa_path, output, _ = self.make_fixture()
                qa = json.loads(qa_path.read_text(encoding="utf-8"))
                if name == "artifact object":
                    qa["artifacts"][0] = 7
                elif name == "artifact family":
                    qa["artifacts"][0]["dataset_family"] = "retmip"
                elif name == "artifact family type":
                    qa["artifacts"][0]["dataset_family"] = []
                elif name == "finding severity":
                    qa["findings"][0]["severity"] = "notice"
                elif name == "finding severity type":
                    qa["findings"][0]["severity"] = {}
                else:
                    qa["artifacts"] = [
                        artifact
                        for artifact in qa["artifacts"]
                        if artifact["dataset_family"] != "laugh_tests"
                    ]
                    qa["summary"]["artifact_count"] = len(qa["artifacts"])
                    qa["summary"]["case_count"] -= 1
                qa_path.write_text(json.dumps(qa), encoding="utf-8")
                with self.assertRaises(ValueError):
                    report_builder.build_report(qa_path, preview, output)

    def test_rejects_invalid_paths_and_unreadable_inputs(self) -> None:
        """Keep QA/output paths contained and fail closed on absent inputs."""
        preview, qa_path, output, _ = self.make_fixture()

        # Missing and malformed QA files beneath the authorized QA directory
        # fail before any report output is written.
        missing_qa = qa_path.with_name("missing.json")
        with self.assertRaisesRegex(ValueError, "missing or empty"):
            report_builder.build_report(missing_qa, preview, output)
        qa_path.write_text("not json", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "invalid artifact-QA JSON"):
            report_builder.build_report(qa_path, preview, output)

        # A syntactically valid but structurally incomplete summary produces the
        # same explicit schema error rather than an incidental exception.
        malformed_preview, malformed_qa, malformed_output, _ = self.make_fixture()
        malformed = json.loads(malformed_qa.read_text(encoding="utf-8"))
        malformed["summary"] = []
        malformed_qa.write_text(json.dumps(malformed), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "invalid artifact-QA JSON"):
            report_builder.build_report(malformed_qa, malformed_preview, malformed_output)

        # The report links both audit representations; a missing human-readable
        # summary must not leave a dead evidence link in an otherwise valid build.
        markdown_preview, markdown_qa, markdown_output, _ = self.make_fixture()
        markdown_qa.with_suffix(".md").unlink()
        with self.assertRaisesRegex(ValueError, "missing or empty artifact-QA Markdown"):
            report_builder.build_report(markdown_qa, markdown_preview, markdown_output)

        stale_preview, stale_qa, stale_output, qa_time = self.make_fixture()
        os.utime(stale_qa.with_suffix(".md"), (qa_time - 1, qa_time - 1))
        with self.assertRaisesRegex(ValueError, "Markdown is older"):
            report_builder.build_report(stale_qa, stale_preview, stale_output)

        # CLI path overrides may select files only within the paired QA and
        # report subtrees of the same preview root.
        path_cases = (
            ("QA", preview / "outside.json", output, "beneath <preview-root>/qa"),
            ("output", qa_path, preview / "outside.md", "beneath <preview-root>/report"),
        )
        for name, selected_qa, selected_output, message in path_cases:
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, message):
                    report_builder.main(
                        [
                            "--qa",
                            str(selected_qa),
                            "--preview-root",
                            str(preview),
                            "--output",
                            str(selected_output),
                        ]
                    )

        # The lower-level collector still owns the absent-root diagnostic used
        # after otherwise valid schema and containment checks.
        expected = {family: {"case_a"} for family in report_builder.FAMILIES}
        generated_at = datetime.fromtimestamp(1_700_000_000, timezone.utc)
        with self.assertRaisesRegex(ValueError, "missing figures root"):
            report_builder.collect_figures(
                preview / "figures/absent", generated_at, expected
            )


if __name__ == "__main__":
    unittest.main()
