"""Focused tests for the seasonal-snow rendered-report integrity checker."""

from __future__ import annotations

import contextlib
import csv
import importlib
import io
import json
import os
import sys
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

# Import the sibling script without making report a production Python package.
REPORT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(REPORT_DIR))
checker = importlib.import_module("check_snow_artifact_qa")


class SnowArtifactReportCheckerTest(unittest.TestCase):
    """Exercise positive and fail-closed report validation paths."""

    def setUp(self) -> None:
        """Create one disposable root for each test procedure."""
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        """Remove every generated report fixture."""
        self.temporary.cleanup()

    def make_fixture(self, name: str = "fixture") -> dict[str, object]:
        """Build one current seven-family report with exact linked images."""
        root = self.root / name
        preview = root / "preview"
        report_dir = preview / "report"
        qa_dir = preview / "qa"
        qmd = root / "source/snow-artifact-qa.qmd"
        report_dir.mkdir(parents=True)
        qa_dir.mkdir(parents=True)
        qmd.parent.mkdir(parents=True)

        # Match the exact canonical 62 seasonal and 56 firn case scopes while
        # retaining stable first records used by focused mutation tests.
        seasonal_records = (
            (("promice", "case_a", "Case A"),)
            + tuple(
                ("promice", f"promice_{index:02d}", f"PROMICE {index:02d}")
                for index in range(2, 52)
            )
            + (("esm_snowmip", "case_b", "Case B"),)
            + tuple(
                ("esm_snowmip", f"esm_{index:02d}", f"ESM {index:02d}")
                for index in range(2, 11)
            )
            + (("laugh_tests", "colbeck", "Colbeck"),)
        )
        firn_records = (
            tuple(
                ("retmip", f"retmip_{index:02d}", f"RetMIP {index:02d}")
                for index in range(1, 6)
            )
            + tuple(
                ("imau", f"imau_{index:02d}", f"IMAU {index:02d}")
                for index in range(1, 4)
            )
            + tuple(
                ("sumup", f"sumup_{index:02d}", f"SUMup {index:02d}")
                for index in range(1, 48)
            )
            + (("research_site", "research", "Research"),)
        )
        records = seasonal_records + firn_records
        sources: list[str] = []
        figures: list[Path] = []
        for family, case, _ in seasonal_records:
            figure = preview / "figures" / family / case / "plot.png"
            figure.parent.mkdir(parents=True)
            figure.write_bytes(b"png")
            figures.append(figure)
            sources.append(f"../figures/{family}/{case}/plot.png")
        firn_rows = []
        for family, case, _ in firn_records:
            figure = preview / "figures/firn" / family / case / "plot.png"
            figure.parent.mkdir(parents=True)
            figure.write_bytes(b"png")
            figures.append(figure)
            sources.append(f"../figures/firn/{family}/{case}/plot.png")
            firn_rows.append(
                {
                    "dataset_family": family,
                    "case_id": case,
                    "figure_file": str(figure.resolve()),
                }
            )

        # The QA timestamp is both a freshness boundary and rendered identity.
        base_time = 1_700_000_000
        generated_at = datetime.fromtimestamp(base_time, timezone.utc).isoformat().replace(
            "+00:00", "Z"
        )
        firn_generated_at = datetime.fromtimestamp(
            base_time + 1, timezone.utc
        ).isoformat().replace("+00:00", "Z")
        qa = qa_dir / "artifact_qa.json"
        qa.write_text(
            json.dumps(
                {
                    "generated_at": generated_at,
                    "families": sorted(checker.FAMILIES),
                    "summary": {"case_count": len(seasonal_records)},
                }
            )
            + "\n",
            encoding="utf-8",
        )
        firn_qa = qa_dir / "firn/artifact_qa.json"
        firn_qa.parent.mkdir(parents=True)
        firn_qa.write_text(
            json.dumps(
                {
                    "generated_at": firn_generated_at,
                    "families": sorted(checker.FIRN_FAMILIES),
                    "summary": {"case_count": len(firn_records)},
                }
            )
            + "\n",
            encoding="utf-8",
        )
        firn_index = qa_dir / "firn_figure_index.csv"
        with firn_index.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=("dataset_family", "case_id", "figure_file")
            )
            writer.writeheader()
            writer.writerows(firn_rows)

        generated = report_dir / "figures.generated.md"
        generated.write_text(
            self.generated_markdown(
                records, sources, generated_at, firn_generated_at
            ),
            encoding="utf-8",
        )
        include = os.path.relpath(generated, qmd.parent).replace(os.sep, "/")
        qmd.write_text(
            "---\n"
            'title: "Snow and firn artifact visual QA"\n'
            "---\n\n"
            "Canonical **tracked** report `prose`.\n\n"
            "```{=html}\n<style>ignored by prose identity</style>\n```\n\n"
            f"{{{{< include {include} >}}}}\n",
            encoding="utf-8",
        )
        report = report_dir / "snow-artifact-qa.html"
        report.write_text(
            self.rendered_html(records, sources, generated_at, firn_generated_at),
            encoding="utf-8",
        )

        fixture: dict[str, object] = {
            "root": root,
            "preview": preview,
            "report": report,
            "generated": generated,
            "qa": qa,
            "firn_qa": firn_qa,
            "firn_index": firn_index,
            "qmd": qmd,
            "figures": figures,
            "sources": sources,
            "records": records,
            "generated_at": generated_at,
            "firn_generated_at": firn_generated_at,
            "base_time": base_time,
        }
        self.refresh_times(fixture)
        return fixture

    @staticmethod
    def generated_markdown(
        records: tuple[tuple[str, str, str], ...],
        sources: list[str],
        generated_at: str,
        firn_generated_at: str,
    ) -> str:
        """Render the deterministic include authority used by a fixture."""
        lines = [
            f"Artifact QA passed at `{generated_at}`.",
            f"Firn artifact QA passed at `{firn_generated_at}`.",
            "",
        ]
        coverage = []
        for (family, case, label), source in zip(records, sources, strict=True):
            if family in {"promice", "esm_snowmip"}:
                anchor = f"coverage-{family}-{case}"
                lines.extend([f'<a href="#{anchor}">Coverage</a>', ""])
                coverage.append((anchor, label))
            lines.extend(
                [
                    f"### {label}",
                    "",
                    f'<a href="{source}"><img src="{source}" loading="lazy"></a>',
                    "",
                ]
            )
        lines.extend(["## Audited per-source channel coverage appendix", ""])
        for anchor, label in coverage:
            lines.extend([f'<span id="{anchor}"></span>', f"#### {label}", ""])
        return "\n".join(lines)

    @staticmethod
    def rendered_html(
        records: tuple[tuple[str, str, str], ...],
        sources: list[str],
        generated_at: str,
        firn_generated_at: str,
    ) -> str:
        """Render the minimal Quarto-like HTML shape checked in production."""
        anchors = "".join(
            f'<a href="#{case.replace("_", "-")}">{label}</a>'
            for _, case, label in records
        )
        figures = []
        coverage = []
        for (family, case, label), source in zip(records, sources, strict=True):
            coverage_link = ""
            if family in {"promice", "esm_snowmip"}:
                anchor = f"coverage-{family}-{case}"
                coverage_link = f'<a href="#{anchor}">Coverage</a>'
                coverage.append(anchor)
            figures.append(
                f'<section id="{case.replace("_", "-")}" class="level3">'
                f'<h3 class="anchored" data-anchor-id="{case.replace("_", "-")}">'
                f"{label}</h3>{coverage_link}"
                f'<a href="{source}"><img src="{source}" loading="lazy"></a>'
                "</section>"
            )
        coverage_spans = "".join(f'<span id="{anchor}"></span>' for anchor in coverage)
        return (
            "<html><body>"
            "<h1>Snow and firn artifact visual QA</h1>"
            "<p>Canonical tracked report prose.</p>"
            f"<p>Artifact QA passed at <code>{generated_at}</code>.</p>"
            f"<p>Firn artifact QA passed at <code>{firn_generated_at}</code>.</p>"
            f'<nav id="TOC">{anchors}</nav>{"".join(figures)}{coverage_spans}'
            "</body></html>\n"
        )

    @staticmethod
    def refresh_times(fixture: dict[str, object]) -> None:
        """Restore a deterministic QA-to-render freshness sequence."""
        base = int(fixture["base_time"])
        os.utime(Path(fixture["qa"]), (base + 10, base + 10))
        os.utime(Path(fixture["firn_qa"]), (base + 11, base + 11))
        os.utime(Path(fixture["firn_index"]), (base + 12, base + 12))
        for figure in fixture["figures"]:
            os.utime(Path(figure), (base + 20, base + 20))
        os.utime(Path(fixture["qmd"]), (base + 30, base + 30))
        os.utime(Path(fixture["generated"]), (base + 40, base + 40))
        os.utime(Path(fixture["report"]), (base + 50, base + 50))

    @staticmethod
    def validate(fixture: dict[str, object]) -> dict[str, int | str]:
        """Run the checker against one fixture."""
        return checker.validate_report(
            Path(fixture["report"]),
            Path(fixture["generated"]),
            Path(fixture["preview"]),
            Path(fixture["qa"]),
            Path(fixture["qmd"]),
        )

    def test_accepts_exact_current_linked_report_and_records_hashes(self) -> None:
        """A current exact report returns deterministic counts and identities."""
        fixture = self.make_fixture()
        result = self.validate(fixture)
        self.assertEqual(result["site_heading_count"], 118)
        self.assertEqual(result["image_count"], 118)
        self.assertEqual(result["broken_image_count"], 0)
        self.assertEqual(result["stale_image_count"], 0)
        self.assertEqual(result["coverage_table_link_count"], 61)
        for key in (
            "qa_sha256",
            "firn_qa_sha256",
            "firn_figure_index_sha256",
            "generated_sha256",
            "qmd_sha256",
            "report_sha256",
        ):
            self.assertRegex(str(result[key]), r"^[0-9a-f]{64}$")

    def test_rejects_missing_or_misplaced_coverage_appendix(self) -> None:
        """Coverage links must target the final post-figure appendix exactly."""
        for mode, message in (
            ("missing-target", "coverage links do not match"),
            ("misplaced", "coverage appendix must follow every figure"),
        ):
            with self.subTest(mode=mode):
                fixture = self.make_fixture(f"coverage-{mode}")
                generated = Path(fixture["generated"])
                text = generated.read_text(encoding="utf-8")
                if mode == "missing-target":
                    text = text.replace(
                        '<span id="coverage-promice-case_a"></span>', "", 1
                    )
                else:
                    heading = "## Audited per-source channel coverage appendix\n"
                    text = text.replace(heading, "", 1).replace(
                        "### Case A", f"{heading}### Case A", 1
                    )
                generated.write_text(text, encoding="utf-8")
                self.refresh_times(fixture)
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(fixture)

    def test_rejects_reordered_nonlazy_and_substituted_images(self) -> None:
        """Rendered image identity must exactly match the generated include."""
        for mode in ("reordered", "nonlazy", "substituted", "anchor"):
            with self.subTest(mode=mode):
                fixture = self.make_fixture(mode)
                report = Path(fixture["report"])
                text = report.read_text(encoding="utf-8")
                sources = list(fixture["sources"])
                if mode == "reordered":
                    first = f'<img src="{sources[0]}" loading="lazy">'
                    second = f'<img src="{sources[1]}" loading="lazy">'
                    text = text.replace(first, "IMAGE_A", 1).replace(second, first, 1)
                    text = text.replace("IMAGE_A", second, 1)
                elif mode == "nonlazy":
                    text = text.replace('loading="lazy"', 'loading="eager"', 1)
                elif mode == "substituted":
                    substitute = (
                        Path(fixture["preview"])
                        / "figures/promice/case_a/substitute.png"
                    )
                    substitute.write_bytes(b"png")
                    replacement = "../figures/promice/case_a/substitute.png"
                    text = text.replace(sources[0], replacement)
                else:
                    text = text.replace(sources[0], "../../outside.png", 1)
                report.write_text(text, encoding="utf-8")
                self.refresh_times(fixture)
                pattern = (
                    "linked-image identities"
                    if mode == "anchor"
                    else "do not match generated include"
                )
                with self.assertRaisesRegex(ValueError, pattern):
                    self.validate(fixture)

    def test_rejects_embedded_remote_absolute_copied_and_traversal_paths(self) -> None:
        """Even matching include/HTML links must use the canonical local path."""
        invalid = (
            "data:image/png;base64,AA==",
            "https://example.com/plot.png",
            "/tmp/plot.png",
            "snow-artifact-qa_files/plot.png",
            "../../outside.png",
        )
        for index, source in enumerate(invalid):
            with self.subTest(source=source):
                fixture = self.make_fixture(f"path-{index}")
                original = list(fixture["sources"])[0]
                if source == "snow-artifact-qa_files/plot.png":
                    copied = Path(fixture["report"]).parent / source
                    copied.parent.mkdir(parents=True)
                    copied.write_bytes(b"png")
                for key in ("generated", "report"):
                    path = Path(fixture[key])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(original, source),
                        encoding="utf-8",
                    )
                self.refresh_times(fixture)
                with self.assertRaisesRegex(ValueError, "broken or noncanonical"):
                    self.validate(fixture)

    def test_rejects_missing_empty_and_stale_figure_files(self) -> None:
        """Every linked figure must be present, nonempty, and postdate QA."""
        for mode in ("missing", "empty", "stale"):
            with self.subTest(mode=mode):
                fixture = self.make_fixture(mode)
                figure = Path(list(fixture["figures"])[0])
                if mode == "missing":
                    figure.unlink()
                elif mode == "empty":
                    figure.write_bytes(b"")
                else:
                    base = int(fixture["base_time"])
                    os.utime(figure, (base - 1, base - 1))
                pattern = "figure files predate" if mode == "stale" else "broken or noncanonical"
                with self.assertRaisesRegex(ValueError, pattern):
                    self.validate(fixture)

    def test_rejects_firn_scope_index_and_freshness_drift(self) -> None:
        """Firn QA, ledger, and nested PNGs remain one exact 56-case input."""
        for mode, message in (
            ("family", "families do not match"),
            ("case-count", "case_count must equal 56"),
            ("index", "does not exactly match"),
            ("stale-figure", "figure files predate"),
        ):
            with self.subTest(mode=mode):
                fixture = self.make_fixture(f"firn-{mode}")
                if mode in {"family", "case-count"}:
                    path = Path(fixture["firn_qa"])
                    qa = json.loads(path.read_text(encoding="utf-8"))
                    if mode == "family":
                        qa["families"][-1] = "unexpected"
                    else:
                        qa["summary"]["case_count"] = 55
                    path.write_text(json.dumps(qa), encoding="utf-8")
                    self.refresh_times(fixture)
                elif mode == "index":
                    path = Path(fixture["firn_index"])
                    with path.open(encoding="utf-8", newline="") as stream:
                        rows = list(csv.DictReader(stream))
                    rows.pop()
                    with path.open("w", encoding="utf-8", newline="") as stream:
                        writer = csv.DictWriter(
                            stream,
                            fieldnames=("dataset_family", "case_id", "figure_file"),
                        )
                        writer.writeheader()
                        writer.writerows(rows)
                    self.refresh_times(fixture)
                else:
                    figure = next(
                        path
                        for path in fixture["figures"]
                        if "/figures/firn/" in Path(path).as_posix()
                    )
                    base = int(fixture["base_time"])
                    os.utime(Path(figure), (base, base))
                with self.assertRaisesRegex(ValueError, message):
                    self.validate(fixture)

    def test_rejects_missing_real_toc_and_missing_site_anchor(self) -> None:
        """Heading self-links cannot replace Quarto's real site-level TOC."""
        for mode in ("no-toc", "missing-anchor"):
            with self.subTest(mode=mode):
                fixture = self.make_fixture(mode)
                report = Path(fixture["report"])
                text = report.read_text(encoding="utf-8")
                if mode == "no-toc":
                    text = text.replace('id="TOC"', 'id="OTHER"')
                else:
                    text = text.replace('<a href="#case-a">Case A</a>', "")
                report.write_text(text, encoding="utf-8")
                self.refresh_times(fixture)
                pattern = "missing nav#TOC" if mode == "no-toc" else "missing from TOC"
                with self.assertRaisesRegex(ValueError, pattern):
                    self.validate(fixture)

    def test_rejects_stale_include_or_rendered_report(self) -> None:
        """The include and HTML must postdate every evidence input they consume."""
        fixture = self.make_fixture("stale-include")
        base = int(fixture["base_time"])
        os.utime(Path(fixture["generated"]), (base + 9, base + 9))
        with self.assertRaisesRegex(ValueError, "generated include predates QA"):
            self.validate(fixture)

        for key in ("qa", "firn_qa", "firn_index", "generated", "qmd"):
            with self.subTest(newer_input=key):
                fixture = self.make_fixture(f"stale-report-{key}")
                input_time = Path(fixture[key]).stat().st_mtime
                os.utime(Path(fixture["report"]), (input_time - 1, input_time - 1))
                with self.assertRaisesRegex(ValueError, "rendered report predates"):
                    self.validate(fixture)

    def test_qmd_prose_preserves_inline_multiplication(self) -> None:
        """Inline-code operators survive Markdown normalization for HTML checks."""
        blocks = checker.qmd_prose_blocks(
            "Runtime met uses `swd * (1 - albedo)` and **audited** inputs.\n"
        )
        self.assertEqual(
            blocks,
            ["Runtime met uses swd * (1 - albedo) and audited inputs."],
        )

    def test_rejects_qa_timestamp_qmd_prose_and_include_identity_drift(self) -> None:
        """QA and tracked QMD content must be identifiable in rendered output."""
        for mode in (
            "include-time",
            "firn-include-time",
            "report-time",
            "firn-report-time",
            "qmd-prose",
            "qmd-include",
        ):
            with self.subTest(mode=mode):
                fixture = self.make_fixture(mode)
                if mode == "include-time":
                    path = Path(fixture["generated"])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(
                            str(fixture["generated_at"]), "2000-01-01T00:00:00Z"
                        ),
                        encoding="utf-8",
                    )
                    pattern = "exact QA timestamp"
                elif mode == "firn-include-time":
                    path = Path(fixture["generated"])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(
                            str(fixture["firn_generated_at"]),
                            "2000-01-01T00:00:01Z",
                        ),
                        encoding="utf-8",
                    )
                    pattern = "exact firn QA timestamp"
                elif mode == "report-time":
                    path = Path(fixture["report"])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(
                            str(fixture["generated_at"]), "2000-01-01T00:00:00Z"
                        ),
                        encoding="utf-8",
                    )
                    pattern = "lacks the QA timestamp"
                elif mode == "firn-report-time":
                    path = Path(fixture["report"])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(
                            str(fixture["firn_generated_at"]),
                            "2000-01-01T00:00:01Z",
                        ),
                        encoding="utf-8",
                    )
                    pattern = "lacks the firn QA timestamp"
                elif mode == "qmd-prose":
                    path = Path(fixture["qmd"])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(
                            "Canonical **tracked** report `prose`.",
                            "Different tracked prose.",
                        ),
                        encoding="utf-8",
                    )
                    pattern = "lacks tracked QMD prose"
                else:
                    path = Path(fixture["qmd"])
                    path.write_text(
                        path.read_text(encoding="utf-8").replace(
                            "figures.generated.md", "other.generated.md"
                        ),
                        encoding="utf-8",
                    )
                    pattern = "does not resolve"
                self.refresh_times(fixture)
                with self.assertRaisesRegex(ValueError, pattern):
                    self.validate(fixture)

    def test_cli_writes_evidence_only_after_success(self) -> None:
        """The first-class command writes its digest ledger beside the HTML."""
        fixture = self.make_fixture("cli")
        output = Path(fixture["report"]).parent / "report_checks.json"
        argv = [
            "--preview-root",
            str(fixture["preview"]),
            "--report",
            str(fixture["report"]),
            "--generated",
            str(fixture["generated"]),
            "--qa",
            str(fixture["qa"]),
            "--qmd",
            str(fixture["qmd"]),
            "--output",
            str(output),
        ]
        with contextlib.redirect_stdout(io.StringIO()):
            checker.main(argv)
        result = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(result["image_count"], 118)
        self.assertEqual(result["qa_generated_at"], fixture["generated_at"])
        self.assertEqual(
            result["firn_qa_generated_at"], fixture["firn_generated_at"]
        )


if __name__ == "__main__":
    unittest.main()
