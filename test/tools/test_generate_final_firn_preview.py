"""Static regression checks for the canonical firn-preview driver."""

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
DRIVER = REPO_ROOT / "test/tools/generate_final_firn_preview.m"
REPORT_README = REPO_ROOT / "report/README.md"


class GenerateFinalFirnPreviewTests(unittest.TestCase):
    """Keep final firn rendering scoped, sealed, and report-compatible."""

    @classmethod
    def setUpClass(cls) -> None:
        """Read the small MATLAB driver once for focused contract checks."""
        cls.source = DRIVER.read_text(encoding="utf-8")
        cls.report_readme = REPORT_README.read_text(encoding="utf-8")

    def test_driver_uses_only_the_four_canonical_firn_families(self) -> None:
        """Seasonal families must never enter the independent firn ledger."""
        self.assertIn(
            'families = ["retmip", "imau", "sumup", "research_site"];',
            self.source,
        )
        self.assertNotIn('"promice"', self.source)
        self.assertIn('fullfile(preview_root, "figures", "firn")', self.source)

    def test_driver_runs_qa_before_replacing_the_figure_tree(self) -> None:
        """An artifact failure must stop before the expensive canonical render."""
        audit = self.source.index(
            "audit = icemodel.verification.auditArtifacts("
        )
        cleanup = self.source.index("rmdir(figure_root, 's')")
        render = self.source.index(
            "summary = icemodel.verification.plotVerificationArtifacts("
        )
        self.assertLess(audit, cleanup)
        self.assertLess(cleanup, render)
        self.assertIn("audit.summary.family_count == 4", self.source)
        self.assertIn("audit.summary.case_count == 56", self.source)
        self.assertIn("audit.summary.artifact_count == 458", self.source)

    def test_driver_seals_exact_ledger_tree_equality(self) -> None:
        """The report input ledger must exactly equal nonempty canonical PNGs."""
        self.assertIn("&& height(summary) == 702", self.source)
        self.assertIn("&& numel(ledger_cases) == 56", self.source)
        self.assertIn(
            "&& isequal(sort(ledger_families), sort(families(:)))",
            self.source,
        )
        self.assertIn("&& all([pngs.bytes] > 0)", self.source)
        self.assertIn(
            "&& isequal(sort(ledger_paths), sort(png_paths))", self.source
        )
        self.assertIn('"FINAL_FIRN_PREVIEW_COMPLETE.json"', self.source)
        self.assertIn('"FINAL_FIRN_PREVIEW_FAILED.txt"', self.source)

    def test_driver_owns_failure_markers_without_staging(self) -> None:
        """Setup and render failures must be sealed without invoking importers."""
        try_start = self.source.index("try\n")
        dependencies = self.source.index("icemodel.dependencies();")
        catch_start = self.source.index("catch err")
        failed_write = self.source.index(
            'writelines(string(getReport(err, "extended", "hyperlinks", "off")),'
        )
        self.assertLess(try_start, dependencies)
        self.assertLess(dependencies, catch_start)
        self.assertLess(catch_start, failed_write)
        self.assertIn(
            '"icemodel-hfc.2.26.8", "evidence"', self.source
        )
        self.assertIn(
            "if isfile(complete_file)\n      delete(complete_file)\n   end",
            self.source[catch_start:],
        )
        for forbidden in (
            "verification.setup.importRetmip(",
            "verification.setup.importImau(",
            "verification.setup.importSumup(",
            "verification.setup.importResearchSites(",
            "verification.setup.previewFirnStaging(",
            "verification.setup.repairRcmArtifactMetadata(",
            "copyfile(",
            "movefile(",
            "tempdir",
        ):
            self.assertNotIn(forbidden, self.source)

    def test_report_sequence_checks_driver_before_rendering(self) -> None:
        """The documented combined workflow must run the focused driver test."""
        static_test = self.report_readme.index(
            "python3 test/tools/test_generate_final_firn_preview.py"
        )
        firn_render = self.report_readme.index(
            "matlab -batch \"run('test/tools/generate_final_firn_preview.m')\""
        )
        report_tests = self.report_readme.index(
            "python3 -m unittest report/test_build_snow_artifact_qa.py"
        )
        report_build = self.report_readme.index(
            "python3 report/build_snow_artifact_qa.py"
        )
        report_check = self.report_readme.index(
            "python3 report/check_snow_artifact_qa.py"
        )
        self.assertLess(static_test, firn_render)
        self.assertLess(firn_render, report_tests)
        self.assertLess(report_tests, report_build)
        self.assertLess(report_build, report_check)


if __name__ == "__main__":
    unittest.main()
