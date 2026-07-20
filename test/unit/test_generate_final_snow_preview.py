"""Static regression checks for the canonical seasonal-preview driver."""

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
DRIVER = (
    REPO_ROOT
    / "icemodel/+icemodel/+verification/+report/generateFinalSnowPreview.m"
)


class GenerateFinalSnowPreviewTests(unittest.TestCase):
    """Keep seasonal refreshes independent of colocated firn previews."""

    @classmethod
    def setUpClass(cls) -> None:
        """Read the small MATLAB driver once for focused contract checks."""
        cls.source = DRIVER.read_text(encoding="utf-8")

    def test_png_inventory_is_scoped_to_declared_seasonal_families(self) -> None:
        """A firn sibling tree must not inflate the 678-image seasonal count."""
        self.assertIn(
            'pngs = dir(fullfile(figure_root, families(1), "**", "*.png"));',
            self.source,
        )
        self.assertIn("for family_idx = 2:numel(families)", self.source)
        self.assertNotIn(
            'pngs = dir(fullfile(figure_root, "**", "*.png"));', self.source
        )

    def test_driver_keeps_the_exact_seasonal_contract(self) -> None:
        """The scope fix must retain accepted families and exact inventory gate."""
        self.assertIn(
            'families = ["promice", "esm_snowmip", "laugh_tests"];',
            self.source,
        )
        self.assertIn("audit.summary.artifact_count == 347", self.source)
        self.assertNotIn("audit.summary.artifact_count == 348", self.source)
        self.assertIn("&& height(summary) == 678", self.source)
        self.assertIn("writetable(summary", self.source)


if __name__ == "__main__":
    unittest.main()
