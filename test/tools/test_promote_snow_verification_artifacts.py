"""Focused standard-library tests for the seasonal artifact promotion tool."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


sys.path.insert(0, str(Path(__file__).parent))
import promote_snow_verification_artifacts as tool  # noqa: E402


class PromotionToolTests(unittest.TestCase):
    """Prove manifest selection, full inventory, pruning, and idempotence."""

    def setUp(self):
        """Create isolated candidate and canonical roots for each test."""
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.candidate = self.root / "candidate"
        self.repo = self.root / "repo"
        self.evidence = self.repo / "sandbox" / "evidence"
        self._write_fixture()

    def tearDown(self):
        """Remove the complete isolated tree."""
        self.temporary.cleanup()

    @staticmethod
    def _write(path: Path, value: str):
        """Write one small fixture file with parent creation."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(value, encoding="utf-8")

    def _manifest(self, root: Path, family: str, payload: dict):
        """Write one minimal family manifest."""
        self._write(root / "eval" / family / "manifest.json", json.dumps(payload))

    @staticmethod
    def _snapshot(root: Path):
        """Return path, bytes, and mtime for every regular file below a root."""
        return {
            path.relative_to(root).as_posix(): (
                path.read_bytes(),
                path.stat().st_mtime_ns,
            )
            for path in root.rglob("*")
            if path.is_file()
        }

    def _write_fixture(self):
        """Build three candidate families and superseded canonical references."""
        # Candidate manifests cover each supported reference field plus inferred
        # atomic ESM met selection.
        self._manifest(self.candidate, "promice", {"cases": [{
            "evaluation_file": "new/observations.mat",
            "colocation": {"promice": {
                "met_files": "promice/met_new.mat",
                "data_files": "promice/data_new.mat",
            }},
        }]})
        self._manifest(self.candidate, "esm_snowmip", {"cases": [{
            "evaluation_file": "obs/observations.mat",
            "obs_file": "",
        }]})
        self._manifest(self.candidate, "laugh_tests", {"cases": [{
            "evaluation_file": "case/evaluation.mat",
            "reference_file": "case/reference.mat",
        }]})
        self._manifest(self.candidate, "retmip", {"cases": [{
            "evaluation_file": "new/observations.mat",
            "colocation": {
                "mar": {
                    "met_files": "retmip/met_new.mat",
                    "data_files": "retmip/data_new.mat",
                    "model_output_files": ["retmip/model_output_new.mat"],
                },
                "retmip": {
                    "model_output_files": [
                        "/Users/example/native/RetMIP_model_output.nc",
                        "C:\\native\\RetMIP_model_output.nc",
                    ],
                },
            },
        }]})
        for relative in [
            "eval/promice/new/observations.mat",
            "input/met/promice/met_new.mat",
            "input/userdata/promice/data_new.mat",
            "eval/esm_snowmip/obs/observations.mat",
            "input/met/esm_snowmip/met_new.mat",
            "eval/laugh_tests/case/evaluation.mat",
            "eval/laugh_tests/case/reference.mat",
            "eval/retmip/new/observations.mat",
            "input/met/retmip/met_new.mat",
            "input/userdata/retmip/data_new.mat",
            "input/userdata/retmip/model_output_new.mat",
            "preview/figures/stale.png",
            "diagnostic.log",
        ]:
            self._write(self.candidate / relative, "candidate:" + relative)

        # Canonical PROMICE owns old evaluation/userdata paths, while a foreign
        # RetMIP manifest shares the old met path and therefore protects it.
        self._manifest(self.repo / "data", "promice", {"cases": [{
            "evaluation_file": "old/observations.mat",
            "colocation": {"promice": {
                "met_files": "promice/met_shared.mat",
                "data_files": "promice/data_old.mat",
                "model_output_files": "retmip/model_output_shared.mat",
            }},
        }]})
        self._manifest(self.repo / "data", "esm_snowmip", {"cases": []})
        self._manifest(self.repo / "data", "laugh_tests", {"cases": []})
        self._manifest(self.repo / "data", "retmip", {"cases": [
            {
                "evaluation_file": "old/observations.mat",
                "colocation": {"mar": {
                    "met_files": "retmip/met_old.mat",
                    "data_files": "retmip/data_old.mat",
                    "model_output_files": [
                        "retmip/model_output_old.mat",
                        "retmip/model_output_shared.mat",
                    ],
                }},
            },
            {
                "colocation": {
                    "promice": {"met_files": "promice/met_shared.mat"},
                },
            },
        ]})
        for relative in [
            "eval/promice/old/observations.mat",
            "input/met/promice/met_shared.mat",
            "input/userdata/promice/data_old.mat",
            "input/met/esm_snowmip/met_old.mat",
            "eval/retmip/old/observations.mat",
            "input/met/retmip/met_old.mat",
            "input/userdata/retmip/data_old.mat",
            "input/userdata/retmip/model_output_old.mat",
            "input/userdata/retmip/model_output_shared.mat",
        ]:
            self._write(self.repo / "data" / relative, "canonical:" + relative)

    def test_inventory_covers_every_candidate_file(self):
        """Every candidate file receives exactly one explicit disposition."""
        rows, selected, pruned = tool.build_inventory(self.candidate, self.repo)
        candidate_count = sum(path.is_file() for path in self.candidate.rglob("*"))
        self.assertEqual(len(rows), candidate_count)
        self.assertEqual(sum(row["classification"] == "promote" for row in rows), len(selected))
        self.assertTrue(any(row["classification"] == "reproducible" for row in rows))
        self.assertTrue(any(row["classification"] == "disposable" for row in rows))
        self.assertIn(Path("eval/promice/old/observations.mat"), pruned)
        self.assertIn(Path("input/userdata/promice/data_old.mat"), pruned)
        self.assertIn(Path("input/met/esm_snowmip/met_old.mat"), pruned)
        self.assertNotIn(Path("input/met/promice/met_shared.mat"), pruned)
        summary = tool.write_evidence(rows, pruned, self.evidence, "test")
        self.assertNotIn("target_families", summary)
        reasons = {
            row["reason"]
            for row in rows
            if row["classification"] == "promote"
        }
        self.assertEqual(
            reasons, {"manifest-selected seasonal-snow artifact"})

    def test_relative_model_output_is_selected_and_absolute_source_is_ignored(self):
        """Only staged relative model outputs join the userdata promotion set."""
        selected = tool.selected_artifacts(self.candidate, ("retmip",))

        self.assertIn(
            Path("input/userdata/retmip/model_output_new.mat"), selected)
        self.assertFalse(any("RetMIP_model_output.nc" in str(path) for path in selected))

    def test_missing_relative_model_output_fails_selection_preflight(self):
        """A manifest-selected relative output must exist in the candidate."""
        manifest_path = self.candidate / "eval/retmip/manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["cases"][0]["colocation"]["mar"]["model_output_files"].append(
            "retmip/missing_model_output.mat")
        self._write(manifest_path, json.dumps(manifest))

        with self.assertRaisesRegex(FileNotFoundError, "missing_model_output"):
            tool.selected_artifacts(self.candidate, ("retmip",))

    def test_escaping_relative_model_output_is_rejected(self):
        """Relative model-output provenance cannot escape candidate userdata."""
        manifest_path = self.candidate / "eval/retmip/manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["cases"][0]["colocation"]["mar"]["model_output_files"] = [
            "../escape.mat",
        ]
        self._write(manifest_path, json.dumps(manifest))

        with self.assertRaisesRegex(ValueError, "unsafe manifest artifact path"):
            tool.selected_artifacts(self.candidate, ("retmip",))

    def test_clone_failure_quietly_falls_back_to_copy(self):
        """A failed native clone is silent and falls back to a verified copy."""
        source = self.root / "source.bin"
        target = self.root / "target" / "artifact.bin"
        self._write(source, "clone fallback payload")
        expected_hash = tool.sha256_file(source)

        with mock.patch.object(
                tool.subprocess, "run",
                return_value=mock.Mock(returncode=1)) as clone_command, \
                mock.patch.object(
                    tool.shutil, "copy2", wraps=tool.shutil.copy2) as copy_fallback:
            tool.clone_or_copy(source, target, expected_hash)

        clone_command.assert_called_once_with(
            ["/bin/cp", "-c", str(source), str(target) + ".promotion.tmp"],
            check=False,
            stdout=tool.subprocess.DEVNULL,
            stderr=tool.subprocess.DEVNULL,
        )
        copy_fallback.assert_called_once_with(
            source, target.with_name(target.name + ".promotion.tmp"))
        self.assertEqual(target.read_bytes(), source.read_bytes())
        self.assertEqual(tool.sha256_file(target), expected_hash)

    def test_copy_replace_failure_removes_temporary_and_preserves_target(self):
        """An interrupted replacement leaves the canonical tree byte-identical."""
        source = self.root / "source.bin"
        target = self.root / "target" / "artifact.bin"
        self._write(source, "replacement payload")
        self._write(target, "canonical payload")
        expected_hash = tool.sha256_file(source)
        before = self._snapshot(self.root)

        with mock.patch.object(
                tool.subprocess, "run",
                return_value=mock.Mock(returncode=1)), \
                mock.patch.object(tool.os, "replace", side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                tool.clone_or_copy(source, target, expected_hash)

        self.assertEqual(self._snapshot(self.root), before)
        self.assertFalse(
            target.with_name(target.name + ".promotion.tmp").exists())

    def test_copy_hash_mismatch_removes_temporary_and_preserves_target(self):
        """A failed copy verification leaves no unowned transaction artifact."""
        source = self.root / "source.bin"
        target = self.root / "target" / "artifact.bin"
        self._write(source, "replacement payload")
        self._write(target, "canonical payload")
        before = self._snapshot(self.root)

        with mock.patch.object(
                tool.subprocess, "run",
                return_value=mock.Mock(returncode=1)):
            with self.assertRaisesRegex(RuntimeError, "copied artifact hash mismatch"):
                tool.clone_or_copy(source, target, "0" * 64)

        self.assertEqual(self._snapshot(self.root), before)
        self.assertFalse(
            target.with_name(target.name + ".promotion.tmp").exists())

    def test_second_apply_is_exact_no_op_and_preserves_foreign_reference(self):
        """A second transaction changes no file and shared foreign data survive."""
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo)
        tool.write_evidence(rows, pruned, self.evidence, "test")
        tool.apply_plan(rows, pruned, self.repo, self.evidence)

        # Selected bytes match, superseded owned files are gone, and the shared
        # foreign met remains untouched.
        for row in rows:
            if row["classification"] == "promote":
                self.assertEqual(
                    tool.sha256_file(Path(row["source_path"])),
                    tool.sha256_file(Path(row["target_path"])),
                )
        self.assertFalse((self.repo / "data/eval/promice/old/observations.mat").exists())
        self.assertFalse((self.repo / "data/input/userdata/promice/data_old.mat").exists())
        self.assertTrue((self.repo / "data/input/met/promice/met_shared.mat").is_file())

        after_rows, _, after_pruned = tool.build_inventory(self.candidate, self.repo)
        promoted = [row for row in after_rows if row["classification"] == "promote"]
        self.assertTrue(all(row["canonical_status"] == "identical" for row in promoted))
        self.assertEqual(after_pruned, [])
        self.assertTrue((self.evidence / tool.TRANSACTION_BACKUP).is_dir())
        self.assertTrue((self.evidence / tool.TRANSACTION_ACTIVE).is_file())
        tool.finalize_transaction(self.repo, self.evidence)
        self.assertFalse((self.evidence / tool.TRANSACTION_BACKUP).exists())
        self.assertFalse((self.evidence / tool.TRANSACTION_ACTIVE).exists())
        self.assertTrue((self.evidence / tool.TRANSACTION_COMPLETE).is_file())

        # Applying the converged plan exercises the public second-pass path and
        # must preserve bytes and mtimes throughout canonical data.
        before_second = self._snapshot(self.repo / "data")
        self.assertFalse(
            tool.apply_plan(after_rows, after_pruned, self.repo, self.evidence))
        self.assertEqual(self._snapshot(self.repo / "data"), before_second)

    def test_initial_noop_apply_seals_empty_transaction(self):
        """A fresh converged promotion emits a seal without touching data."""
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo)
        tool.write_evidence(rows, pruned, self.evidence, "test")
        tool.apply_plan(rows, pruned, self.repo, self.evidence)
        after_rows, _, after_pruned = tool.build_inventory(self.candidate, self.repo)
        tool.finalize_transaction(self.repo, self.evidence)
        noop_evidence = self.repo / "sandbox" / "noop-evidence"
        before = self._snapshot(self.repo / "data")

        self.assertFalse(
            tool.apply_plan(after_rows, after_pruned, self.repo, noop_evidence))

        self.assertEqual(self._snapshot(self.repo / "data"), before)
        seal = json.loads(
            (noop_evidence / tool.TRANSACTION_COMPLETE).read_text(encoding="utf-8"))
        self.assertEqual(seal["state"], "complete")
        self.assertEqual(seal["operations"], [])

    def test_firn_family_selection_protects_seasonal_and_prunes_owned_stale(self):
        """Explicit firn selection replaces only owned, unshared target paths."""
        families = ("retmip",)
        before = self._snapshot(self.repo / "data")
        rows, selected, pruned = tool.build_inventory(
            self.candidate, self.repo, families)

        self.assertEqual(selected, {
            Path("eval/retmip/manifest.json"),
            Path("eval/retmip/new/observations.mat"),
            Path("input/met/retmip/met_new.mat"),
            Path("input/userdata/retmip/data_new.mat"),
            Path("input/userdata/retmip/model_output_new.mat"),
        })
        self.assertIn(Path("eval/retmip/old/observations.mat"), pruned)
        self.assertIn(Path("input/met/retmip/met_old.mat"), pruned)
        self.assertIn(Path("input/userdata/retmip/data_old.mat"), pruned)
        self.assertIn(Path("input/userdata/retmip/model_output_old.mat"), pruned)
        self.assertNotIn(Path("input/userdata/retmip/model_output_shared.mat"), pruned)
        self.assertNotIn(Path("eval/promice/old/observations.mat"), pruned)
        self.assertNotIn(Path("input/met/esm_snowmip/met_old.mat"), pruned)
        summary = tool.write_evidence(
            rows, pruned, self.evidence, "test", families)
        self.assertEqual(summary["target_families"], ["retmip"])

        tool.apply_plan(rows, pruned, self.repo, self.evidence)
        after = self._snapshot(self.repo / "data")
        protected = (
            "eval/promice/manifest.json",
            "eval/promice/old/observations.mat",
            "eval/esm_snowmip/manifest.json",
            "eval/laugh_tests/manifest.json",
            "input/met/esm_snowmip/met_old.mat",
            "input/met/promice/met_shared.mat",
            "input/userdata/promice/data_old.mat",
            "input/userdata/retmip/model_output_shared.mat",
        )
        for relative in protected:
            self.assertEqual(after[relative], before[relative])
        tool.finalize_transaction(self.repo, self.evidence)

    def test_apply_rolls_back_changed_new_and_pruned_paths(self):
        """A post-mutation failure restores the exact pre-transaction inventory."""
        families = ("retmip",)
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo, families)
        before = self._snapshot(self.repo / "data")
        real_write = tool.write_transaction

        def fail_before_commit(path, transaction):
            """Fail after mutation while leaving recovery writes available."""
            if transaction.get("state") == "committed":
                raise RuntimeError("injected post-promotion failure")
            return real_write(path, transaction)

        with mock.patch.object(
                tool, "write_transaction", side_effect=fail_before_commit):
            with self.assertRaisesRegex(RuntimeError, "injected"):
                tool.apply_plan(rows, pruned, self.repo, self.evidence)

        self.assertEqual(self._snapshot(self.repo / "data"), before)
        self.assertFalse((self.evidence / tool.TRANSACTION_BACKUP).exists())
        self.assertFalse((self.evidence / tool.TRANSACTION_ACTIVE).exists())
        self.assertTrue((self.evidence / tool.TRANSACTION_RECOVERED).is_file())

    def test_interrupted_apply_is_recoverable_from_sealed_disk_state(self):
        """Process loss after all mutations retains enough state for recovery."""
        families = ("retmip",)
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo, families)
        before = self._snapshot(self.repo / "data")
        real_write = tool.write_transaction

        def interrupt_before_commit(path, transaction):
            """Simulate SIGKILL semantics after postchecks but before commit."""
            if transaction.get("state") == "committed":
                raise KeyboardInterrupt("injected process loss")
            return real_write(path, transaction)

        with mock.patch.object(
                tool, "write_transaction", side_effect=interrupt_before_commit):
            with self.assertRaisesRegex(KeyboardInterrupt, "process loss"):
                tool.apply_plan(rows, pruned, self.repo, self.evidence)

        self.assertNotEqual(self._snapshot(self.repo / "data"), before)
        self.assertTrue((self.evidence / tool.TRANSACTION_ACTIVE).is_file())
        self.assertTrue((self.evidence / tool.TRANSACTION_BACKUP).is_dir())
        transaction = json.loads(
            (self.evidence / tool.TRANSACTION_ACTIVE).read_text(encoding="utf-8"))
        manifest_flags = [
            operation["target_relative"].endswith("manifest.json")
            for operation in transaction["operations"]
        ]
        first_manifest = manifest_flags.index(True)
        self.assertTrue(all(manifest_flags[first_manifest:]))
        tool.recover_transaction(self.repo, self.evidence)
        self.assertEqual(self._snapshot(self.repo / "data"), before)
        self.assertFalse((self.evidence / tool.TRANSACTION_ACTIVE).exists())
        self.assertFalse((self.evidence / tool.TRANSACTION_BACKUP).exists())
        self.assertTrue((self.evidence / tool.TRANSACTION_RECOVERED).is_file())

    def test_apply_rejects_stale_canonical_plan_before_backup_or_mutation(self):
        """Canonical drift after inventory fails before transaction preparation."""
        families = ("retmip",)
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo, families)
        manifest = self.repo / "data/eval/retmip/manifest.json"
        self._write(manifest, "external canonical mutation")

        with self.assertRaisesRegex(RuntimeError, "stale promotion target plan"):
            tool.apply_plan(rows, pruned, self.repo, self.evidence)

        self.assertEqual(manifest.read_text(encoding="utf-8"),
                         "external canonical mutation")
        self.assertFalse((self.evidence / tool.TRANSACTION_ACTIVE).exists())
        self.assertFalse((self.evidence / tool.TRANSACTION_BACKUP).exists())

    def test_recover_finishes_interrupted_finalize_before_active_unlink(self):
        """A sealed COMPLETE state finishes cleanup without rolling back."""
        families = ("retmip",)
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo, families)
        tool.apply_plan(rows, pruned, self.repo, self.evidence)
        accepted = self._snapshot(self.repo / "data")
        active_path = self.evidence / tool.TRANSACTION_ACTIVE
        real_unlink = Path.unlink

        def interrupt_active_unlink(path, *args, **kwargs):
            """Stop after COMPLETE is sealed but before active is removed."""
            if path == active_path:
                raise KeyboardInterrupt("injected finalize interruption")
            return real_unlink(path, *args, **kwargs)

        with mock.patch.object(Path, "unlink", new=interrupt_active_unlink):
            with self.assertRaisesRegex(KeyboardInterrupt, "finalize"):
                tool.finalize_transaction(self.repo, self.evidence)

        self.assertTrue((self.evidence / tool.TRANSACTION_COMPLETE).is_file())
        self.assertTrue(active_path.is_file())
        self.assertTrue((self.evidence / tool.TRANSACTION_BACKUP).is_dir())
        recovered = tool.recover_transaction(self.repo, self.evidence)
        self.assertEqual(recovered["state"], "complete")
        self.assertEqual(self._snapshot(self.repo / "data"), accepted)
        self.assertFalse(active_path.exists())
        self.assertFalse((self.evidence / tool.TRANSACTION_BACKUP).exists())

    def test_recover_finishes_interrupted_finalize_after_active_unlink(self):
        """A COMPLETE seal also retires a backup left after active removal."""
        families = ("retmip",)
        rows, _, pruned = tool.build_inventory(self.candidate, self.repo, families)
        tool.apply_plan(rows, pruned, self.repo, self.evidence)
        accepted = self._snapshot(self.repo / "data")
        backup_root = self.evidence / tool.TRANSACTION_BACKUP
        real_rmtree = tool.shutil.rmtree

        def interrupt_backup_cleanup(path, *args, **kwargs):
            """Stop after active removal but before backup retirement."""
            if Path(path) == backup_root:
                raise KeyboardInterrupt("injected backup cleanup interruption")
            return real_rmtree(path, *args, **kwargs)

        with mock.patch.object(
                tool.shutil, "rmtree", side_effect=interrupt_backup_cleanup):
            with self.assertRaisesRegex(KeyboardInterrupt, "backup cleanup"):
                tool.finalize_transaction(self.repo, self.evidence)

        self.assertTrue((self.evidence / tool.TRANSACTION_COMPLETE).is_file())
        self.assertFalse((self.evidence / tool.TRANSACTION_ACTIVE).exists())
        self.assertTrue(backup_root.is_dir())
        recovered = tool.recover_transaction(self.repo, self.evidence)
        self.assertEqual(recovered["state"], "complete")
        self.assertEqual(self._snapshot(self.repo / "data"), accepted)
        self.assertFalse(backup_root.exists())

    def test_family_selection_rejects_empty_and_escaping_names(self):
        """Target-family selection cannot be empty or escape the evaluation root."""
        with self.assertRaisesRegex(ValueError, "at least one"):
            tool.selected_artifacts(self.candidate, ())
        with self.assertRaisesRegex(ValueError, "unsafe target family"):
            tool.selected_artifacts(self.candidate, ("../retmip",))


if __name__ == "__main__":
    unittest.main()
