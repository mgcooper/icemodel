#!/usr/bin/env python3
"""Inventory and transactionally promote a manifest-selected artifact candidate.

The default selection remains PROMICE, ESM-SnowMIP, and Laugh Tests for the
seasonal-snow workflow. Callers may instead name another explicit family set,
such as the four firn families, without introducing a second promotion path.
Every candidate file is inventoried so diagnostics and backups have an explicit
durable disposition.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path, PureWindowsPath


TARGET_FAMILIES = ("promice", "esm_snowmip", "laugh_tests")
EVAL_FIELDS = {"evaluation_file", "obs_file", "reference_file"}
MET_FIELDS = {"met_files"}
DATA_FIELDS = {"data_files"}
MODEL_OUTPUT_FIELDS = {"model_output_files"}
TRANSACTION_ACTIVE = "promotion_transaction.json"
TRANSACTION_COMPLETE = "promotion_transaction.COMPLETE.json"
TRANSACTION_RECOVERED = "promotion_transaction.RECOVERED.json"
TRANSACTION_BACKUP = "transaction_backup"


def sha256_file(path: Path) -> str:
    """Return the streaming SHA-256 digest of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_strings(value):
    """Yield scalar strings from decoded JSON containers."""
    if isinstance(value, str) and value.strip():
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from iter_strings(item)
    elif isinstance(value, dict):
        for item in value.values():
            yield from iter_strings(item)


def safe_relative(value: str) -> Path:
    """Normalize a manifest path and reject absolute or escaping references."""
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"unsafe manifest artifact path: {value}")
    return path


def manifest_references(manifest_path: Path, family: str) -> set[Path]:
    """Return candidate-root-relative artifacts named by one family manifest."""
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    references: set[Path] = set()

    # Field names, rather than arbitrary .mat strings, separate staged outputs
    # from raw source files. Relative model outputs are staged userdata, while
    # absolute RetMIP comparison paths remain external source provenance.
    def visit(value, key=""):
        if key in EVAL_FIELDS:
            for item in iter_strings(value):
                references.add(Path("eval") / family / safe_relative(item))
            return
        if key in MET_FIELDS:
            for item in iter_strings(value):
                references.add(Path("input") / "met" / safe_relative(item))
            return
        if key in DATA_FIELDS:
            for item in iter_strings(value):
                references.add(Path("input") / "userdata" / safe_relative(item))
            return
        if key in MODEL_OUTPUT_FIELDS:
            for item in iter_strings(value):
                if Path(item).is_absolute() or PureWindowsPath(item).is_absolute():
                    continue
                references.add(Path("input") / "userdata" / safe_relative(item))
            return
        if isinstance(value, dict):
            for child_key, child in value.items():
                visit(child, child_key)
        elif isinstance(value, list):
            for child in value:
                visit(child, key)

    visit(manifest)
    return references


def selected_artifacts(
        candidate_root: Path,
        target_families=TARGET_FAMILIES,
) -> set[Path]:
    """Build the exact manifest-backed artifact selection for named families."""
    families = tuple(dict.fromkeys(target_families))
    if not families:
        raise ValueError("at least one target family is required")
    selected: set[Path] = set()
    for family in families:
        if not family or Path(family).parts != (family,):
            raise ValueError(f"unsafe target family: {family}")
        manifest = candidate_root / "eval" / family / "manifest.json"
        if not manifest.is_file():
            raise FileNotFoundError(f"missing candidate manifest: {manifest}")
        selected.add(manifest.relative_to(candidate_root))
        selected.update(manifest_references(manifest, family))

    # Atomic ESM-SnowMIP manifests resolve their met files by naming convention,
    # so its family-owned directory joins the selection only when requested.
    if "esm_snowmip" in families:
        esm_met = candidate_root / "input" / "met" / "esm_snowmip"
        selected.update(path.relative_to(candidate_root) for path in esm_met.glob("*.mat"))
    missing = [str(path) for path in selected if not (candidate_root / path).is_file()]
    if missing:
        raise FileNotFoundError("missing selected candidate artifact(s): " + ", ".join(missing))
    return selected


def candidate_disposition(
        relative: Path,
        selected: set[Path],
        target_families=TARGET_FAMILIES,
) -> tuple[str, str]:
    """Classify one candidate file as promoted, reproducible, or disposable."""
    if relative in selected:
        if tuple(target_families) == TARGET_FAMILIES:
            return "promote", "manifest-selected seasonal-snow artifact"
        return "promote", "manifest-selected verification artifact"
    if relative.parts and relative.parts[0] == "preview":
        return (
            "reproducible",
            "regenerated from promoted canonical roots by repository QA/report tools",
        )
    if relative.parts and relative.parts[0] in {"eval", "input"}:
        return "disposable", "unreferenced candidate artifact or filesystem metadata"
    return (
        "disposable",
        "temporary diagnostic, backup, log, cache, or replay wrapper; durable "
        "result is in repository source, tests, Beads, or ExecPlan",
    )


def canonical_status(source_hash: str, target: Path) -> tuple[str, str]:
    """Compare one selected source with its canonical target."""
    if not target.exists():
        return "new", ""
    if not target.is_file() or target.is_symlink():
        raise ValueError(f"canonical target is not a regular file: {target}")
    target_hash = sha256_file(target)
    return ("identical" if target_hash == source_hash else "changed"), target_hash


def all_manifest_references(eval_root: Path, excluded: set[str]) -> set[Path]:
    """Collect references owned by canonical families outside this promotion."""
    references: set[Path] = set()
    if not eval_root.is_dir():
        return references
    for manifest in eval_root.glob("*/manifest.json"):
        family = manifest.parent.name
        if family not in excluded:
            references.update(manifest_references(manifest, family))
    return references


def prune_plan(
        candidate_root: Path,
        repo_root: Path,
        selected: set[Path],
        target_families=TARGET_FAMILIES,
) -> list[Path]:
    """Find superseded target-family files without pruning foreign references."""
    data_root = repo_root / "data"
    old: set[Path] = set()
    families = tuple(dict.fromkeys(target_families))
    for family in families:
        manifest = data_root / "eval" / family / "manifest.json"
        if manifest.is_file():
            old.update(manifest_references(manifest, family))

        # Evaluation directories are family-owned, so stale files not selected
        # by the replacement manifest are safe candidates for exact pruning.
        family_root = data_root / "eval" / family
        if family_root.is_dir():
            old.update(
                path.relative_to(data_root)
                for path in family_root.rglob("*")
                if path.is_file()
            )

    if "esm_snowmip" in families:
        esm_root = data_root / "input" / "met" / "esm_snowmip"
        if esm_root.is_dir():
            old.update(path.relative_to(data_root) for path in esm_root.glob("*.mat"))

    foreign = all_manifest_references(data_root / "eval", set(families))
    return sorted(path for path in old - selected - foreign if (data_root / path).is_file())


def build_inventory(
        candidate_root: Path,
        repo_root: Path,
        target_families=TARGET_FAMILIES,
):
    """Hash and classify every candidate file and return the promotion plan."""
    families = tuple(dict.fromkeys(target_families))
    selected = selected_artifacts(candidate_root, families)
    rows = []
    data_root = repo_root / "data"
    for source in sorted(path for path in candidate_root.rglob("*") if path.is_file()):
        relative = source.relative_to(candidate_root)
        source_hash = sha256_file(source)
        classification, reason = candidate_disposition(
            relative, selected, families)
        target = data_root / relative if classification == "promote" else None
        status, target_hash = canonical_status(source_hash, target) if target else ("", "")
        rows.append({
            "source_relative": relative.as_posix(),
            "source_path": str(source),
            "bytes": source.stat().st_size,
            "sha256": source_hash,
            "classification": classification,
            "canonical_status": status,
            "target_relative": (Path("data") / relative).as_posix() if target else "",
            "target_path": str(target) if target else "",
            "target_sha256_before": target_hash,
            "reason": reason,
        })
    pruned = prune_plan(candidate_root, repo_root, selected, families)
    return rows, selected, pruned


def write_evidence(
        rows,
        pruned,
        evidence_dir: Path,
        mode: str,
        target_families=TARGET_FAMILIES,
) -> dict:
    """Write the complete CSV inventory and compact machine-readable summary."""
    evidence_dir.mkdir(parents=True, exist_ok=True)
    inventory_path = evidence_dir / "candidate_disposition_inventory.csv"
    with inventory_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    counts = {}
    for row in rows:
        key = row["classification"]
        counts[key] = counts.get(key, 0) + 1
    statuses = {}
    for row in rows:
        if row["canonical_status"]:
            key = row["canonical_status"]
            statuses[key] = statuses.get(key, 0) + 1
    summary = {
        "mode": mode,
        "candidate_file_count": len(rows),
        "candidate_byte_count": sum(row["bytes"] for row in rows),
        "classification_counts": counts,
        "canonical_status_counts": statuses,
        "prune_count": len(pruned),
        "prune_paths": [str(Path("data") / path) for path in pruned],
        "inventory_sha256": sha256_file(inventory_path),
    }
    if tuple(target_families) != TARGET_FAMILIES:
        summary["target_families"] = list(dict.fromkeys(target_families))
    (evidence_dir / "promotion_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def clone_or_copy(source: Path, target: Path, expected_hash: str) -> None:
    """Copy through a temporary file, preferring native APFS clone-on-write."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".promotion.tmp")
    temporary.unlink(missing_ok=True)
    try:
        try:
            result = subprocess.run(
                ["/bin/cp", "-c", str(source), str(temporary)],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            cloned = result.returncode == 0
        except OSError:
            cloned = False
        if not cloned:
            shutil.copy2(source, temporary)
        if sha256_file(temporary) != expected_hash:
            raise RuntimeError(f"copied artifact hash mismatch: {target}")
        os.replace(temporary, target)
    finally:
        # Recovery owns sealed targets, so never leave an untracked copy artifact.
        temporary.unlink(missing_ok=True)


def backup_file(target: Path, repo_root: Path, backup_root: Path) -> Path | None:
    """Hard-link one existing target into the transaction backup."""
    if not target.exists():
        return None
    relative = target.relative_to(repo_root)
    backup = backup_root / relative
    backup.parent.mkdir(parents=True, exist_ok=True)
    backup.unlink(missing_ok=True)
    os.link(target, backup)
    return backup


def optional_file_hash(path: Path) -> str | None:
    """Hash a regular file, returning None only when the path is absent."""
    if not path.exists():
        return None
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"transaction path is not a regular file: {path}")
    return sha256_file(path)


def transaction_digest(transaction: dict) -> str:
    """Seal the immutable repository and operation fields of a transaction."""
    sealed = {
        "schema_version": transaction["schema_version"],
        "repo_root": transaction["repo_root"],
        "operations": transaction["operations"],
    }
    payload = json.dumps(
        sealed, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_transaction(path: Path, transaction: dict) -> None:
    """Atomically persist and flush one transaction state document."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(transaction, stream, indent=2)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def validate_transaction(transaction: dict, repo_root: Path) -> list[dict]:
    """Reject a transaction whose root, schema, or sealed plan has drifted."""
    if transaction.get("schema_version") != 1:
        raise RuntimeError("unsupported promotion transaction schema")
    if transaction.get("repo_root") != str(repo_root.resolve()):
        raise RuntimeError("promotion transaction repository root changed")
    if transaction.get("plan_sha256") != transaction_digest(transaction):
        raise RuntimeError("promotion transaction seal mismatch")
    operations = transaction.get("operations")
    if not isinstance(operations, list):
        raise RuntimeError("promotion transaction operations are invalid")
    return operations


def recover_transaction(repo_root: Path, evidence_dir: Path) -> dict:
    """Restore pre-state or finish cleanup from a sealed COMPLETE transaction."""
    repo_root = repo_root.resolve()
    active_path = evidence_dir / TRANSACTION_ACTIVE
    complete_path = evidence_dir / TRANSACTION_COMPLETE
    recovered_path = evidence_dir / TRANSACTION_RECOVERED
    backup_root = evidence_dir / TRANSACTION_BACKUP
    if complete_path.is_file():
        # COMPLETE is written only after caller-level convergence. A crash during
        # final cleanup therefore finishes forward instead of rolling back a
        # transaction whose accepted postconditions are already durably sealed.
        complete = json.loads(complete_path.read_text(encoding="utf-8"))
        operations = validate_transaction(complete, repo_root)
        if complete.get("state") != "complete":
            raise RuntimeError("completed promotion transaction state is invalid")
        if active_path.is_file():
            active = json.loads(active_path.read_text(encoding="utf-8"))
            validate_transaction(active, repo_root)
            if active["plan_sha256"] != complete["plan_sha256"]:
                raise RuntimeError(
                    "complete transaction seal does not match active plan")
        for operation in operations:
            target = repo_root / safe_relative(operation["target_relative"])
            if optional_file_hash(target) != operation["after_sha256"]:
                raise RuntimeError(
                    f"completed transaction hash mismatch: {target}")
        active_path.unlink(missing_ok=True)
        shutil.rmtree(backup_root, ignore_errors=True)
        return complete
    if not active_path.is_file():
        if recovered_path.is_file():
            return json.loads(recovered_path.read_text(encoding="utf-8"))
        raise FileNotFoundError("no unfinished promotion transaction to recover")

    transaction = json.loads(active_path.read_text(encoding="utf-8"))
    operations = validate_transaction(transaction, repo_root)
    if recovered_path.is_file():
        recovered = json.loads(recovered_path.read_text(encoding="utf-8"))
        validate_transaction(recovered, repo_root)
        if recovered["plan_sha256"] != transaction["plan_sha256"]:
            raise RuntimeError("recovered transaction seal does not match active plan")
    elif transaction.get("state") not in {"preparing", "applying", "committed"}:
        raise RuntimeError("promotion transaction state cannot be recovered")
    elif transaction["state"] != "preparing":
        # Verify the complete backup set before changing any canonical path.
        for operation in operations:
            before = operation["before_sha256"]
            if before is None:
                continue
            backup = backup_root / safe_relative(operation["target_relative"])
            if optional_file_hash(backup) != before:
                raise RuntimeError(f"missing or changed transaction backup: {backup}")

        # The sealed operation list is the durable touched-set ledger. Restoring
        # it in reverse also removes targets that did not exist before apply.
        for operation in reversed(operations):
            relative = safe_relative(operation["target_relative"])
            target = repo_root / relative
            target.unlink(missing_ok=True)
            if operation["before_sha256"] is not None:
                backup = backup_root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                os.link(backup, target)

    # Both a never-mutated preparing state and a restored applying state must
    # exactly equal every sealed precondition before recovery is published.
    for operation in operations:
        target = repo_root / safe_relative(operation["target_relative"])
        if optional_file_hash(target) != operation["before_sha256"]:
            raise RuntimeError(f"transaction recovery hash mismatch: {target}")
    recovered = dict(transaction)
    recovered["state"] = "recovered"
    write_transaction(recovered_path, recovered)
    active_path.unlink(missing_ok=True)
    shutil.rmtree(backup_root, ignore_errors=True)
    return recovered


def finalize_transaction(repo_root: Path, evidence_dir: Path) -> dict:
    """Seal a converged transaction complete and retire its rollback backup."""
    repo_root = repo_root.resolve()
    active_path = evidence_dir / TRANSACTION_ACTIVE
    complete_path = evidence_dir / TRANSACTION_COMPLETE
    backup_root = evidence_dir / TRANSACTION_BACKUP
    if complete_path.is_file():
        complete = json.loads(complete_path.read_text(encoding="utf-8"))
        operations = validate_transaction(complete, repo_root)
        if active_path.is_file():
            active = json.loads(active_path.read_text(encoding="utf-8"))
            validate_transaction(active, repo_root)
            if active["plan_sha256"] != complete["plan_sha256"]:
                raise RuntimeError(
                    "complete transaction seal does not match active plan")
    else:
        if not active_path.is_file():
            raise FileNotFoundError("no committed promotion transaction to finalize")
        complete = json.loads(active_path.read_text(encoding="utf-8"))
        operations = validate_transaction(complete, repo_root)
        if complete.get("state") != "committed":
            raise RuntimeError("promotion transaction is not committed")

    # Finalization occurs only after caller-level convergence, but it repeats
    # each operation-level postcondition before discarding rollback material.
    for operation in operations:
        target = repo_root / safe_relative(operation["target_relative"])
        if optional_file_hash(target) != operation["after_sha256"]:
            raise RuntimeError(f"transaction finalization hash mismatch: {target}")
    complete["state"] = "complete"
    write_transaction(complete_path, complete)
    active_path.unlink(missing_ok=True)
    shutil.rmtree(backup_root, ignore_errors=True)
    return complete


def apply_plan(rows, pruned, repo_root: Path, evidence_dir: Path) -> bool:
    """Apply a sealed plan and retain rollback material until convergence."""
    repo_root = repo_root.resolve()
    evidence_dir.mkdir(parents=True, exist_ok=True)
    active_path = evidence_dir / TRANSACTION_ACTIVE
    complete_path = evidence_dir / TRANSACTION_COMPLETE
    backup_root = evidence_dir / TRANSACTION_BACKUP
    if active_path.exists() or backup_root.exists():
        raise RuntimeError(
            "unfinished promotion transaction exists; recover it before apply")

    selected = [row for row in rows if row["classification"] == "promote"]
    selected.sort(key=lambda row: row["source_relative"].endswith("manifest.json"))
    writes = []
    for row in selected:
        source = Path(row["source_path"])
        target = Path(row["target_path"])
        before = row["target_sha256_before"] or None
        if sha256_file(source) != row["sha256"]:
            raise RuntimeError(f"stale promotion source plan: {source}")
        if optional_file_hash(target) != before:
            raise RuntimeError(f"stale promotion target plan: {target}")
        if row["canonical_status"] != "identical":
            writes.append({
                "kind": "write",
                "target_relative": target.resolve().relative_to(repo_root).as_posix(),
                "source_path": str(source),
                "before_sha256": before,
                "after_sha256": row["sha256"],
            })

    prunes = []
    for relative in pruned:
        target = repo_root / "data" / relative
        before = optional_file_hash(target)
        if before is None:
            raise RuntimeError(f"stale promotion prune plan: {target}")
        prunes.append({
            "kind": "prune",
            "target_relative": target.relative_to(repo_root).as_posix(),
            "source_path": "",
            "before_sha256": before,
            "after_sha256": None,
        })

    # Stale payloads disappear before replacement manifests publish the new
    # references. The operation order itself is sealed for recovery.
    payload_writes = [
        operation for operation in writes
        if not operation["target_relative"].endswith("manifest.json")
    ]
    manifest_writes = [
        operation for operation in writes
        if operation["target_relative"].endswith("manifest.json")
    ]
    operations = payload_writes + prunes + manifest_writes
    if not operations:
        # A fresh already-converged tree still needs a durable transaction seal
        # for downstream gates; repeated applies preserve the existing seal.
        recovered_path = evidence_dir / TRANSACTION_RECOVERED
        if not complete_path.exists() and not recovered_path.exists():
            transaction = {
                "schema_version": 1,
                "state": "complete",
                "repo_root": str(repo_root),
                "operations": [],
            }
            transaction["plan_sha256"] = transaction_digest(transaction)
            write_transaction(complete_path, transaction)
        return False
    if complete_path.exists() or (
            evidence_dir / TRANSACTION_RECOVERED).exists():
        raise RuntimeError(
            "promotion evidence is already sealed; use a fresh evidence directory")

    transaction = {
        "schema_version": 1,
        "state": "preparing",
        "repo_root": str(repo_root),
        "operations": operations,
    }
    transaction["plan_sha256"] = transaction_digest(transaction)
    write_transaction(active_path, transaction)
    try:
        # Seal every pre-existing byte before the first mutation. A process loss
        # after this point leaves both the immutable plan and complete backups.
        for operation in operations:
            if operation["before_sha256"] is None:
                continue
            target = repo_root / safe_relative(operation["target_relative"])
            backup = backup_file(target, repo_root, backup_root)
            if optional_file_hash(backup) != operation["before_sha256"]:
                raise RuntimeError(f"transaction backup hash mismatch: {target}")
        transaction["state"] = "applying"
        write_transaction(active_path, transaction)

        for operation in operations:
            target = repo_root / safe_relative(operation["target_relative"])
            if operation["kind"] == "write":
                clone_or_copy(
                    Path(operation["source_path"]),
                    target,
                    operation["after_sha256"],
                )
            else:
                target.unlink()

        # Rehash every selected target and require every approved prune absent.
        for row in selected:
            target = Path(row["target_path"])
            if sha256_file(target) != row["sha256"]:
                raise RuntimeError(f"post-promotion hash mismatch: {target}")
        if any((repo_root / "data" / relative).exists() for relative in pruned):
            raise RuntimeError("one or more superseded artifacts survived pruning")
        transaction["state"] = "committed"
        write_transaction(active_path, transaction)
    except Exception:
        try:
            recover_transaction(repo_root, evidence_dir)
        except Exception as recovery_error:
            raise RuntimeError(
                "promotion failed and automatic recovery also failed") \
                from recovery_error
        raise
    return True


def parse_args() -> argparse.Namespace:
    """Parse the explicit candidate, repository, and evidence roots."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--repo-root", default=Path.cwd(), type=Path)
    parser.add_argument("--evidence-dir", required=True, type=Path)
    parser.add_argument(
        "--families",
        nargs="+",
        default=list(TARGET_FAMILIES),
        help="manifest families to promote (default: seasonal-snow families)",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true")
    mode.add_argument(
        "--recover",
        action="store_true",
        help="reconcile an interrupted sealed transaction instead of planning",
    )
    return parser.parse_args()


def main() -> None:
    """Run inventory-only dry mode or the bounded promotion transaction."""
    args = parse_args()
    candidate_root = args.candidate_root.resolve()
    repo_root = args.repo_root.resolve()
    evidence_dir = args.evidence_dir.resolve()
    if args.recover:
        recovered = recover_transaction(repo_root, evidence_dir)
        print(json.dumps(recovered, sort_keys=True))
        return
    rows, _, pruned = build_inventory(candidate_root, repo_root, args.families)
    summary = write_evidence(
        rows,
        pruned,
        evidence_dir,
        "apply" if args.apply else "dry_run",
        args.families,
    )
    marker = evidence_dir / ("COMPLETE.json" if args.apply else "DRY_RUN_COMPLETE.json")
    if args.apply:
        transaction_started = apply_plan(rows, pruned, repo_root, evidence_dir)
        # Recompute after mutation so completion proves the selected tree is an
        # idempotent exact match and the second dry run has zero changed/new rows.
        after_rows, _, after_pruned = build_inventory(
            candidate_root, repo_root, args.families)
        after = write_evidence(
            after_rows,
            after_pruned,
            evidence_dir,
            "complete",
            args.families,
        )
        expected_statuses = {
            "identical": after["classification_counts"]["promote"],
        }
        if after["canonical_status_counts"] != expected_statuses:
            raise RuntimeError("promotion did not converge to an all-identical selected set")
        if after["prune_count"] != 0:
            raise RuntimeError("promotion did not converge to a zero-prune second plan")
        if transaction_started:
            finalize_transaction(repo_root, evidence_dir)
        summary = after
    marker.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
