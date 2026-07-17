"""Validate the rendered snow and firn QA report using only the standard library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlparse


FAMILIES = {"promice", "esm_snowmip", "laugh_tests"}
FIRN_FAMILIES = {"retmip", "imau", "sumup", "research_site"}
SEASONAL_CASE_COUNT = 62
FIRN_CASE_COUNT = 56
COVERAGE_CASE_COUNT = 61


class ReportParser(HTMLParser):
    """Collect linked images, anchors, headings, and visible report text."""

    def __init__(self) -> None:
        super().__init__()
        self.images: list[tuple[str, str]] = []
        self.linked_images: list[tuple[str, str, str]] = []
        self.links: list[str] = []
        self.ids: list[str] = []
        self.toc_anchors: list[str] = []
        self.h3_ids: list[str] = []
        self.text: list[str] = []
        self.in_toc = False
        self.current_href = ""

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        """Record only attributes needed for report-integrity checks."""
        values = dict(attrs)
        if values.get("id"):
            self.ids.append(values["id"] or "")
        if tag == "nav" and values.get("id") == "TOC":
            self.in_toc = True
        if tag == "a" and values.get("href"):
            self.current_href = values["href"] or ""
            self.links.append(self.current_href)
            if self.in_toc:
                self.toc_anchors.append(self.current_href)
        if tag == "img" and values.get("src"):
            source = values["src"] or ""
            loading = values.get("loading") or ""
            self.images.append((source, loading))
            self.linked_images.append((self.current_href, source, loading))
        # Quarto 1.9 owns level-three IDs on sections; retain the direct-H3
        # fallback for older valid renders without double-counting either shape.
        elif tag == "section" and values.get("id"):
            classes = (values.get("class") or "").split()
            if "level3" in classes:
                self.h3_ids.append(values["id"] or "")
        elif tag == "h3":
            heading_id = values.get("id") or values.get("data-anchor-id") or ""
            if heading_id and (not self.h3_ids or self.h3_ids[-1] != heading_id):
                self.h3_ids.append(heading_id)

    def handle_endtag(self, tag: str) -> None:
        """Leave the site-level TOC when its navigation element closes."""
        if tag == "a":
            self.current_href = ""
        if self.in_toc and tag == "nav":
            self.in_toc = False

    def handle_data(self, data: str) -> None:
        """Keep visible text for QA timestamp and QMD-prose identity checks."""
        if data.strip():
            self.text.append(data)


def require(condition: bool, message: str) -> None:
    """Raise a validation error without relying on optimizable assertions."""
    if not condition:
        raise ValueError(message)


def read_nonempty(path: Path, label: str) -> bytes:
    """Read one required nonempty report input."""
    require(path.is_file(), f"missing {label}: {path}")
    content = path.read_bytes()
    require(bool(content), f"empty {label}: {path}")
    return content


def digest(content: bytes) -> str:
    """Return the SHA-256 identity of one validated input."""
    return hashlib.sha256(content).hexdigest()


def normalized_text(text: str) -> str:
    """Collapse whitespace for comparisons across Markdown and rendered HTML."""
    punctuation = str.maketrans({"\u2018": "'", "\u2019": "'", "\u201c": '"', "\u201d": '"'})
    value = " ".join(text.translate(punctuation).replace("\u00a0", " ").split())
    return re.sub(r"\s+([,.;:!?])", r"\1", value)


def qmd_prose_blocks(qmd_text: str) -> list[str]:
    """Extract human prose that must survive the Quarto render."""
    lines = qmd_text.splitlines()
    in_frontmatter = bool(lines and lines[0].strip() == "---")
    in_fence = False
    blocks: list[str] = []
    current: list[str] = []

    # Ignore YAML, executable/raw blocks, and the include directive itself.
    for line in lines[1:] if in_frontmatter else lines:
        stripped = line.strip()
        if in_frontmatter:
            if stripped == "---":
                in_frontmatter = False
            continue
        if stripped.startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence or stripped.startswith("{{< include "):
            continue
        if not stripped:
            if current:
                blocks.append(" ".join(current))
                current = []
            continue
        current.append(stripped.lstrip("# "))
    if current:
        blocks.append(" ".join(current))

    # Reduce the small inline-Markdown subset used by the tracked report prose.
    plain: list[str] = []
    for block in blocks:
        block = re.sub(r"\[([^]]+)\]\([^)]+\)", r"\1", block)
        block = re.sub(r"`([^`]*)`", r"\1", block)
        # Remove the paired emphasis delimiter used by the QMD, but preserve a
        # literal multiplication operator that originated inside inline code.
        block = block.replace("**", "")
        value = normalized_text(block)
        if value:
            plain.append(value)
    return plain


def load_qa(
    qa_path: Path,
    expected_families: set[str],
    expected_case_count: int,
    label: str,
) -> tuple[dict, datetime, bytes]:
    """Load the QA timestamp and case count used to authorize the report."""
    content = read_nonempty(qa_path, f"{label} artifact-QA JSON")
    try:
        qa = json.loads(content)
        generated_at = datetime.fromisoformat(
            str(qa["generated_at"]).replace("Z", "+00:00")
        )
        families = qa["families"]
        case_count = qa["summary"]["case_count"]
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid artifact-QA JSON: {qa_path}") from error
    require(
        generated_at.tzinfo is not None and generated_at.utcoffset() is not None,
        "artifact-QA generated_at must include a timezone",
    )
    require(
        isinstance(families, list)
        and len(families) == len(expected_families)
        and set(families) == expected_families,
        f"{label} artifact-QA families do not match the accepted scope",
    )
    require(
        type(case_count) is int and case_count == expected_case_count,
        f"{label} artifact-QA case_count must equal {expected_case_count}",
    )
    return qa, generated_at, content


def load_firn_index(index_path: Path, preview: Path) -> tuple[bytes, set[Path]]:
    """Validate the exact canonical firn ledger/tree and return its paths."""
    content = read_nonempty(index_path, "firn figure index")
    try:
        reader = csv.DictReader(content.decode("utf-8").splitlines())
        rows = list(reader)
    except UnicodeDecodeError as error:
        raise ValueError(f"invalid firn figure index: {index_path}") from error
    required = {"dataset_family", "case_id", "figure_file"}
    require(
        reader.fieldnames is not None and required.issubset(reader.fieldnames),
        "firn figure index lacks required columns",
    )
    require(bool(rows), "firn figure index contains no rows")

    root = (preview / "figures/firn").resolve()
    paths: list[Path] = []
    cases: dict[str, set[str]] = {family: set() for family in FIRN_FAMILIES}
    for row in rows:
        family = row.get("dataset_family", "")
        case = row.get("case_id", "")
        filename = Path(row.get("figure_file", ""))
        require(
            family in FIRN_FAMILIES and bool(case) and filename.is_absolute(),
            "firn figure index contains an invalid family, case, or path",
        )
        resolved = filename.resolve()
        require(
            resolved.is_relative_to(root),
            "firn figure index contains a path outside figures/firn",
        )
        relative = resolved.relative_to(root)
        require(
            len(relative.parts) >= 3 and relative.parts[:2] == (family, case),
            "firn figure index path disagrees with its family or case",
        )
        paths.append(resolved)
        cases[family].add(case)
    tree = {
        path.resolve()
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() == ".png"
    }
    require(len(paths) == len(set(paths)), "firn figure index contains duplicate paths")
    require(set(paths) == tree, "firn figure index does not exactly match the PNG tree")
    require(
        sum(len(values) for values in cases.values()) == FIRN_CASE_COUNT
        and all(cases.values()),
        "firn figure index case map must contain exactly 4 families and 56 cases",
    )
    return content, set(paths)


def validate_qmd_include(qmd: Path, qmd_text: str, generated: Path) -> None:
    """Require one QMD include that resolves to the checked generated Markdown."""
    matches = re.findall(r"\{\{<\s*include\s+([^>]+?)\s*>\}\}", qmd_text)
    require(len(matches) == 1, "QMD must contain exactly one include directive")
    include = matches[0].strip().strip("\"").strip("'")
    target = (qmd.parent / include).resolve()
    require(target == generated.resolve(), "QMD include does not resolve to generated Markdown")


def validate_report(
    report: Path,
    generated: Path,
    preview: Path,
    qa_path: Path,
    qmd: Path,
) -> dict[str, int | str]:
    """Validate one rendered report against its exact QA, QMD, and figure inputs."""
    preview = preview.resolve()
    report_root = (preview / "report").resolve()
    qa_root = (preview / "qa").resolve()
    firn_qa_path = qa_root / "firn/artifact_qa.json"
    firn_index_path = qa_root / "firn_figure_index.csv"

    # Keep every report-side input in its canonical preview subtree.
    require(report.resolve().is_relative_to(report_root), "report is outside preview/report")
    require(
        generated.resolve().is_relative_to(report_root),
        "generated include is outside preview/report",
    )
    require(qa_path.resolve().is_relative_to(qa_root), "QA JSON is outside preview/qa")
    require(
        firn_qa_path.resolve().is_relative_to(qa_root),
        "firn QA JSON is outside preview/qa",
    )
    require(
        firn_index_path.resolve().is_relative_to(qa_root),
        "firn figure index is outside preview/qa",
    )

    # Read the complete evidence set only after its paths pass containment.
    report_content = read_nonempty(report, "rendered report")
    generated_content = read_nonempty(generated, "generated include")
    qmd_content = read_nonempty(qmd, "QMD source")
    qa, generated_at, qa_content = load_qa(
        qa_path, FAMILIES, SEASONAL_CASE_COUNT, "seasonal"
    )
    firn_qa, firn_generated_at, firn_qa_content = load_qa(
        firn_qa_path, FIRN_FAMILIES, FIRN_CASE_COUNT, "firn"
    )
    firn_index_content, firn_index_paths = load_firn_index(
        firn_index_path, preview
    )
    report_text = report_content.decode("utf-8")
    generated_text = generated_content.decode("utf-8")
    qmd_text = qmd_content.decode("utf-8")
    validate_qmd_include(qmd, qmd_text, generated)

    # The render must postdate every textual input and preserve the audited time.
    qa_time = generated_at.timestamp()
    firn_qa_time = firn_generated_at.timestamp()
    require(
        generated.stat().st_mtime
        >= max(
            qa_path.stat().st_mtime,
            firn_qa_path.stat().st_mtime,
            firn_index_path.stat().st_mtime,
        ),
        "generated include predates QA JSON or firn figure index",
    )
    newest_input = max(
        generated.stat().st_mtime,
        qa_path.stat().st_mtime,
        firn_qa_path.stat().st_mtime,
        firn_index_path.stat().st_mtime,
        qmd.stat().st_mtime,
    )
    require(report.stat().st_mtime >= newest_input, "rendered report predates an input")
    generated_at_text = str(qa["generated_at"])
    firn_generated_at_text = str(firn_qa["generated_at"])
    timestamp_count = 1 if generated_at_text != firn_generated_at_text else 2
    require(
        generated_text.count(generated_at_text) == timestamp_count,
        "generated include does not identify the exact QA timestamp once",
    )
    require(
        generated_text.count(firn_generated_at_text) == timestamp_count,
        "generated include does not identify the exact firn QA timestamp once",
    )

    # Parse both representations so their ordered linked-image identities agree.
    parser = ReportParser()
    parser.feed(report_text)
    generated_parser = ReportParser()
    generated_parser.feed(generated_text)
    rendered_text = normalized_text(" ".join(parser.text))
    require(generated_at_text in rendered_text, "rendered report lacks the QA timestamp")
    require(
        firn_generated_at_text in rendered_text,
        "rendered report lacks the firn QA timestamp",
    )
    for block in qmd_prose_blocks(qmd_text):
        require(block in rendered_text, "rendered report lacks tracked QMD prose")

    generated_lines = generated_text.splitlines()
    expected_sites = sum(line.startswith("### ") for line in generated_lines)
    expected_images = generated_parser.images
    expected_figures = len(expected_images)
    sources = [source for source, _ in parser.images]
    require(
        expected_sites == SEASONAL_CASE_COUNT + FIRN_CASE_COUNT,
        "generated site headings do not match the 118-case QA scope",
    )
    require(
        len(parser.h3_ids) == expected_sites,
        f"HTML site headings {len(parser.h3_ids)} != {expected_sites}",
    )
    require(len(set(parser.h3_ids)) == len(parser.h3_ids), "duplicate site heading ids")
    require(
        parser.images == expected_images,
        "HTML image sources or lazy-loading attributes do not match generated include",
    )
    require(
        parser.linked_images == generated_parser.linked_images,
        "HTML linked-image identities do not match generated include",
    )
    require(len(sources) == expected_figures, "HTML image count does not match generated include")
    require(len(set(sources)) == len(sources), "duplicate image links")
    require(parser.in_toc is False, "unclosed nav#TOC")

    # Every PROMICE and ESM-SnowMIP case links forward to one unique table in
    # the final appendix; Quarto must preserve the exact fragment graph.
    coverage_heading = "## Audited per-source channel coverage appendix"
    require(
        generated_text.count(coverage_heading) == 1,
        "generated include must contain one coverage appendix",
    )
    appendix_start = generated_text.index(coverage_heading)
    require(
        appendix_start > generated_text.rfind('loading="lazy"'),
        "coverage appendix must follow every figure",
    )
    generated_coverage_links = [
        unquote(urlparse(href).fragment)
        for href in generated_parser.links
        if not urlparse(href).path
        and unquote(urlparse(href).fragment).startswith("coverage-")
    ]
    generated_coverage_ids = [
        anchor for anchor in generated_parser.ids if anchor.startswith("coverage-")
    ]
    rendered_coverage_links = [
        unquote(urlparse(href).fragment)
        for href in parser.links
        if not urlparse(href).path
        and unquote(urlparse(href).fragment).startswith("coverage-")
    ]
    rendered_coverage_ids = [
        anchor for anchor in parser.ids if anchor.startswith("coverage-")
    ]
    require(
        len(generated_coverage_links) == COVERAGE_CASE_COUNT,
        "generated coverage-link count does not match the 61-case table scope",
    )
    require(
        generated_coverage_links == generated_coverage_ids,
        "generated coverage links do not match final appendix targets",
    )
    require(
        rendered_coverage_links == generated_coverage_links
        and rendered_coverage_ids == generated_coverage_ids,
        "rendered coverage links or targets do not match generated include",
    )

    # Resolve only canonical lazy links to current nonempty post-QA figure files.
    broken: list[str] = []
    stale: list[str] = []
    linked_firn: set[Path] = set()
    figures_root = (preview / "figures").resolve()
    for source, loading in parser.images:
        parsed = urlparse(source)
        decoded = unquote(parsed.path)
        parts = PurePosixPath(decoded).parts
        seasonal_path = (
            len(parts) == 5
            and parts[:2] == ("..", "figures")
            and parts[2] in FAMILIES
            and parts[3] not in {"", ".", ".."}
            and parts[4].lower().endswith(".png")
        )
        firn_path = (
            len(parts) == 6
            and parts[:3] == ("..", "figures", "firn")
            and parts[3] in FIRN_FAMILIES
            and parts[4] not in {"", ".", ".."}
            and parts[5].lower().endswith(".png")
        )
        canonical = (
            not parsed.scheme
            and not parsed.netloc
            and not parsed.query
            and not parsed.fragment
            and not Path(decoded).is_absolute()
            and (seasonal_path or firn_path)
            and loading == "lazy"
        )
        target = (report.parent / decoded).resolve()
        if (
            not canonical
            or not target.is_relative_to(figures_root)
            or not target.is_file()
            or target.stat().st_size == 0
        ):
            broken.append(source)
        elif target.stat().st_mtime < (firn_qa_time if firn_path else qa_time):
            stale.append(source)
        elif firn_path:
            linked_firn.add(target)
    require(
        all(href == source for href, source, _ in parser.linked_images),
        "linked figure anchors must equal their image sources",
    )
    require(not broken, f"broken or noncanonical image links: {broken[:5]}")
    require(not stale, f"figure files predate artifact QA: {stale[:5]}")
    require(
        linked_firn == firn_index_paths,
        "rendered firn image inventory does not match the figure index",
    )

    # Count only links inside a real Quarto site-level TOC.
    toc_targets = {href[1:] for href in parser.toc_anchors if href.startswith("#")}
    missing_toc = sorted(set(parser.h3_ids) - toc_targets)
    require(bool(parser.toc_anchors), "missing nav#TOC")
    require(not missing_toc, f"site headings missing from TOC: {missing_toc[:5]}")

    # Persist exact input/output identities with the structural check counts.
    return {
        "qa_generated_at": generated_at_text,
        "qa_sha256": digest(qa_content),
        "firn_qa_generated_at": firn_generated_at_text,
        "firn_qa_sha256": digest(firn_qa_content),
        "firn_figure_index_sha256": digest(firn_index_content),
        "generated_sha256": digest(generated_content),
        "qmd_sha256": digest(qmd_content),
        "report_sha256": digest(report_content),
        "site_heading_count": expected_sites,
        "image_count": expected_figures,
        "unique_image_count": len(set(sources)),
        "broken_image_count": len(broken),
        "stale_image_count": len(stale),
        "missing_toc_count": len(missing_toc),
        "embedded_image_count": sum(source.startswith("data:") for source in sources),
        "coverage_table_link_count": len(generated_coverage_links),
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse canonical defaults plus isolated-preview overrides."""
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview-root", type=Path, default=repo_root / "data/preview")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--generated", type=Path)
    parser.add_argument("--qa", type=Path)
    parser.add_argument("--qmd", type=Path, default=repo_root / "report/snow-artifact-qa.qmd")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)

    # Derive dependent defaults from an overridden preview root.
    if args.report is None:
        args.report = args.preview_root / "report/snow-artifact-qa.html"
    if args.generated is None:
        args.generated = args.preview_root / "report/figures.generated.md"
    if args.qa is None:
        args.qa = args.preview_root / "qa/artifact_qa.json"
    if args.output is None:
        args.output = args.report.parent / "report_checks.json"
    return args


def main(argv: list[str] | None = None) -> None:
    """Validate the report and write a compact reproducible evidence ledger."""
    args = parse_args(argv)
    result = validate_report(
        args.report,
        args.generated,
        args.preview_root,
        args.qa,
        args.qmd,
    )

    # Write evidence only after every fail-closed validation succeeds.
    require(
        args.output.resolve().parent == args.report.resolve().parent,
        "checker output must be beside the rendered report",
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
