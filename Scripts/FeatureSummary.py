#!/usr/bin/env python3
"""
genbank_plastome_tool.py

Unified command-line tool for GenBank/plastome analysis.

Goals
- Preserve legacy outputs where possible
- Merge duplicate parsing/stat functions
- Support single file, multiple files, directory, recursive scanning
- Export CSV/TSV
- Optional reports:
  * feature summary
  * plastome stats
  * gene position/length table
  * exon/intron tables

Dependencies:
    pip install biopython

Recommended:
    Python 3.10+
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqUtils import GC123

GENBANK_SUFFIXES = {".gb", ".gbk", ".gbf", ".genbank"}

TARGET_FEATURES = ["tRNA", "intron", "CDS", "gene", "rRNA", "exon"]


@dataclass
class Region:
    label: str
    start: int
    end: int
    length: int


@dataclass
class RegionSet:
    lsc: Region | None = None
    ssc: Region | None = None
    ira: Region | None = None
    irb: Region | None = None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    prog_stem = Path(sys.argv[0]).stem.lower()

    parser = argparse.ArgumentParser(
        prog=Path(sys.argv[0]).name,
        description="Unified GenBank/plastome analysis tool",
    )

    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="One or more GenBank files and/or directories",
    )
    parser.add_argument(
        "-r", "--recursive",
        action="store_true",
        help="Recursively scan directories",
    )
    parser.add_argument(
        "--glob",
        default="*",
        help="Optional file glob pattern inside directories (default: *)",
    )

    parser.add_argument(
        "--compat",
        choices=[
            "auto", "feature-summary", "plastome-stats",
            "gene-position", "intron-exon", "all"
        ],
        default="auto",
        help="Legacy compatibility preset",
    )

    # Toggle pairs
    parser.add_argument("--infer-introns", dest="infer_introns", action="store_true", default=True)
    parser.add_argument("--no-infer-introns", dest="infer_introns", action="store_false")

    parser.add_argument("--regions", dest="compute_regions", action="store_true", default=True)
    parser.add_argument("--no-regions", dest="compute_regions", action="store_false")

    parser.add_argument("--cds-stats", dest="compute_cds_stats", action="store_true", default=True)
    parser.add_argument("--no-cds-stats", dest="compute_cds_stats", action="store_false")

    parser.add_argument("--codon-at", dest="compute_codon_at", action="store_true", default=True)
    parser.add_argument("--no-codon-at", dest="compute_codon_at", action="store_false")

    parser.add_argument("--gene-position-table", action="store_true")
    parser.add_argument("--exon-intron-tables", action="store_true")

    parser.add_argument("-o", "--output", help="Primary output file or prefix")
    parser.add_argument("--feature-summary-out")
    parser.add_argument("--plastome-stats-out")
    parser.add_argument("--gene-position-out")
    parser.add_argument("--gene-intron-summary-out")
    parser.add_argument("--exon-table-out")
    parser.add_argument("--intron-table-out")

    parser.add_argument(
        "--export",
        choices=["csv", "tsv", "both"],
        default="tsv",
        help="Export format. Default: tsv",
    )
    parser.add_argument("--encoding", default="utf-8")
    parser.add_argument("--strict", action="store_true", help="Stop on first file error")
    parser.add_argument("--skip-errors", action="store_true", help="Skip invalid files/records")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args(argv)

    if args.strict and args.skip_errors:
        parser.error("--strict and --skip-errors are mutually exclusive")

    # Unified-script behavior:
    # The script name must not silently determine which analyses are run.
    # With --compat auto (default), run all reports.
    # Legacy/single-report behavior remains available through an explicit
    # --compat feature-summary|plastome-stats|gene-position|intron-exon option.
    if args.compat == "auto":
        args.compat = "all"

    # Presets
    if args.compat == "feature-summary":
        args.compute_regions = False
        args.compute_cds_stats = False
        args.compute_codon_at = False
        args.gene_position_table = False
        args.exon_intron_tables = False
    elif args.compat == "plastome-stats":
        args.compute_regions = True
        args.compute_cds_stats = True
        args.compute_codon_at = True
    elif args.compat == "gene-position":
        args.gene_position_table = True
        args.compute_regions = False
        args.compute_cds_stats = False
        args.compute_codon_at = False
    elif args.compat == "intron-exon":
        args.exon_intron_tables = True

    return args


def collect_input_paths(inputs: Sequence[str], recursive: bool, pattern: str) -> list[Path]:
    results: list[Path] = []

    for item in inputs:
        p = Path(item)
        if p.is_file():
            results.append(p)
        elif p.is_dir():
            iterator = p.rglob(pattern) if recursive else p.glob(pattern)
            for child in iterator:
                if child.is_file() and child.suffix.lower() in GENBANK_SUFFIXES:
                    results.append(child)
        else:
            # allow shell-expanded patterns to have failed on Windows; try glob from cwd
            for child in Path(".").glob(item):
                if child.is_file() and child.suffix.lower() in GENBANK_SUFFIXES:
                    results.append(child)

    # Stable order for reproducible outputs/tests
    unique_sorted = sorted({x.resolve() for x in results})
    return [Path(x) for x in unique_sorted]


def load_records(path: Path):
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found: {path}")
    return records


def qualifier_first(feature, keys=("gene", "locus_tag", "product", "note")) -> str:
    q = getattr(feature, "qualifiers", {}) or {}
    for k in keys:
        if k in q and q[k]:
            return str(q[k][0])
    return "unknown"


def strict_gc_percent(seq) -> float:
    """Keep legacy behavior: only A/T/G/C counted; ambiguous bases ignored."""
    s = str(seq).upper()
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    total = a + t + g + c
    return 0.0 if total == 0 else (g + c) / total * 100.0


def feature_length(feature) -> int:
    try:
        return int(len(feature.location))
    except Exception:
        start = int(feature.location.start)
        end = int(feature.location.end)
        return max(0, end - start)


def summarize_feature_counts(records, infer_introns: bool = True) -> dict[str, int | str]:
    counter: dict[str, int] = {}
    total_seq = []

    for record in records:
        total_seq.append(str(record.seq))
        for feature in record.features:
            counter[feature.type] = counter.get(feature.type, 0) + 1

    if infer_introns and counter.get("intron", 0) == 0:
        counter["intron"] = len(infer_introns_from_join_features(records))

    gc = strict_gc_percent("".join(total_seq))
    out: dict[str, int | str] = {}
    for ft in TARGET_FEATURES:
        out[ft] = counter.get(ft, 0)
    out["GC Content"] = f"{gc:.2f}%"
    return out


def infer_introns_from_join_features(records):
    inferred = set()

    for record in records:
        for feature in record.features:
            if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
                continue
            location = feature.location
            if not isinstance(location, CompoundLocation):
                continue

            parts = sorted(location.parts, key=lambda x: int(x.start))
            if len(parts) < 2:
                continue

            fname = qualifier_first(feature)
            strand = parts[0].strand

            for prev_part, next_part in zip(parts, parts[1:]):
                intron_start = int(prev_part.end)
                intron_end = int(next_part.start)
                if intron_end > intron_start:
                    inferred.add(
                        (record.id, fname, strand, intron_start, intron_end)
                    )
    return inferred


def feature_rows(records, source_file: Path):
    rows = []
    for record in records:
        for feature in record.features:
            rows.append({
                "file": source_file.name,
                "record_id": record.id,
                "feature_type": feature.type,
                "name": qualifier_first(feature),
                "start": int(feature.location.start) + 1,  # 1-based user-facing
                "end": int(feature.location.end),
                "strand": feature.location.strand,
                "length": feature_length(feature),
            })
    return rows


def detect_regions(record) -> RegionSet:
    """
    Annotation-first region detector.

    Priority:
    1) explicit labels/notes for LSC/SSC/IRa/IRb
    2) repeat_region heuristics for IRs
    3) infer LSC/SSC from two IRs if possible
    """
    rs = RegionSet()
    repeat_candidates = []

    for feature in record.features:
        label = " ".join(
            [qualifier_first(feature, ("gene", "product", "note"))]
        ).lower()
        ftype = feature.type.lower()
        start = int(feature.location.start)
        end = int(feature.location.end)
        region = Region(label=feature.type, start=start, end=end, length=max(0, end - start))

        text = f"{ftype} {label}"
        if "large single copy" in text or " lsc" in text:
            rs.lsc = Region("LSC", start, end, region.length)
        elif "small single copy" in text or " ssc" in text:
            rs.ssc = Region("SSC", start, end, region.length)
        elif "ira" in text or "ir a" in text:
            rs.ira = Region("IRa", start, end, region.length)
        elif "irb" in text or "ir b" in text:
            rs.irb = Region("IRb", start, end, region.length)
        elif ftype == "repeat_region":
            repeat_candidates.append(region)

    if (rs.ira is None or rs.irb is None) and len(repeat_candidates) >= 2:
        repeat_candidates.sort(key=lambda r: r.length, reverse=True)
        a, b = sorted(repeat_candidates[:2], key=lambda r: r.start)
        rs.irb = rs.irb or Region("IRb", a.start, a.end, a.length)
        rs.ira = rs.ira or Region("IRa", b.start, b.end, b.length)

    if rs.lsc is None or rs.ssc is None:
        if rs.ira and rs.irb:
            inferred_lsc, inferred_ssc = infer_sc_regions_from_irs(
                len(record.seq), rs.irb, rs.ira
            )
            rs.lsc = rs.lsc or inferred_lsc
            rs.ssc = rs.ssc or inferred_ssc

    return rs


def infer_sc_regions_from_irs(genome_len: int, ir1: Region, ir2: Region):
    # assume non-overlapping IRs on a circular genome
    left, right = sorted([ir1, ir2], key=lambda x: x.start)
    gap1_start, gap1_end = left.end, right.start
    gap2_start, gap2_end = right.end, left.start

    len1 = max(0, gap1_end - gap1_start)
    len2 = max(0, genome_len - gap2_start + gap2_end)

    region1 = Region("SC1", gap1_start, gap1_end, len1)
    region2 = Region("SC2", gap2_start, gap2_end, len2)

    if region1.length >= region2.length:
        lsc = Region("LSC", region1.start, region1.end, region1.length)
        ssc = Region("SSC", region2.start, region2.end, region2.length)
    else:
        lsc = Region("LSC", region2.start, region2.end, region2.length)
        ssc = Region("SSC", region1.start, region1.end, region1.length)
    return lsc, ssc


def extract_circular_region(seq: str, start: int, end: int) -> str:
    if end >= start:
        return seq[start:end]
    return seq[start:] + seq[:end]


def compute_plastome_stats(records, source_file: Path,
                           compute_regions: bool = True,
                           compute_cds_stats: bool = True,
                           compute_codon_at: bool = True):
    rows = []

    for record in records:
        genome_seq = str(record.seq)
        total_gc = strict_gc_percent(genome_seq)
        row = {
            "file": source_file.name,
            "record_id": record.id,
            "genome_bp": len(genome_seq),
            "gc_total_pct": f"{total_gc:.2f}",
        }

        # basic gene class counts
        protein = sum(1 for f in record.features if f.type == "CDS")
        trna = sum(1 for f in record.features if f.type == "tRNA")
        rrna = sum(1 for f in record.features if f.type == "rRNA")
        row["protein_coding_genes"] = protein
        row["tRNA"] = trna
        row["rRNA"] = rrna
        row["total_genes"] = protein + trna + rrna

        if compute_regions:
            rs = detect_regions(record)
            for key, region in [
                ("lsc", rs.lsc), ("ssc", rs.ssc), ("ira", rs.ira), ("irb", rs.irb)
            ]:
                row[f"{key}_bp"] = region.length if region else ""
                if region:
                    region_seq = extract_circular_region(genome_seq, region.start, region.end)
                    row[f"{key}_gc_pct"] = f"{strict_gc_percent(region_seq):.2f}"
                else:
                    row[f"{key}_gc_pct"] = ""

            if rs.ira and rs.irb:
                row["ir_mean_bp"] = (rs.ira.length + rs.irb.length) / 2
            else:
                row["ir_mean_bp"] = ""

        if compute_cds_stats:
            cds_total_bp = sum(feature_length(f) for f in record.features if f.type == "CDS")
            row["cds_total_bp"] = cds_total_bp
            row["cds_pct_of_genome"] = f"{(cds_total_bp / len(genome_seq) * 100):.2f}" if genome_seq else ""

        if compute_codon_at:
            at1, at2, at3 = compute_codon_at_percent(record)
            row["AT1_pct"] = f"{at1:.2f}" if at1 is not None else ""
            row["AT2_pct"] = f"{at2:.2f}" if at2 is not None else ""
            row["AT3_pct"] = f"{at3:.2f}" if at3 is not None else ""

        rows.append(row)

    return rows


def compute_codon_at_percent(record):
    parts = []

    for feature in record.features:
        if feature.type != "CDS":
            continue
        try:
            frag = str(feature.location.extract(record.seq)).upper()
        except Exception:
            continue
        # Keep only strict DNA letters to avoid GC123 ambiguity issue
        frag = "".join(x for x in frag if x in "ATGC")
        if len(frag) < 3:
            continue
        usable = len(frag) - (len(frag) % 3)
        if usable:
            parts.append(frag[:usable])

    cds_concat = "".join(parts)
    if len(cds_concat) < 3:
        return None, None, None

    _, gc1, gc2, gc3 = GC123(cds_concat)
    return 100 - gc1, 100 - gc2, 100 - gc3


def build_exon_intron_tables(records, source_file: Path, infer_introns: bool = True):
    gene_rows = []
    exon_rows = []
    intron_rows = []

    inferred_introns = infer_introns_from_join_features(records) if infer_introns else set()

    inferred_lookup = {}
    for record_id, gene_name, strand, start, end in inferred_introns:
        inferred_lookup.setdefault((record_id, gene_name), []).append((start, end, strand))

    for record in records:
        for feature in record.features:
            if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
                continue

            gene_name = qualifier_first(feature)
            exons = []
            if isinstance(feature.location, CompoundLocation):
                parts = sorted(feature.location.parts, key=lambda x: int(x.start))
                for i, part in enumerate(parts, start=1):
                    exon_rows.append({
                        "file": source_file.name,
                        "record_id": record.id,
                        "gene": gene_name,
                        "parent_type": feature.type,
                        "exon_no": i,
                        "start": int(part.start) + 1,
                        "end": int(part.end),
                        "strand": part.strand,
                        "length": int(part.end) - int(part.start),
                    })
                    exons.append((int(part.start), int(part.end)))
            else:
                exon_rows.append({
                    "file": source_file.name,
                    "record_id": record.id,
                    "gene": gene_name,
                    "parent_type": feature.type,
                    "exon_no": 1,
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                    "strand": feature.location.strand,
                    "length": feature_length(feature),
                })
                exons.append((int(feature.location.start), int(feature.location.end)))

            introns = []
            key = (record.id, gene_name)
            for idx, (s, e, strand) in enumerate(sorted(inferred_lookup.get(key, [])), start=1):
                intron_rows.append({
                    "file": source_file.name,
                    "record_id": record.id,
                    "gene": gene_name,
                    "parent_type": feature.type,
                    "intron_no": idx,
                    "start": s + 1,
                    "end": e,
                    "strand": strand,
                    "length": e - s,
                    "source": "inferred_join",
                })
                introns.append((s, e))

            gene_rows.append({
                "file": source_file.name,
                "record_id": record.id,
                "gene": gene_name,
                "parent_type": feature.type,
                "exon_count": len(exons),
                "intron_count": len(introns),
                "gene_start": int(feature.location.start) + 1,
                "gene_end": int(feature.location.end),
                "gene_length": feature_length(feature),
            })

    return gene_rows, exon_rows, intron_rows


def write_delimited(rows, output_path: Path, fmt: str, encoding: str = "utf-8"):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "csv":
        delimiter = ","
    elif fmt == "tsv":
        delimiter = "\t"
    else:
        raise ValueError(f"Unsupported export format: {fmt}")

    rows = list(rows)
    if not rows:
        return

    fieldnames = list(rows[0].keys())
    with open(output_path, "w", newline="", encoding=encoding) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def choose_output_path(primary: str | None, explicit: str | None, default_name: str, ext: str) -> Path:
    if explicit:
        return Path(explicit)
    if primary:
        p = Path(primary)
        if p.suffix:
            return p
        return p / f"{default_name}.{ext}"
    return Path(f"{default_name}.{ext}")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    paths = collect_input_paths(args.input, args.recursive, args.glob)

    if not paths:
        raise SystemExit("No GenBank input files found.")

    feature_rows_all = []
    plastome_rows_all = []
    genepos_rows_all = []
    gene_summary_all = []
    exon_rows_all = []
    intron_rows_all = []

    for path in paths:
        try:
            records = load_records(path)
        except Exception as exc:
            msg = f"[ERROR] {path}: {exc}"
            if args.strict:
                raise SystemExit(msg)
            if not args.quiet:
                print(msg, file=sys.stderr)
            if args.skip_errors:
                continue
            raise

        if args.compat in {"feature-summary", "all"}:
            fs = summarize_feature_counts(records, infer_introns=args.infer_introns)
            feature_rows_all.append({
                "source_file": path.name,
                "tRNA": fs["tRNA"],
                "intron": fs["intron"],
                "CDS": fs["CDS"],
                "gene": fs["gene"],
                "rRNA": fs["rRNA"],
                "exon": fs["exon"],
                "GC Content": fs["GC Content"],
            })

        if args.compat in {"plastome-stats", "all"} or args.compute_regions or args.compute_cds_stats or args.compute_codon_at:
            plastome_rows_all.extend(
                compute_plastome_stats(
                    records, path,
                    compute_regions=args.compute_regions,
                    compute_cds_stats=args.compute_cds_stats,
                    compute_codon_at=args.compute_codon_at,
                )
            )

        if args.compat in {"gene-position", "all"} or args.gene_position_table:
            genepos_rows_all.extend(feature_rows(records, path))

        if args.compat in {"intron-exon", "all"} or args.exon_intron_tables:
            gene_rows, exon_rows, intron_rows = build_exon_intron_tables(
                records, path, infer_introns=args.infer_introns
            )
            gene_summary_all.extend(gene_rows)
            exon_rows_all.extend(exon_rows)
            intron_rows_all.extend(intron_rows)

    formats = ["csv", "tsv"] if args.export == "both" else [args.export]

    for fmt in formats:
        if feature_rows_all:
            out = choose_output_path(args.output, args.feature_summary_out, "feature_summary", fmt)
            write_delimited(feature_rows_all, out, fmt, args.encoding)

        if plastome_rows_all:
            out = choose_output_path(args.output, args.plastome_stats_out, "plastome_stats", fmt)
            write_delimited(plastome_rows_all, out, fmt, args.encoding)

        if genepos_rows_all:
            out = choose_output_path(args.output, args.gene_position_out, "gene_position_length", fmt)
            write_delimited(genepos_rows_all, out, fmt, args.encoding)

        if gene_summary_all:
            out = choose_output_path(args.output, args.gene_intron_summary_out, "gene_intron_summary", fmt)
            write_delimited(gene_summary_all, out, fmt, args.encoding)

        if exon_rows_all:
            out = choose_output_path(args.output, args.exon_table_out, "exon_table", fmt)
            write_delimited(exon_rows_all, out, fmt, args.encoding)

        if intron_rows_all:
            out = choose_output_path(args.output, args.intron_table_out, "intron_table", fmt)
            write_delimited(intron_rows_all, out, fmt, args.encoding)

    if not args.quiet:
        print(f"Done. Processed {len(paths)} file(s).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
