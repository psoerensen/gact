"""Extract common autosomal SNVs from IGSR's 2019-03-12 GRCh38 release."""

import argparse
from collections import Counter
import gzip
import hashlib
from itertools import groupby
import json
from pathlib import Path
import re


IGSR = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
)
UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
RSID = re.compile(rb"rs[0-9]+")
BASES = {b"A", b"C", b"G", b"T"}


def md5(path):
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source(path, url, expected=None):
    actual = md5(path)
    if expected is not None and actual != expected:
        raise ValueError(f"MD5 mismatch: {path}")
    return {"file": path.name, "url": url, "md5": actual}


def reference_sequence(path, chrom):
    with gzip.open(path, "rb") as stream:
        if stream.readline().strip() != f">chr{chrom}".encode():
            raise ValueError(f"Unexpected FASTA chromosome: {path}")
        return stream.read().translate(None, b"\r\n").upper()


def rsid_rows(path, chrom):
    with gzip.open(path, "rb") as stream:
        if stream.readline().rstrip() != b"#CHROM\tPOS\tID\tREF\tALT":
            raise ValueError(f"Unexpected rsID lookup header: {path}")
        for line in stream:
            label, position, ids, ref, alt = line.rstrip().split(b"\t")
            if label != f"chr{chrom}".encode():
                raise ValueError(f"Unexpected rsID lookup chromosome: {path}")
            yield int(position), (ref, alt), ids


def prepare(args):
    chrom = args.chromosome
    panel_name = f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    panel = args.panel / panel_name
    fasta = args.fasta / f"chr{chrom}.fa.gz"
    lookup = args.rsids / f"gnomad_rsid_chr{chrom}.tsv.bgz"
    vcf_sums = dict((row[0].removeprefix("./"), row[2]) for row in (
        line.split() for line in
        (args.panel / "20190312_biallelic_SNV_and_INDEL_MANIFEST.txt").read_text().splitlines()
    ) if len(row) == 3)
    fasta_sums = dict((row[1], row[0]) for row in (
        line.split() for line in (args.fasta / "md5sum.txt").read_text().splitlines()
    ) if len(row) == 2)
    sources = [
        source(panel, f"{IGSR}/{panel.name}", vcf_sums[panel.name]),
        source(fasta, f"{UCSC}/{fasta.name}", fasta_sums[fasta.name]),
        source(lookup, "https://gnomad.broadinstitute.org/downloads#v4"),
    ]
    sequence = reference_sequence(fasta, chrom)
    annotations = groupby(rsid_rows(lookup, chrom), key=lambda row: row[0])
    position, rows = next(annotations, (None, ()))
    matched = {}
    counts = Counter()
    args.out.mkdir(parents=True, exist_ok=True)
    out = args.out / f"chr{chrom}.tsv.gz"
    receipt = args.out / f"chr{chrom}.json"
    if out.exists() or receipt.exists():
        raise FileExistsError(f"Output already exists for chromosome {chrom}")
    staged = args.out / f"chr{chrom}.tsv.gz.part"
    try:
        with gzip.open(panel, "rb") as variants, gzip.open(staged, "wb", compresslevel=1) as output:
            output.write(b"rsids\tchr\tpos\tref\talt\taf\n")
            is_grch38 = False
            for line in variants:
                if line.startswith(b"##reference="):
                    is_grch38 = b"GRCh38" in line
                if line.startswith(b"#"):
                    continue
                if not is_grch38:
                    raise ValueError(f"Expected a GRCh38 VCF: {panel}")
                label, pos, rsid, ref, alt, _, status, info, _ = line.split(b"\t", 8)
                coordinate = int(pos)
                if label != str(chrom).encode() or not 0 < coordinate <= len(sequence):
                    raise ValueError(f"Unexpected chromosome or position: {label!r}:{pos!r}")
                counts["input"] += 1
                if status != b"PASS" or ref not in BASES or alt not in BASES or ref == alt:
                    counts["not_pass_biallelic_snv"] += 1
                    continue
                fields = dict(field.split(b"=", 1) for field in info.split(b";") if b"=" in field)
                ac, an = int(fields[b"AC"]), int(fields[b"AN"])
                if not 0 <= ac <= an or an == 0:
                    raise ValueError(f"Invalid allele counts: {label!r}:{pos!r}")
                af = ac / an
                if min(af, 1 - af) < args.min_maf:
                    counts["below_min_maf"] += 1
                    continue
                if sequence[coordinate - 1:coordinate] != ref:
                    counts["reference_mismatch"] += 1
                    continue
                while position is not None and position < coordinate:
                    position, rows = next(annotations, (None, ()))
                    matched = {}
                if position == coordinate and not matched:
                    for _, alleles, ids in rows:
                        matched.setdefault(alleles, set()).update(
                            token for token in ids.split(b";") if RSID.fullmatch(token)
                        )
                candidates = matched.get((ref, alt), set()) if position == coordinate else set()
                if RSID.fullmatch(rsid):
                    candidates = candidates | {rsid}
                if len(candidates) != 1:
                    counts["unresolved_rsid"] += 1
                    continue
                output.write(b"\t".join((next(iter(candidates)), label, pos, ref, alt,
                                         repr(af).encode())) + b"\n")
                counts["retained"] += 1
        if not counts["retained"]:
            raise ValueError(f"No eligible variants on chromosome {chrom}")
        receipt_data = {
            "chromosome": chrom, "build": "GRCh38", "min_maf": args.min_maf, "counts": counts,
            "sources": sources, "output": {"file": out.name, "md5": md5(staged)},
        }
        staged.replace(out)
        receipt.write_text(json.dumps(receipt_data, indent=2) + "\n", encoding="utf-8")
    finally:
        staged.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chromosome", type=int, choices=range(1, 23))
    for directory in ("panel", "fasta", "rsids", "out"):
        parser.add_argument(f"--{directory}", type=Path, required=True)
    parser.add_argument("--min-maf", type=float, default=0.01)
    args = parser.parse_args()
    if not 0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must be between 0 and 0.5")
    prepare(args)
