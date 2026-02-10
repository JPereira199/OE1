#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
import sys

# -----------------------------
# Developer mode
# -----------------------------
test = False  # set True if you want hardcoded paths

if test:
    print("Developer mode activated!")
    inputs = {
        "fasta_in": "/home/jpereira/OEs/Results/OE1/test/ToxoPasteur1_mini/Data/filter_sizes/reads.size_filtered.fasta",
        "paf_in": "/home/jpereira/OEs/Results/OE1/test/ToxoPasteur1_mini/Data/map_contamination/query_x_ref-contamination.paf",
    }
    outputs = {
        "fasta_out": "/home/jpereira/OEs/Results/OE1/test/ToxoPasteur1_mini/Data/decontamination/decont.fasta",
        "stats_out": "/home/jpereira/OEs/Results/OE1/test/ToxoPasteur1_mini/Data/decontamination/decont-stats.tsv",
        "align_stats_out": "/home/jpereira/OEs/Results/OE1/test/ToxoPasteur1_mini/Data/decontamination/alignments-stats.tsv",
    }
    params = {
        "min_alignment_length": 700,
        "max_de": 0.30,
        "min_query_coverage": 0.70,  # fraction (0-1)
    }
else:
    parser = argparse.ArgumentParser(
        description="""
        Discard sequences from a FASTA file if they align well enough to any
        reference in a PAF file, based on:
          - min alignment length (bp) on the query
          - max de
          - min query coverage (aligned/query_len)
        Also writes a per-alignment table with de, alignment length and query coverage.
        """
    )
    parser.add_argument("--fasta_in", required=True, help="Input (filtered) FASTA file.")
    parser.add_argument("--paf_in", required=True, help="Input PAF alignment file.")
    parser.add_argument("--fasta_out", required=True, help="Output (decontaminated) FASTA file.")
    parser.add_argument("--stats_out", required=True, help="Output TSV summary stats file.")
    parser.add_argument(
        "--align_stats_out",
        required=True,
        help="Output TSV with one row per alignment (de, aln_len, qcov, etc.).",
    )

    parser.add_argument(
        "--min_alignment_length",
        type=int,
        default=700,
        help="Minimum alignment length on the query (bp). Default: 700.",
    )
    parser.add_argument(
        "--max_de",
        type=float,
        default=0.30,
        help="Maximum de value allowed for an alignment to be considered. Default: 0.30.",
    )
    parser.add_argument(
        "--min_query_coverage",
        type=float,
        default=0.70,
        help="Minimum query coverage as a fraction (0-1). Default: 0.70. "
             "If >1, interpreted as percent (e.g., 70 -> 0.70).",
    )

    args = parser.parse_args()

    min_query_cov = args.min_query_coverage
    if min_query_cov > 1.0:
        min_query_cov = min_query_cov / 100.0

    inputs = {"fasta_in": args.fasta_in, "paf_in": args.paf_in}
    outputs = {
        "fasta_out": args.fasta_out,
        "stats_out": args.stats_out,
        "align_stats_out": args.align_stats_out,
    }
    params = {
        "min_alignment_length": args.min_alignment_length,
        "max_de": args.max_de,
        "min_query_coverage": min_query_cov,
    }


def parse_paf_line(line: str):
    """
    Parse one PAF line and return key fields:
      qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, de
    """
    fields = line.strip().split()  # supports tabs/spaces
    if len(fields) < 12:
        raise ValueError(f"PAF line has fewer than 12 fields: {line!r}")

    qname = fields[0]
    qlen = int(fields[1])
    qstart = int(fields[2])
    qend = int(fields[3])
    strand = fields[4]
    tname = fields[5]
    tlen = int(fields[6])
    tstart = int(fields[7])
    tend = int(fields[8])
    nmatch = int(fields[9])
    alen = int(fields[10])   # alignment block length on query (includes gaps)
    mapq = int(fields[11])

    de_val = None
    print("START")
    for f in fields[12:]:
        print(f)
        if f.startswith("de:f:"):
            try:
                de_val = float(f.split(":", 2)[2])
            except Exception:
                de_val = None
            break
    print(f"de_val: {de_val}")
    print("END")

    return (qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, de_val)


def main(inputs: dict, outputs: dict, params: dict):
    os.makedirs(os.path.dirname(outputs["fasta_out"]), exist_ok=True)
    os.makedirs(os.path.dirname(outputs["stats_out"]), exist_ok=True)
    os.makedirs(os.path.dirname(outputs["align_stats_out"]), exist_ok=True)

    min_aln_len = params["min_alignment_length"]
    max_de = params["max_de"]
    min_qcov = params["min_query_coverage"]

    # best qcov among alignments that pass aln_len + de filters
    best_qcov = {}

    # 1) Parse PAF and write per-alignment table
    with open(inputs["paf_in"], "r") as paf_file, open(outputs["align_stats_out"], "w") as out_tab:
        out_tab.write(
            "qname\tqlen\tqstart\tqend\tstrand\taln_len\tqcov\t"
            "tname\ttlen\ttstart\ttend\tmapq\tde\t"
            "pass_min_aln_len\tpass_max_de\tpass_min_qcov\tpasses_all\n"
        )

        for line in paf_file:
            if not line.strip() or line.startswith("#"):
                continue

            (qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend,
             nmatch, alen, mapq, de_val) = parse_paf_line(line)

            if qlen <= 0:
                continue

            aln_len = abs(qend - qstart)          # on-query span
            qcov = aln_len / qlen                 # fraction

            pass_len = aln_len >= min_aln_len
            pass_de = (de_val is not None) and (de_val <= max_de)
            pass_cov = qcov >= min_qcov
            passes_all = pass_len and pass_de and pass_cov

            # write row
            out_tab.write(
                f"{qname}\t{qlen}\t{qstart}\t{qend}\t{strand}\t{aln_len}\t{qcov:.6f}\t"
                f"{tname}\t{tlen}\t{tstart}\t{tend}\t{mapq}\t"
                f"{'' if de_val is None else f'{de_val:.6f}'}\t"
                f"{int(pass_len)}\t{int(pass_de)}\t{int(pass_cov)}\t{int(passes_all)}\n"
            )

            # update best_qcov only for alignments that pass len+de (and then we apply cov threshold later)
            if pass_len and pass_de:
                if (qname not in best_qcov) or (qcov > best_qcov[qname]):
                    best_qcov[qname] = qcov

    # 2) Determine contaminated sequences: any passing alignment with qcov >= min_qcov
    contaminated = set()
    total_sequences = 0

    with open(inputs["fasta_in"], "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_sequences += 1
            seq_id = record.id
            if best_qcov.get(seq_id, 0.0) >= min_qcov:
                contaminated.add(seq_id)

    # 3) Write decontaminated FASTA
    retained_count = 0
    with open(inputs["fasta_in"], "r") as fasta_file, open(outputs["fasta_out"], "w") as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id not in contaminated:
                retained_count += 1
                SeqIO.write(record, out_fasta, "fasta")

    # 4) Summary stats
    contaminated_count = len(contaminated)
    with open(outputs["stats_out"], "w") as stats_file:
        stats_file.write(
            "Total\tContaminated\tRetained\tMinAlignmentLength\tMaxDe\tMinQueryCoverage\n"
        )
        stats_file.write(
            f"{total_sequences}\t{contaminated_count}\t{retained_count}\t"
            f"{min_aln_len}\t{max_de}\t{min_qcov}\n"
        )


if __name__ == "__main__":
    main(inputs, outputs, params)
