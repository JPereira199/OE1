#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import importlib.util
import pandas as pd
from utils.blast_utils import blastn_dotplot, find_unmapped

def str2bool(v: str) -> bool:
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ("yes", "y", "true", "t", "1"):
        return True
    if v in ("no", "n", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("boolean value expected")


def write_id_list(ids: set[str], out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for x in sorted(ids):
            fh.write(f"{x}\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run blastn_dotplot and report unmapped query/subject IDs.")
    p.add_argument("--fasta1", required=True, type=Path, help="Query FASTA for dotplot (passed to blastn_dotplot).")
    p.add_argument("--fasta2", required=True, type=Path, help="Subject FASTA for dotplot (passed to blastn_dotplot).")
    p.add_argument("--output_base_dir", required=True, type=Path, help="Base output directory.")
    p.add_argument("--word_size", type=int, default=20)
    p.add_argument("--qnames", choices=("False", "FALSE", "false", "True", "TRUE", "true"), default="True")
    p.add_argument("--snames", choices=("False", "FALSE", "false", "True", "TRUE", "true"), default="True")
    p.add_argument("--annot_size", type=float, default=12)
    p.add_argument("--title", type=str, default="Validation Seqs vs AutoBlocks")

    return p.parse_args()



def main():
    # Make matplotlib safe for Snakemake/headless runs (blast_utils uses plt.show()).
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.show = lambda *a, **k: None  # prevent GUI show

    args = parse_args()
    #bu = load_blast_utils(args.blast_utils_path)

    out_workdir = args.output_base_dir
    out_workdir.mkdir(parents=True, exist_ok=True)

    # 1) Dotplot (this also writes blast TSV + dotplot.svg inside out_workdir)
    dot_df = blastn_dotplot(
        title=str(args.title),
        fasta1=str(args.fasta1),
        fasta2=str(args.fasta2),
        output_base_dir=str(args.output_base_dir),
        word_size=int(args.word_size),
        qnames=str2bool(args.qnames),
        snames=str2bool(args.snames),
        annot_size=float(args.annot_size),
    )

    # Save the dotplot table (handy as a tracked output)
    dotplot_table_tsv = out_workdir / "dotplot_table.tsv"
    dot_df.to_csv(dotplot_table_tsv, sep="\t", index=False)

    # 2) Unmapped IDs
    unm_q, unm_s = find_unmapped(
        df=dot_df,
        fasta_queries=str(args.fasta1),
        fasta_subjects=str(args.fasta2),
    )

    unmapped_queries_txt = out_workdir / "unmapped_queries.txt"
    unmapped_subjects_txt = out_workdir / "unmapped_subjects.txt"
    write_id_list(unm_q, unmapped_queries_txt)
    write_id_list(unm_s, unmapped_subjects_txt)

    # Small summary TSV (nice for downstream reports)
    summary_tsv = out_workdir / "unmapped_summary.tsv"
    pd.DataFrame(
        [{
            "work_name": args.title,
            "num_alignments": int(len(dot_df)),
            "unmapped_queries": int(len(unm_q)),
            "unmapped_subjects": int(len(unm_s)),
        }]
    ).to_csv(summary_tsv, sep="\t", index=False)

    # Ensure we don't keep figures open
    plt.close("all")


if __name__ == "__main__":
    main()
