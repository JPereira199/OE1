#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import shutil

import gzip

def open_text_maybe_gzip(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")

def open_bin_maybe_gzip(path: Path):
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


def detect_seq_format(path: Path) -> str:
    """
    Detect FASTA vs FASTQ from the first non-empty line.
    Works for .fa/.fasta/.fastq and .gz versions.
    """
    with open_text_maybe_gzip(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.startswith("@"):
                return "fastq"
            break
    raise ValueError(f"Could not detect format (FASTA/FASTQ) from file header: {path}")



def run_nanofilt_fastq(in_fastq: Path, out_fastq: Path, min_mean_q: float, threads: int = 1) -> None:
    """
    Supports .fastq and .fastq.gz.
    Note: NanoFilt is stream-based; threads is kept for interface consistency but not used.
    """
    nanofilt_cmd = ["NanoFilt", "-q", str(min_mean_q)]

    out_fastq.parent.mkdir(parents=True, exist_ok=True)

    # Pick a decompressor that exists
    pigz = shutil.which("pigz")
    gzip_bin = shutil.which("gzip")
    decompressor = None
    if in_fastq.suffix == ".gz":
        if pigz:
            decompressor = [pigz, "-dc", str(in_fastq)]
        elif gzip_bin:
            decompressor = [gzip_bin, "-dc", str(in_fastq)]
        else:
            raise RuntimeError("Neither 'pigz' nor 'gzip' found in PATH; cannot decompress .gz input.")

    with out_fastq.open("wb") as fout:
        if decompressor:
            # decompressor | NanoFilt > out_fastq
            p1 = subprocess.Popen(decompressor, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p2 = subprocess.Popen(nanofilt_cmd, stdin=p1.stdout, stdout=fout, stderr=subprocess.PIPE)

            assert p1.stdout is not None
            p1.stdout.close()

            err2 = p2.communicate()[1]
            err1 = p1.communicate()[1]

            if p1.returncode != 0:
                raise RuntimeError(f"Decompression failed:\n{err1.decode(errors='replace')}")
            if p2.returncode != 0:
                raise RuntimeError(f"NanoFilt failed:\n{err2.decode(errors='replace')}")
        else:
            # Plain FASTQ
            with in_fastq.open("rb") as fin:
                p = subprocess.run(nanofilt_cmd, stdin=fin, stdout=fout, stderr=subprocess.PIPE)
            if p.returncode != 0:
                raise RuntimeError(f"NanoFilt failed:\n{p.stderr.decode(errors='replace')}")


def length_filter_fasta(in_fasta: Path, out_fasta: Path, min_len: int) -> tuple[list[int], list[int]]:
    lens_before: list[int] = []
    lens_after: list[int] = []

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    with open_text_maybe_gzip(in_fasta) as fin, out_fasta.open("w") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            L = len(rec.seq)
            lens_before.append(L)
            if L >= min_len:
                lens_after.append(L)
                SeqIO.write(rec, fout, "fasta")

    return lens_before, lens_after





def run_seqkit_fq2fa(in_fastq: Path, out_fasta: Path) -> None:
    cmd = ["seqkit", "fq2fa", str(in_fastq)]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False)
    if p.returncode != 0:
        raise RuntimeError(f"seqkit fq2fa failed:\n{p.stderr.decode(errors='replace')}")
    out_fasta.write_bytes(p.stdout)



def add_annotations(data: list[int], ax) -> None:
    if not data:
        ax.text(
            0.5, 0.5, "No sequences",
            transform=ax.transAxes, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"),
        )
        return

    max_val = max(data)
    min_val = min(data)
    mean_val = sum(data) / len(data)
    median_val = float(np.median(data))

    textstr = f"n: {len(data)}\nMax: {max_val}\nMin: {min_val}\nMean: {mean_val:.0f}\nMedian: {median_val:.0f}"
    ax.text(
        0.88, 0.90, textstr,
        transform=ax.transAxes, fontsize=11,
        verticalalignment="top", horizontalalignment="left",
        bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"),
    )


def plot_histograms(lens_before: list[int], lens_after: list[int], out_fig: Path, title_prefix: str) -> None:
    out_fig.parent.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    sns.histplot(lens_before, ax=axes[0], bins=50)
    axes[0].set_title(f"{title_prefix} — Before length filter")

    sns.histplot(lens_after, ax=axes[1], bins=50)
    axes[1].set_title(f"{title_prefix} — After length filter")

    add_annotations(lens_before, axes[0])
    add_annotations(lens_after, axes[1])

    plt.tight_layout()
    plt.savefig(out_fig)
    plt.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter long reads by quality (FASTQ only), convert to FASTA, filter by length, and plot histograms."
    )
    p.add_argument("--input", required=True, type=Path, help="Input FASTA or FASTQ.")
    p.add_argument("--min_len", required=True, type=int, help="Minimum read length (bp).")
    p.add_argument("--out_fig", required=True, type=Path, help="Output histogram figure (e.g., .pdf/.png).")
    p.add_argument("--out_fasta", required=True, type=Path, help="Output filtered FASTA.")

    # Only used if input is FASTQ
    p.add_argument("--min_mean_q", type=int, default=10, help="Minimum mean read quality (FASTQ only).")
    p.add_argument("--threads", type=int, default=1, help="Kept for interface consistency (NanoFilt is stream-based).")

    # Keep intermediates (useful for debugging)
    p.add_argument("--keep_intermediate", action="store_true", help="Keep intermediate FASTQ/FASTA files.")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    in_path: Path = args.input
    out_fig: Path = args.out_fig
    out_fasta: Path = args.out_fasta
    min_len: int = args.min_len
    min_mean_q: float = args.min_mean_q

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_fig.parent.mkdir(parents=True, exist_ok=True)

    fmt = detect_seq_format(in_path)

    # Prepare a FASTA to length-filter
    if fmt == "fastq":
        if min_mean_q <= 0:
            raise ValueError("Input is FASTQ, but --min_mean_q was not set (> 0).")

        workdir = out_fasta.parent
        tmp_fastq = workdir / f"{in_path.stem}.q{min_mean_q}.fastq"
        tmp_fasta = workdir / f"{in_path.stem}.q{min_mean_q}.fasta"

        print(f"[FASTQ] Filtering by mean Q >= {min_mean_q} with NanoFilt...")
        run_nanofilt_fastq(in_path, tmp_fastq, min_mean_q=min_mean_q, threads=args.threads)

        print("[FASTQ] Converting FASTQ -> FASTA with seqkit fq2fa...")
        run_seqkit_fq2fa(tmp_fastq, tmp_fasta)

        fasta_for_len_filter = tmp_fasta
        title_prefix = f"{in_path.name} (Q>={min_mean_q})"

        if not args.keep_intermediate:
            # We'll delete tmp_fastq after conversion succeeded
            try:
                tmp_fastq.unlink(missing_ok=True)
            except Exception:
                pass

    else:
        print("[FASTA] Skipping quality filtering and conversion.")
        fasta_for_len_filter = in_path
        title_prefix = in_path.name

    # Length filter + histograms
    print(f"[LEN] Filtering reads with length >= {min_len} bp...")
    lens_before, lens_after = length_filter_fasta(fasta_for_len_filter, out_fasta, min_len=min_len)

    print(f"[PLOT] Writing histogram figure to: {out_fig}")
    plot_histograms(lens_before, lens_after, out_fig, title_prefix=title_prefix)

    # Cleanup tmp_fasta unless requested
    if fmt == "fastq" and (not args.keep_intermediate):
        try:
            Path(fasta_for_len_filter).unlink(missing_ok=True)
        except Exception:
            pass

    print("[DONE]")
    print(f"  Output FASTA: {out_fasta}")
    print(f"  Output FIG : {out_fig}")
    print(f"  Reads before length filter: {len(lens_before)}")
    print(f"  Reads after  length filter: {len(lens_after)}")


if __name__ == "__main__":
    main()
