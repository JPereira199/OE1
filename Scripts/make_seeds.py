from utils.bedtools_utils import bbtools_kmercountexact
from pathlib import Path
from Bio import SeqIO
import argparse
import re


test = False
if test:
    print("Warning! Developer mode activated")
    input_fasta_dir = Path('/home/jpereira/OEs/Results/OE1/NamSeqs/Data/mcl_clustering/sequneces_clusters/')
    params_min_fasta_size = 35
    params_kmer_size = 35
    params_min_kmer_counts = 35
    params_threads = 5
    output_seeds_dir = Path('/home/jpereira/OEs/Results/OE1/NamSeqs/Data/make_seeds/')
else:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-fasta-dir', required=True, help="directory with the sequences useed to make seeds")    
    parser.add_argument('--params-min-fasta-size', default=30, type=int, help='Minimum number of reads reqeried in each fasta')
    parser.add_argument('--params-kmer-size', default=30, type=int, help='Size of K-mer')
    parser.add_argument('--params-min-kmer-counts', default=30, type=int, help='Min number of K-mer apparences')
    parser.add_argument('--params-threads', default=5, type=int )
    parser.add_argument('--output-seeds-dir', required=True, help='Output directory with the seeds tables in tsv format')
    args = parser.parse_args()

    input_fasta_dir        = Path(args.input_fasta_dir)
    params_min_fasta_size  = int(args.params_min_fasta_size)
    params_kmer_size       = int(args.params_kmer_size)
    params_min_kmer_counts = int(args.params_min_kmer_counts)
    params_threads         = int(args.params_threads)
    output_seeds_dir       = Path(args.output_seeds_dir)

output_seeds_dir.mkdir(exist_ok=True)


# Loop over files in the input directory
for fasta in input_fasta_dir.iterdir():
    
    if not fasta.is_file():
        continue

    if not re.search(r"\.(fasta|fa|fastq)$", fasta.name):
        print(f"❌ Skipping non-FASTA/Q file: {fasta}")
        continue
    
    # Count the number of entries in the sequence file
    q = sum(1 for _ in SeqIO.parse(fasta, "fasta"))  # or "fastq" depending on format
    
    if q < params_min_fasta_size:
        print(f"⚠️ Skipping {fasta.name}: only {q} entries (min is {params_min_fasta_size})")
        continue

    # Remove file extension using regex
    file_name = re.sub(r"\.(fasta|fa|fastq)$", "", fasta.name)
    
    output_seed_tsv = output_seeds_dir / f"{file_name}.tsv"
    
    bbtools_kmercountexact(
        infile=fasta,
        outfile=output_seed_tsv,
        k=params_kmer_size,
        min_counts=params_min_kmer_counts,
        khist=None,
        make_fasta=False,
        threads=params_threads,
        extra=None,
        exe="kmercountexact.sh")
    
    # Ensurig that kmercountexact.sh saved the K-mers counts in a tsv table
    
    