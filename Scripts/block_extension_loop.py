############################################ Loading Libraries and Input Data ##################################

#!/usr/bin/env python3
from utils.blast_utils import retrieve_nhits_alignments, blastn_subject, makeblast_db, blastn, default_blast_columns
from utils.block_utils import align_seed_neighbours, find_block_boundaries, extend_seeds_to_blocks
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import numpy as np
import argparse
import re

## TO DO
# Make setting to allow dynamic adjustments of parameters according to the block aligment noise. A noiser block  
# should have laxer parameters, meanwhile an unfinished block could bennefit with stringent parameters

test = False
if test:
    input_reads_dir = Path('/home/jpereira/OEs/Results/OE1/toxoBrasil06/Data/mcl_clustering/sequneces_clusters/')
    input_seeds_dir = Path('/home/jpereira/OEs/Results/OE1/toxoBrasil06/Data/make_seeds/')
    params_threads  = 10
    params_min_hits                    = 30
    params_max_aln_extension           = 5000
    params_seed_map_coverage_threshold = 0.8
    params_window_radius               = 15
    params_signal_threshold            = 0.10
    params_aln_w_step_size             = 300
    params_try_mean_signal             = True
    params_mean_signal_threshold       = 0.06
    params_plot                        = True
    work_name = None
    output_data_dir = Path('/home/jpereira/OEs/Results/OE1/toxoBrasil06/Data/seed_extension')
    output_vis_dir = Path('/home/jpereira/OEs/Results/OE1/toxoBrasil06/Visuals/seed_extension')
    
else:
    parser = argparse.ArgumentParser(description="Determine full recombinant-tract boundaries")
    parser.add_argument("--input-reads-dir", type=Path, required=True,help="FASTA or FASTQ file of sequencing reads to search")
    parser.add_argument("--input-seeds-dir", type=Path, required=True)
    parser.add_argument("--params-threads", default=10)
    parser.add_argument("--params-min-hits", required=True, type=float)
    parser.add_argument("--params-max-aln-extension", required=True, type=int)
    parser.add_argument("--params-seed-map-coverage-threshold", type=float, default=0.8, help="Fraction of reads that must cover an MSA column to be valid")
    parser.add_argument("--params-window-radius", type=int, default=15)
    parser.add_argument("--params-signal-threshold", type=float, default=0.1)
    parser.add_argument("--params-aln-w-step-size", type=int, default=300)
    parser.add_argument("--params-try-mean-signal", type=bool, default=True)
    parser.add_argument("--params-mean-signal-threshold", type=float, default=0.06)
    parser.add_argument("--params-plot", type=bool, default=True)
    parser.add_argument("--work-name", type=str, default=None)
    parser.add_argument("--output-data-dir",type=Path, required=True)
    parser.add_argument("--output-vis-dir", type=Path, required=True)

    args = parser.parse_args()
    
    
    input_reads_dir                    = args.input_reads_dir
    input_seeds_dir                    = args.input_seeds_dir
    params_threads                     = args.params_threads
    params_min_hits                    = args.params_min_hits
    params_max_aln_extension           = args.params_max_aln_extension
    params_seed_map_coverage_threshold = args.params_seed_map_coverage_threshold
    params_window_radius               = args.params_window_radius
    params_signal_threshold            = args.params_signal_threshold
    params_aln_w_step_size             = args.params_aln_w_step_size
    params_try_mean_signal             = args.params_try_mean_signal
    params_mean_signal_threshold       = args.params_mean_signal_threshold
    params_plot                        = args.params_plot
    work_name                          = args.work_name
    output_data_dir                    = args.output_data_dir
    output_vis_dir                     = args.output_vis_dir

if work_name:
    output_data_dir = output_data_dir / work_name
    output_vis_dir = output_vis_dir / work_name

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

################### Looping througt the files present in seeds and reads directories ####################

output_data_dir.mkdir(exist_ok=True, parents=True)
output_vis_dir.mkdir(exist_ok=True, parents=True)

fasta_names_set = {f.stem for f in input_reads_dir.glob("*.fasta")}
seeds_names_set = {f.stem for f in input_seeds_dir.glob("*.tsv")}

common_names = list(fasta_names_set & seeds_names_set)
if len(common_names) == 0:
    raise ValueError("There aren't files with common names in the seeds directory or reads directory")

for file_name in common_names:
    input_reads_fasta = input_reads_dir / f"{file_name}.fasta" 
    input_seeds_tsv = input_seeds_dir / f"{file_name}.tsv" 
    
    with open(input_reads_fasta) as fh:
        if not any(fh):
            print(f"Warning! Empty file. Skipping file {file_name}, {input_reads_fasta} is empty")
            continue
    
    with open(input_seeds_tsv) as fh:
        if not any(fh):
            print(f"Warning! Empty file. Skipping file {file_name}, {input_seeds_tsv} is empty")
            continue
    
    output_cl_data_dir = output_data_dir / file_name
    output_cl_vis_dir = output_vis_dir / file_name
    
    output_cl_data_dir.mkdir(exist_ok=True)
    output_cl_vis_dir.mkdir(exist_ok=True)

    
    seeds_df = pd.read_csv(input_seeds_tsv, sep='\t', header=None)
    seeds_df.columns = ['seq', 'count']
    seeds_df = seeds_df.sort_values('count', ascending=False)
    seeds_df = seeds_df.reset_index(drop=True).reset_index()
    seeds_df = seeds_df.rename(columns={'index' : 'qseqid' })

    kmers_fasta = output_data_dir / file_name / "temp.kmers.fasta"
    with open(kmers_fasta, 'w') as fasta:
        for i, row in seeds_df.iterrows():
            fasta.write(f">k:{row['qseqid']}_c:{row['count']}\n")
            fasta.write(f"{row['seq']}\n")

    blast_db_dir = output_data_dir / file_name / 'input_reads_db'
    blast_mapped_seeds_tsv = output_data_dir / file_name / 'blastn.mapped_seeds.tsv'

    # Map seeds to raw reads
    blast_db_out = makeblast_db(input_reads_fasta, db_out=blast_db_dir, remove_old_db=True,show_command=False)
    blastn(blast_input_seqs=kmers_fasta,
        blast_db_file=blast_db_out,
        blast_output_table_tsv=blast_mapped_seeds_tsv,
        word_size=10, show_command=False)

    blast_df = pd.read_csv(blast_mapped_seeds_tsv, sep='\t', header=None)
    blast_df.columns = default_blast_columns

    # Filter alignments with low identity and with low coverage
    blast_df = blast_df[(blast_df['pident'] > 90) & (blast_df['length']/blast_df['qlen'] > 0.95) ]

    hhits_df, lhits_df = retrieve_nhits_alignments(blast_df, n=params_min_hits)

    # Get the unique query IDs and their hit counts
    seed_set = set(hhits_df['qseqid'].astype(str))
    seeds_hits_df = hhits_df.groupby('qseqid')['sseqid'].count().sort_values(ascending=False).to_frame().reset_index()

    # Extract sequences from the FASTA file that match the seed IDs
    seeds_id = []
    seeds_seq = []

    for record in SeqIO.parse(kmers_fasta, "fasta"):
        if record.id in seed_set:
            seeds_id.append(record.id)
            seeds_seq.append(str(record.seq))  # Ensure sequences are strings

    # Create a DataFrame from the extracted records
    seeds_df = pd.DataFrame({
        'qseqid': seeds_id,
        'seq': seeds_seq
    })

    #### - - Extending seeds to blocks - - - - - - - - - - - - - - - - - -
    
    results = extend_seeds_to_blocks(
        output_data_dir=output_cl_data_dir,
        output_vis_dir=output_cl_vis_dir,
        seeds_df=seeds_df,
        seeds_hits_df=seeds_hits_df,
        blast_df=blast_df,
        input_reads_fasta=input_reads_fasta,
        align_seed_neighbours=align_seed_neighbours,
        find_block_boundaries=find_block_boundaries,
        blastn_subject=blastn_subject,
        default_blast_columns=default_blast_columns,
        params_aln_w_step_size=params_aln_w_step_size,
        params_window_radius=params_window_radius,
        params_signal_threshold=params_signal_threshold,
        params_mean_signal_threshold=params_mean_signal_threshold,
        params_max_aln_extension=params_max_aln_extension,
        params_seed_map_coverage_threshold=params_seed_map_coverage_threshold,
        threads=params_threads, max_iter=100, gap_open=1.0, gap_extend=0.1, make_plots=True,
    )
    print("Blocks FASTA:", results["blocks_fasta"])
    print("Mapped seeds:", len(results["mapped_seeds"]))
    print("Blocks found:", len(results["blocks"]))
    