#!/usr/bin/env python3

import argparse
from genericpath import isfile
from Bio import SeqIO
import os


test = False
if test:
    print('Developer mode activated!')
    inputs = {
        'fasta_in' : '/home/jpereira/OEs/OE1.v2/Data/filter_sizes/toxo1.size_filtered.fasta',
        'paf_in'   : '/home/jpereira/OEs/OE1.v2/Data/map_contamination/query_x_ref-contamination.paf'}
    outputs = {
        'fasta_out' : '/home/jpereira/OEs/OE1.v2/Data/decontamination/test-decont.fasta',
        'stats_out' : '/home/jpereira/OEs/OE1.v2/Data/decontamination/test-decont-stats.tsv'}
    params = {
        'min_coverage' : 70}
else:
    parser = argparse.ArgumentParser(
        description="""
        Discard sequences from a FASTA file if they align well (≥ min_coverage) 
        to any reference in the given PAF alignment file.
        """)
    parser.add_argument("--fasta_in", help="Path to the input (filtered) FASTA file.")
    parser.add_argument("--paf_in", help="Path to the PAF alignment file.")
    parser.add_argument("--fasta_out", help="Path to the output (decontaminated) FASTA file.")
    parser.add_argument("--stats_out", help="Path to the output TSV stats file.")
    parser.add_argument("--min_coverage",type=float,default=70.0,
        help="Minimum coverage percentage (0-100) to consider a sequence contaminated. Default: 70.")
    args = parser.parse_args()

    inputs = {
        'fasta_in' : args.fasta_in,
        'paf_in'   : args.paf_in}

    outputs = {
        'fasta_out' : args.fasta_out,
        'stats_out' : args.stats_out}

    params = {
        'min_coverage' : args.min_coverage}



def parse_paf_line(line):
    """
    Parse one line from a PAF file and return:
      - qname (query/contig ID)
      - qlen (query length)
      - qstart (start of the alignment on the query)
      - qend (end of the alignment on the query)
    The PAF format has many fields, but these four are enough to estimate coverage.
    """
    fields = line.strip().split('\t')
    qname = fields[0]
    qlen = int(fields[1])
    qstart = int(fields[2])
    qend = int(fields[3])
    matching_bases = int(fields[9])
    return qname, qlen, qstart, qend, matching_bases



def main(inputs: dict, outputs: dict, params: dict):

    # Make a directory for the output files
    os.makedirs(os.path.dirname(outputs['fasta_out']), exist_ok=True)
    
    # Store the highest coverage found for each query
    coverage_dict = {}

    # 1. Parse the PAF file to identify contaminated sequences
    with open(inputs['paf_in'], 'r') as paf_file:
        for line in paf_file:
            # Skip empty or comment lines
            if not line.strip() or line.startswith('#'):
                continue
            qname, qlen, qstart, qend, matching_bases = parse_paf_line(line)
            
            # Calculate coverage: (aligned length / total query length) * 100
            #aligned_len = abs(qend - qstart)
            coverage = (matching_bases / qlen) * 100 if qlen > 0 else 0
            
            # Keep track of the maximum coverage observed for this sequence
            if qname not in coverage_dict or coverage > coverage_dict[qname]:
                coverage_dict[qname] = coverage

    # 2. Read the FASTA and separate contaminated vs. clean sequences
    contaminated = set()
    total_sequences = 0

    with open(inputs['fasta_in'], 'r') as fasta_file:
        # We won't write out here yet, just counting
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_sequences += 1
            seq_id = record.id
            max_cov = coverage_dict.get(seq_id, 0.0)
            
            # If coverage >= threshold, mark as contaminated
            if max_cov >= params['min_coverage']:
                contaminated.add(seq_id)

    # 3. Write out the decontaminated sequences to a new FASTA
    retained_count = 0    
    with open(inputs['fasta_in'], 'r') as fasta_file, open(outputs['fasta_out'], 'w') as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id not in contaminated:
                retained_count += 1
                SeqIO.write(record, out_fasta, "fasta")

    # 4. Write stats to TSV
    contaminated_count = len(contaminated)
    with open(outputs['stats_out'], 'w') as stats_file:
        stats_file.write("Total\tContaminated\tRetained\tMinCoverage\n")
        stats_file.write(f"{total_sequences}\t{contaminated_count}\t{retained_count}\t{params['min_coverage']}\n")

if __name__ == "__main__":
    main(inputs, outputs, params)
