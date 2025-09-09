################################################ Helper functions used in the script ################################################
import pandas as pd
def total_coverage_generic(cols: pd.DataFrame, start_coord_col_name: str, end_coord_col_name: str, seq_len_col_name: str) -> float:
    """
    Calculates the total coverage fraction of a sequence by merging overlapping intervals.

    Parameters:
        cols (pd.DataFrame): A DataFrame sorted by 'seqid', then by start and end coordinates.
        start_coord_col_name (str): Name of the column with the start coordinate.
        end_coord_col_name (str): Name of the column with the end coordinate.
        seq_len_col_name (str): Name of the column with the sequence length.

    Assumptions:
        - The DataFrame contains only one 'seqid' or is grouped by 'seqid'.
        - Coordinates are sorted beforehand.

    Returns:
        float: The fraction of the sequence covered by the merged intervals.
    """
    
    # Handle empty DataFrame
    if cols.empty:
        return 0.0

    sstart = cols[start_coord_col_name].values
    send = cols[end_coord_col_name].values
    total_len = cols[seq_len_col_name].iloc[0]
    
    total_cover = 0
    current_start, current_end = sstart[0], send[0]

    for i in range(1, len(sstart)):
        if sstart[i] <= current_end:
            current_end = max(current_end, send[i])
        else:
            total_cover += current_end - current_start
            current_start, current_end = sstart[i], send[i]

    total_cover += current_end - current_start
    return total_cover / total_len if total_len > 0 else 0.0

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

# Script objectives:
# - Use the parameters to filter STR sequences from an input FASTA using the STR table.
#   * If multiple tandem repeats are found by TRF software for the same sequence,
#     each tandem is analyzed individually and the sequence is discarded if ANY of the
#     tandems exceed the user-defined parameters.
#
# - Create a FASTA file with the filtered (passed) sequences.
# - Create a FASTA file with the discarded sequences.


################################################ Loading libraries and arguments ################################################

import pandas as pd
import argparse
from Bio import SeqIO
import numpy as np

test = False
if test:
    print('Warning! Developing mode activated! ')
    input_fasta = '/home/jpereira/OEs/OE1.v2/Data/check_gc/decont.low_gc.fasta'
    input_str_tsv = '/home/jpereira/OEs/OE1.v2/Data/short_tandem_repeats/short_tandem_repeats.tsv'
    params_period_size = 0
    params_total_tandem_coverage = 0
    params_tandem_coverage = 0.85
    output_passed_sequences = '/home/jpereira/OEs/OE1.v2/Data/short_tandem_repeats/str_filtered.fasta'
    output_discarded_sequences = '/home/jpereira/OEs/OE1.v2/Data/short_tandem_repeats/str_filtered.fasta'
    output_repeats_tsv = '/home/jpereira/OEs/OE1.v2/Data/short_tandem_repeats/filter_short_tandem_repeats.tsv'
    output_repeats_fasta = '/home/jpereira/OEs/OE1.v2/Data/short_tandem_repeats/repeats.fasta'

else:
    parser = argparse.ArgumentParser(
        description="""
        Discard sequences from a FASTA file if they have short tandem repeats (STR).
        """)
    parser.add_argument("--input-fasta", help="Path to the input FASTA file.")
    parser.add_argument("--input-str-table", help="Path to a TSV table with short tandem repeats present in the FASTA file.")
    
    parser.add_argument("--param-period-size", type=float, default=0, help="Minimum period size of an STR required to discard the sequence.")
    parser.add_argument("--param-tandem-coverage", type=float, default=0.85, help="Coverage ratio of an STR found in the sequence (0-1) required to discard it.")
    parser.add_argument("--param-total-tandem-coverage", type=float, default=0.85, help="Coverage ratio of all STR found in the sequence (0-1) required to discard the sequence.")
    
    parser.add_argument("--output-passed-sequences", help="Output FASTA file of sequences WITHOUT STR", default="passed_sequences.fasta")
    parser.add_argument("--output-discarded-sequences", help="Output FASTA file of sequences WITH STR", default="discarded_sequences.fasta")
    parser.add_argument("--output-short-repeats", help="Output FASTA file with the short repeats in each discarded sequence")
    parser.add_argument("--output-str-table", help="A table with the STR data of sequences that didn't pass the filter")
    args = parser.parse_args()

    input_fasta = args.input_fasta
    input_str_tsv = args.input_str_table
        
    params_tandem_coverage = args.param_tandem_coverage,
    params_total_tandem_coverage = args.param_total_tandem_coverage
    params_period_size = args.param_period_size
    
    output_passed_sequences = args.output_passed_sequences
    output_discarded_sequences = args.output_discarded_sequences
    output_repeats_fasta = args.output_short_repeats
    output_repeats_tsv = args.output_str_table
    
    
###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

################################################ Loading inputs and filtering tandems ################################################

# Loading input FASTA sequences 
infasta_dict = {'seqid' : [], 'qlen' : []}
with open(input_fasta, 'r') as fasta:
    for record in SeqIO.parse(fasta,'fasta'):
        infasta_dict['seqid'].append(record.id)
        infasta_dict['qlen'].append(len(record.seq))
infasta_df = pd.DataFrame(infasta_dict)

# Loading input STR table
str_df = pd.read_csv(input_str_tsv, sep='\t')

# Add len sizes to str_df 
str_df = pd.merge(str_df, infasta_df, how='left', on='seqid')

# Evaluate each tandem line indivualy
# The dataframe is grouped by the column row in order to calculate the coverage of each tandem indiviauly.
# If were grouped by 'seqid' all tandem coverge would be summed
str_df['row'] = np.linspace(1, str_df.shape[0], num=str_df.shape[0])
tandem_coverage = str_df.sort_values(by=['seqid','rstart', 'rend']).groupby('row').apply( lambda x: total_coverage_generic(x, 'rstart', 'rend', 'qlen')).reset_index(name='tandem_coverage_fraction')
str_df = pd.merge(str_df, tandem_coverage, on='row', how='left')

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -

###################################### Saving filtered sequences, repeats and filtered str-table ######################################


#Extract relevant STR and save them as a TSV 
filtered_str_df = str_df[str_df['tandem_coverage_fraction'] > params_tandem_coverage]
filtered_str_df.to_csv(output_repeats_tsv, sep='\t', index=False)

#Write relevant STR sequneces as individual FASTA entries
repeats_df = filtered_str_df[['seqid', 'repeat_sequence']]

ids_counts = {}
with open(output_repeats_fasta, 'w') as fasta:
    for _, row in repeats_df.iterrows():
        seqid = row['seqid']
        sequence = row['repeat_sequence']

        # Ensure unique identifiers
        ids_counts[seqid] = ids_counts.get(seqid, 0) + 1
        fasta_id = f"{seqid}_{ids_counts[seqid]}"

        fasta.write(f">{fasta_id}\n{sequence}\n")

# Create a set for fast lookup
tandem_ids = set(repeats_df['seqid'])

# Separate sequences from the input FASTA based on whether they have tandem repeats
with open(input_fasta, 'r') as in_fasta, \
     open(output_discarded_sequences, 'w') as discard_fasta, \
     open(output_passed_sequences, 'w') as filter_fasta:

    for record in SeqIO.parse(in_fasta, 'fasta'):
        if record.id in tandem_ids:
            SeqIO.write(record, discard_fasta, 'fasta')
        else:
            SeqIO.write(record, filter_fasta, 'fasta')

###### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -