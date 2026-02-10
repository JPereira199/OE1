import numpy as np
import pandas as pd
import os, subprocess, shutil

default_blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

def makeblast_db(seqs_path: str, db_out: str, show_command: bool = False, log_file: str = None, remove_old_db = False) -> str:
    # Remove the output directory if it exists
    if remove_old_db and os.path.exists(db_out):
        shutil.rmtree(db_out)
    elif os.path.exists(db_out):
        print(f'Warning output db directory already exist: {db_out}')


    # Create a fresh output directory
    os.makedirs(db_out, exist_ok=True)

    # Set the full path to the database output (same basename as input FASTA)
    db_file = os.path.join(db_out, os.path.basename(seqs_path))

    # Construct the makeblastdb command
    command = f"makeblastdb -in {seqs_path} -dbtype nucl -out {db_file}"
    
    if show_command:
        print("Running:", command)
    
    # Open the log file if provided
    log = open(log_file, "a") if log_file else None
    try:
        subprocess.run(command, shell=True, check=True, stdout=log, stderr=log)
        print('makeblastdb ran successfully.')
    except subprocess.CalledProcessError as e:
        print('Error in running makeblastdb:', e)
        raise
    finally:
        if log:
            log.close()
    
    return db_file


def blastn(blast_input_seqs: str, blast_db_file: str, blast_output_table_tsv: str,
           num_threads: int = 10, reward: int = 1, gap_extend: int = 2, gapopen: int = 5, penalty: int = -2, 
           word_size: int = 25, show_command: bool = True, log_file: str = None) -> None:
    command = (
        f"blastn -query {blast_input_seqs} -db {blast_db_file} -out {blast_output_table_tsv} "
        f"-num_threads {num_threads} -reward {reward} -gapextend {gap_extend} -gapopen {gapopen} -penalty {penalty} "
        f"-word_size {word_size} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
    )

    if show_command:
        print("Running:", command)

    log = open(log_file, "a") if log_file else None
    try:
        subprocess.run(command, shell=True, check=True, stdout=log, stderr=log)
        print('blastn ran successfully.')
    except subprocess.CalledProcessError as e:
        print('Error in running blastn:', e)
    finally:
        if log:
            log.close()


def blastn_subject(blast_input_seqs: str, blast_subject_seqs: str, blast_output_table_tsv: str,
           reward: int = 1, gap_extend: int = 2, gapopen: int = 5, penalty: int = -2, 
           word_size: int = 25, show_command: bool = True, log_file: str = None) -> None:
    command = (
        f"blastn -query {blast_input_seqs} -subject {blast_subject_seqs} -out {blast_output_table_tsv} "
        f"-reward {reward} -gapextend {gap_extend} -gapopen {gapopen} -penalty {penalty} "
        f"-word_size {word_size} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
    )

    if show_command:
        print("Running:", command)

    log = open(log_file, "a") if log_file else None
    try:
        subprocess.run(command, shell=True, check=True, stdout=log, stderr=log)
        print('blastn_subject ran successfully.')
    except subprocess.CalledProcessError as e:
        print('Error in running blastn_subject:', e)
    finally:
        if log:
            log.close()


def alignment_absolute_start_end(blast_df: pd.DataFrame) -> pd.DataFrame:
    """
    Process a BLAST output table loaded as a pandas DataFrame with the columns:
    'qstart', 'qend', 'sstart', and 'send'.

    Creates the columns: 'a.qstart', 'a.qend', 'a.sstart', 'a.send' to ensure that the
    start positions are always smaller than the end positions for each alignment.

    Also adds a 'dir' column to indicate the direction of the subject alignment:
    'F' for forward (sstart < send) and 'R' for reverse (sstart > send).
    """
    # Create a copy to avoid modifying the original slice
    blast_df = blast_df.copy()

    # Ensure the absolute start and end positions for query and subject alignments
    blast_df['a.qstart'] = np.minimum(blast_df['qstart'], blast_df['qend'])
    blast_df['a.qend']   = np.maximum(blast_df['qstart'], blast_df['qend'])
    blast_df['a.sstart'] = np.minimum(blast_df['sstart'], blast_df['send'])
    blast_df['a.send']   = np.maximum(blast_df['sstart'], blast_df['send'])

    # Determine alignment direction based on subject coordinates
    blast_df['dir'] = np.where(blast_df['sstart'] < blast_df['send'], 'F', 'R')

    return blast_df


# Check if there are differences when ordering by a.qstart, or by a.start
def total_coverage(cols: pd.DataFrame, use: str) -> float:
    
    """
    Use as input a blast output table loaded as a pandas DataFrame sorted by the column a.qstart (or a.sstart). 
    Also asumes the table has one subject, or is grouped by its subjects. 
    Returns a series with the total coverage of each subject by its queries.
    """

    # Reset index to avoid multi-level index issues
    cols = cols.reset_index(drop=True)
    
    # Handle query vs. subject-based coverage
    if use == 'query':
        sstart = cols['qstart'].values
        send = cols['qend'].values
        total_len = cols['qlen'].iloc[0]  # Assuming same qlen for all rows
    elif use == 'subject':
        sstart = cols['a.sstart'].values
        send = cols['a.send'].values
        total_len = cols['slen'].iloc[0]  # Assuming same slen for all rows
    else:
        raise ValueError(f"Unknown parameter 'use': {use}")
    
    # Initialize variables to track total coverage
    total_cover = 0
    current_start, current_end = sstart[0], send[0]
    
    # Loop through each interval and merge overlapping ones
    for i in range(1, len(sstart)):
        if sstart[i] <= current_end:  # Overlapping interval
            current_end = max(current_end, send[i])  # Merge intervals
        else:
            # Add the length of the current merged interval to total coverage
            total_cover += current_end - current_start
            # Start a new interval
            current_start, current_end = sstart[i], send[i]
    
    # Add the last interval
    total_cover += current_end - current_start
    
    # Calculate the total length of the query
    total_len = cols['slen'].iloc[0]  # Assuming qlen is the same for all rows
    
    # Return the fraction of coverage
    return total_cover / total_len

def total_coverage_temp(cols: pd.DataFrame, use: str) -> float:
    
    """
    Use as input a blast output table loaded as a pandas DataFrame sorted by the column a.qstart (or a.sstart). 
    Also asumes the table has one subject, or is grouped by its subjects. 
    Returns a series with the total coverage of each subject by its queries.
    """
    
    # Reset index to avoid multi-level index issues
    cols = cols.reset_index(drop=True)
    
    # Handle query vs. subject-based coverage
    if use == 'query':
        sstart = cols['a.qstart'].values
        send = cols['a.qend'].values
        total_len = cols['qlen'].iloc[0]  # Assuming same qlen for all rows
    elif use == 'subject':
        sstart = cols['a.sstart'].values
        send = cols['a.send'].values
        total_len = cols['slen'].iloc[0]  # Assuming same slen for all rows
    else:
        raise ValueError(f"Unknown parameter 'use': {use}")

    # Initialize coverage tracking
    total_cover = 0
    current_start, current_end = sstart[0], send[0]

    # Loop to merge overlapping intervals
    for i in range(1, len(sstart)):
        if sstart[i] <= current_end:
            current_end = max(current_end, send[i])  # Extend interval
        else:
            total_cover += current_end - current_start  # Add interval length
            current_start, current_end = sstart[i], send[i]  # Start new interval

    # Add the last interval
    total_cover += current_end - current_start

    # Calculate and return the fraction of coverage
    return total_cover / total_len

import numpy as np
import pandas as pd

def total_coverage_merged(cols: pd.DataFrame, use: str, inclusive: bool = True) -> float:
    """
    Compute fraction of covered bases by merging overlapping HSP intervals.

    Parameters
    ----------
    cols : DataFrame
        BLAST(-like) hits for a single subject or single query (already grouped).
        Accepts either 'qstart/qend' or 'a.qstart/a.qend' style; same for 'sstart/send'.
    use : {'query','subject'}
        Measure coverage on the query or the subject axis.
    inclusive : bool
        If True, treat coordinates as 1-based inclusive (BLAST style): length = end - start + 1.
        If False, treat as half-open: length = end - start.

    Returns
    -------
    float
        Fraction of covered length in [0,1]. Returns 0.0 for empty input.
    """
    if cols is None or len(cols) == 0:
        return 0.0

    # Normalize column names (support with/without 'a.' prefix)
    def col(*names):
        for n in names:
            if n in cols.columns:
                return n
        raise KeyError(f"None of the columns {names} found in DataFrame.")

    if use == "query":
        start_col = col("qstart", "a.qstart")
        end_col   = col("qend",   "a.qend")
        len_col   = col("qlen")
    elif use == "subject":
        start_col = col("sstart", "a.sstart")
        end_col   = col("send",   "a.send")
        len_col   = col("slen")
    else:
        raise ValueError(f"Unknown parameter 'use': {use}")

    # Extract and normalize intervals
    starts = cols[start_col].to_numpy(dtype=np.int64, copy=False)
    ends   = cols[end_col].to_numpy(dtype=np.int64, copy=False)

    # Ensure start <= end per HSP (handles reverse-strand hits)
    s = np.minimum(starts, ends)
    e = np.maximum(starts, ends)

    # Length convention
    add_one = 1 if inclusive else 0

    # Sort by start
    order = np.argsort(s, kind="mergesort")
    s, e = s[order], e[order]

    # Optionally clip to the sequence length (defensive)
    total_len = int(cols[len_col].iloc[0])
    if total_len <= 0:
        return 0.0
    s = np.clip(s, 1, total_len)
    e = np.clip(e, 1, total_len)

    # Merge overlaps
    total_cover = 0
    cur_s, cur_e = int(s[0]), int(e[0])

    for i in range(1, len(s)):
        if s[i] <= cur_e:  # overlap or touching
            if e[i] > cur_e:
                cur_e = int(e[i])
        else:
            total_cover += (cur_e - cur_s + add_one)
            cur_s, cur_e = int(s[i]), int(e[i])

    # Add last interval
    total_cover += (cur_e - cur_s + add_one)

    # Fraction
    frac = total_cover / total_len
    # Keep it in [0,1] just in case
    return max(0.0, min(1.0, float(frac)))


def select_internal_aligns(df: pd.DataFrame, border: int) -> pd.DataFrame:
    """
    Selects BLAST alignments where both the query and subject alignments are 
    internal, based on a user-defined border.
    
    Parameters:
        df (pd.DataFrame): The DataFrame containing BLAST alignment results.
        border (int): The border threshold for filtering.
        
    Returns:
        pd.DataFrame: Filtered DataFrame with only internal alignments.
    """
    
    # Create Boolean masks for the conditions
    query_internal = (df['qstart'] > border) & (df['qend'] < (df['qlen'] - border))
    subject_internal = (df['a.sstart'] > border) & (df['a.send'] < (df['slen'] - border))
    
    # Apply both masks to filter the DataFrame
    internal_df = df[query_internal & subject_internal]
    
    return internal_df


import pandas as pd

def fuse_contained_intervals(
    df: pd.DataFrame,
    start_col: str = 'qstart',
    end_col: str   = 'qend',
    group_cols: list | None = None
) -> pd.DataFrame:
    """
    Remove intervals that are fully contained within larger ones,
    keeping only the 'outermost' alignments per group.

    Parameters:
    -----------
    df : pd.DataFrame
        Your BLAST (or other) alignment table.
    start_col : str
        Name of the start‐coordinate column (default: 'qstart').
    end_col : str
        Name of the end‐coordinate column (default: 'qend').
    group_cols : list of str, optional
        Columns to group by. If None, defaults to ['qseqid','sseqid']
        (but only uses those that actually exist in `df`).

    Returns:
    --------
    pd.DataFrame
        A filtered DataFrame containing only the non‐contained alignments,
        sorted by the group columns and start position.
    """
    # Handle empty
    if df.empty:
        return df.copy()

    # Pick sensible defaults for grouping
    if group_cols is None:
        group_cols = [c for c in ('qseqid','sseqid') if c in df.columns]
    req = set(group_cols) | {start_col, end_col}
    missing = req - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing!r}")

    # Sort so that within each group, longer (outer) intervals come first:
    #  - same start → larger end first
    sorted_df = df.sort_values(
        by=group_cols + [start_col, end_col],
        ascending=[True]*len(group_cols) + [True, False],
        kind="mergesort"
    )

    keep_idx: list[int] = []
    # For each group, walk through sorted intervals and drop those contained
    for _, grp in sorted_df.groupby(group_cols, sort=False):
        seen: list[tuple[int,int]] = []
        for idx, row in grp.iterrows():
            s, e = row[start_col], row[end_col]
            # if it's contained in any "seen" interval, skip it
            if any(s >= ss and e <= ee for ss, ee in seen):
                continue
            # otherwise keep it, and record its bounds
            keep_idx.append(idx)
            seen.append((s, e))

    # Rebuild result, sorted by group + start
    result = df.loc[keep_idx].copy()
    result = result.sort_values(by=group_cols + [start_col], kind="mergesort")
    return result.reset_index(drop=True)


#from intervaltree import IntervalTree

import pandas as pd

def fuse_contained_intervals(
    df: pd.DataFrame,
    start_col: str,
    end_col: str,
    group_cols: list[str] | None = None,
    keep_highest_score_col: str | None = None
) -> pd.DataFrame:
    """
    Removes rows where the interval [start_col, end_col] is fully contained in another,
    within each group defined by `group_cols`.

    Parameters
    ----------
    df : pd.DataFrame
        The input table (e.g. BLAST hits or genomic alignments).
    start_col : str
        Column name for interval start.
    end_col : str
        Column name for interval end.
    group_cols : list of str or None
        Grouping columns. If None, uses ['qseqid', 'sseqid'] if they exist.
    keep_highest_score_col : str or None
        If set, breaks ties by keeping the row with the highest value in this column (e.g. 'bitscore').

    Returns
    -------
    pd.DataFrame
        A filtered DataFrame with redundant (contained) intervals removed.
    """
    if df.empty:
        return df.copy()

    # Infer group columns if not provided
    if group_cols is None:
        group_cols = [col for col in ['qseqid', 'sseqid'] if col in df.columns]

    required = set([start_col, end_col] + group_cols)
    if missing := required - set(df.columns):
        raise KeyError(f"Missing columns: {missing}")

    # Determine sort order
    sort_cols = group_cols + [start_col, end_col]
    ascending = [True] * len(group_cols) + [True, False]
    if keep_highest_score_col:
        sort_cols.insert(len(group_cols), keep_highest_score_col)
        ascending.insert(len(group_cols), False)

    sorted_df = df.sort_values(by=sort_cols, ascending=ascending, kind='mergesort')

    keep_indices = []

    for _, group in sorted_df.groupby(group_cols, sort=False):
        seen_intervals = []
        for idx, row in group.iterrows():
            s, e = row[start_col], row[end_col]
            if any(s >= ss and e <= ee for ss, ee in seen_intervals):
                continue
            keep_indices.append(idx)
            seen_intervals.append((s, e))

    return df.loc[keep_indices].sort_values(by=sort_cols).reset_index(drop=True)


def jumping_windows(coords_units: dict, radius=15, q=0) -> dict:
    """A recursive function that adds neighborhood coordinate appearances 
    to the coordinate with most appearances within a radius.
    
    Args:
        coords_units (dict): A dictionary where keys are coordinates (int) 
                             and values are counts (int).
        radius (int, optional): Search radius for neighboring coordinates. Defaults to 3.
        q (int, optional): Internal variable, used in recursion. Defaults to 0.
    """
    
    # End function if q is greater than the size of coords_units
    if q >= len(coords_units.keys()):
        return coords_units
    
    # Sort coordinates by their number of appearances (value), descending
    coords_order = sorted(coords_units.keys(), reverse=True)
    
    # Check for coordinates in neighborhood of current q
    for r in range(-radius, radius + 1):
        if r != 0:
            current_coord = coords_order[q]
            neighbor_coord = current_coord + r
            # Add neighbor counts to current coord
            coords_units[current_coord] += coords_units.get(neighbor_coord, 0)
            # Remove neighbor coord if it exists
            coords_units.pop(neighbor_coord, None)
    
    # Recurse to next coordinate
    q += 1
    coords_units = jumping_windows(coords_units, radius=radius, q=q)
    return coords_units



def unnoise_coords( df: pd.DataFrame, radious: int = 15) -> dict: 

    coords = df['qstart'].to_list()
    coords.extend(df['qend'])
    coords.extend([1,df['qlen'].iloc[0]])

    coords_df = pd.DataFrame(data={ 'coords': coords, 'units' : np.ones(len(coords))})
    coords_df = coords_df.groupby('coords').sum().sort_values('units', ascending=False)        
    coords_dict = coords_df['units'].to_dict()

    unnoised_coords = jumping_windows(coords_dict, radious)
    unnoised_df = pd.DataFrame.from_dict(unnoised_coords, orient='index', columns=['instances'])
    unnoised_df = unnoised_df.reset_index().rename(columns={'index' : 'coords'})

    return(unnoised_df)


def retrieve_single_alignments(blast_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Count the number of alignments per query
    counts = blast_df.groupby('qseqid')['sseqid'].count()
    
    # Identify those with only one match
    single_hit_qseqids = counts[counts == 1].index

    # Retrieve alignments with only one match
    retrieve_df = blast_df[blast_df['qseqid'].isin(single_hit_qseqids)]

    # Retrieve alignments with more than one match
    filtered_df = blast_df[~blast_df['qseqid'].isin(single_hit_qseqids)]

    return retrieve_df, filtered_df

def retrieve_nhits_alignments(blast_df: pd.DataFrame, n: int = 10) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Count the number of alignments per query
    counts = blast_df.groupby('qseqid')['sseqid'].count()
    
    # Identify those with only one match
    single_hit_qseqids = counts[counts >= n].index

    # Retrieve alignments with only one match
    retrieve_df = blast_df[blast_df['qseqid'].isin(single_hit_qseqids)]

    # Retrieve alignments with more than one match
    filtered_df = blast_df[~blast_df['qseqid'].isin(single_hit_qseqids)]

    return retrieve_df, filtered_df


from Bio import SeqIO
import pandas as pd

def find_unmapped(df: pd.DataFrame, fasta_queries: str, fasta_subjects: str):
    """
    Identify which query and subject sequences had no BLASTn hits.

    Parameters
    ----------
    df : pd.DataFrame
        BLASTn output table with at least 'qseqid' and 'sseqid' columns.
    fasta_queries : str
        Path to the FASTA file of query sequences.
    fasta_subjects : str
        Path to the FASTA file of subject sequences.

    Returns
    -------
    unmapped_queries : set of str
        IDs from the query FASTA that do not appear in df['qseqid'].
    unmapped_subjects : set of str
        IDs from the subject FASTA that do not appear in df['sseqid'].
    """
    # Read all query IDs
    query_ids = {record.id for record in SeqIO.parse(fasta_queries, "fasta")}
    # Read all subject IDs
    subject_ids = {record.id for record in SeqIO.parse(fasta_subjects, "fasta")}

    # IDs that had at least one alignment
    aligned_queries = set(df['qseqid'])
    aligned_subjects = set(df['sseqid'])

    # Compute the unmapped sets
    unmapped_queries = query_ids - aligned_queries
    unmapped_subjects = subject_ids - aligned_subjects

    return unmapped_queries, unmapped_subjects

# Example usage:
# import pandas as pd
# df = pd.read_csv("blast_results.tsv", sep="\t", header=None,
#                  names=[...,'qlen','slen'])
# unm_q, unm_s = find_unmapped(df, "queries.fasta", "subjects.fasta")
# print("Unmapped queries:", unm_q)
# print("Unmapped subjects:", unm_s)



from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt

def blastn_dotplot(
    fasta1: str,
    fasta2: str,
    output_base_dir: str,
    word_size: int = 20,
    qnames: bool = False,
    snames: bool = False,
    annot_size: float = 12,
    title: str | None = "Blastn dotplot",
    work_name: str | None = None,
) -> pd.DataFrame:
    """
    Run two BLASTn searches (auto vs. auto, nam vs. auto), load results,
    compute cumulative positions, and plot a dot-plot of alignments.

    Parameters
    ----------
    work_name : str
        Identifier for naming output files.
    fasta1 : str
        Path to a FASTA file with one or more sequences.
    fasta2 : str
        Path to a FASTA file with one or more sequences.
    output_base_dir : str
        Base directory under which BLAST output TSVs and visuals will be saved.
    word_size_auto : int, optional
        BLASTn word size for the auto-auto comparison (default=20).
    word_size_nam : int, optional
        BLASTn word size for the Nam-vs-auto comparison (default=20).
    """
    # Prepare output paths
    base = Path(output_base_dir)
    if work_name:
        workdir = base / work_name
    else:
        workdir = base
        
    workdir.mkdir(exist_ok=True, parents=True)
    blast_tsv  = workdir / f"blastn.fasta1_vs_fasta2.tsv"

    # Run BLASTn
    blastn_subject(
        blast_input_seqs=fasta1,
        blast_subject_seqs=fasta2,
        blast_output_table_tsv=blast_tsv,
        word_size=word_size
    )

    # Load results
    df = pd.read_csv(blast_tsv, sep='\t', header=None)
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

    # Sort and compute absolute start/end
    df['minus_qstart'] = -df['qstart']
    df = df.sort_values(['qlen', 'minus_qstart'], ascending=False)
    df = alignment_absolute_start_end(df)

    # Compute cumulative offsets for each sequence
    q_order = df['qseqid'].unique()
    s_order = df['sseqid'].unique()
    qlen_mean = df.groupby('qseqid')['qlen'].mean().reindex(q_order)
    slen_mean = df.groupby('sseqid')['slen'].mean().reindex(s_order)

    cumsum_q = qlen_mean.cumsum().shift(1).fillna(0).reset_index()
    cumsum_s = slen_mean.cumsum().shift(1).fillna(0).reset_index()
    cumsum_q.columns = ['qseqid', 'cumsum_qlen']
    cumsum_s.columns = ['sseqid', 'cumsum_slen']

    df = df.merge(cumsum_q, on='qseqid').merge(cumsum_s, on='sseqid')

    # Prepare plot scaling
    total_q = cumsum_q['cumsum_qlen'].iat[-1]
    total_s = cumsum_s['cumsum_slen'].iat[-1]
    fig_scale_q = total_q / (total_q + total_s)
    fig_scale_s = total_s / (total_q + total_s)

    # Add cumulative positions to endpoints
    df['x_start'] = df['cumsum_qlen'] + df['qstart']
    df['x_end']   = df['cumsum_qlen'] + df['qend']
    df['y_start'] = df['cumsum_slen'] + df['sstart']
    df['y_end']   = df['cumsum_slen'] + df['send']

    # Plot
    fig, ax = plt.subplots(figsize=(24 * fig_scale_q, 24 * fig_scale_s))
    ax.set_xlim(0, df['x_end'].max())
    ax.set_ylim(0, df['y_end'].max())

    for _, row in df.iterrows():
        ax.plot(
            [row['x_start'], row['x_end']],
            [row['y_start'], row['y_end']],
            linestyle='-', color='blue'
        )

    # Sequence boundaries
    for pos in df['cumsum_qlen'].unique():
        ax.axvline(pos, linestyle='--', linewidth=0.5, color='blue')
    for pos in df['cumsum_slen'].unique():
        ax.axhline(pos, linestyle='--', linewidth=0.5, color='red')

    if qnames:
        # Extract one row per query sequence: cumsum start + its length
        qs = (
            df[['qseqid', 'cumsum_qlen', 'qlen']]
            .drop_duplicates(subset='qseqid')
            .reset_index(drop=True)
        )
            
        # Annotate each query at the midpoint of its block
        for _, row in qs.iterrows():
            start = row.cumsum_qlen
            length = row.qlen
            midpoint = start + length / 2

            ax.annotate(
                row.qseqid,
                xy=(midpoint, 0),
                xytext=(0, -30),              # offset below the axis
                textcoords='offset points',
                rotation=90,                  # vertical text
                va='top',
                ha='center',
                fontsize=annot_size
            )
        
    if snames:                
        # Extract one row per query sequence: cumsum start + its length
        ss = (
            df[['sseqid', 'cumsum_slen', 'slen']]
            .drop_duplicates(subset='sseqid')
            .reset_index(drop=True)
        )


        # Annotate each query at the midpoint of its block
        for _, row in ss.iterrows():
            start = row.cumsum_slen
            length = row.slen
            midpoint = start + length / 2

            ax.annotate(
                row.sseqid,
                xy=(0, midpoint),
                xytext=(-100, 0),              # offset below the axis
                textcoords='offset points',
                rotation=0,                  # vertical text
                va='top',
                ha='center',
                fontsize=annot_size
                
            )
    
    # Adjust axis labels with padding
    ax.set_xlabel('Query Sequence', labelpad=0)
    ax.set_ylabel('Subject Sequence Position')
    ax.set_title(f'Dot Plot of BLASTN Alignments: {title}')
    plt.tight_layout()
    plt.savefig( workdir / 'dotplot.svg')
    
    print(f"Dotplot Saved in: {workdir / 'dotplot.svg'}")
    plt.show()
    return(df)


import pandas as pd

def unmapped_subject_regions(
    df: pd.DataFrame,
    subject_id_col: str = "sseqid",
    subject_len_col: str = "slen",
    sstart_col: str = "a.sstart",   # use "sstart" if that's your column
    send_col: str = "a.send",       # use "send" if that's your column
    min_pident: float | None = None,
    min_aln_len: int | None = None,
    pident_col: str = "pident",
    length_col: str = "length",
    merge_touching: bool = True,    # treat [1,10] and [11,20] as continuous
    pad: int = 0                    # expand mapped regions by N bases before complement
) -> pd.DataFrame:
    """
    Return a table of unmapped intervals per subject.
    Output columns: [sseqid, gap_start, gap_end, gap_len, n_mapped_intervals]
    Coordinates are 1-based inclusive like BLAST.
    """
    df = df.copy()

    # Optional filtering
    if (min_pident is not None) and (pident_col in df):
        df = df[df[pident_col] >= min_pident]
    if (min_aln_len is not None) and (length_col in df):
        df = df[df[length_col] >= min_aln_len]

    # Normalize subject intervals (handle minus strand)
    s_lo = df[[sstart_col, send_col]].min(axis=1)
    s_hi = df[[sstart_col, send_col]].max(axis=1)
    if pad > 0:
        # We'll clamp later using slen
        df["_s_lo"] = s_lo - pad
        df["_s_hi"] = s_hi + pad
    else:
        df["_s_lo"] = s_lo
        df["_s_hi"] = s_hi

    # We need subject lengths. If they are not present per-row, merge them in.
    # Expect either a per-row 'slen' or an auxiliary table you can join beforehand.
    if subject_len_col not in df.columns:
        raise ValueError(f"Column '{subject_len_col}' with subject lengths is required.")

    out_rows = []

    for sseqid, sub in df.groupby(subject_id_col, sort=False):
        slen = int(sub[subject_len_col].iloc[0])

        # Collect and clamp intervals to [1, slen]
        ivals = []
        for lo, hi in zip(sub["_s_lo"], sub["_s_hi"]):
            lo = max(1, int(lo))
            hi = min(slen, int(hi))
            if lo <= hi:
                ivals.append((lo, hi))

        if not ivals:
            # No mapped region => whole subject is unmapped
            out_rows.append({
                subject_id_col: sseqid,
                "gap_start": 1,
                "gap_end": slen,
                "gap_len": slen,
                "n_mapped_intervals": 0
            })
            continue

        # Sort & merge mapped intervals
        ivals.sort(key=lambda x: x[0])
        merged = []
        cur_lo, cur_hi = ivals[0]
        for lo, hi in ivals[1:]:
            # If touching is considered merged, use lo <= cur_hi + 1
            if (lo <= cur_hi + (1 if merge_touching else 0)):
                cur_hi = max(cur_hi, hi)
            else:
                merged.append((cur_lo, cur_hi))
                cur_lo, cur_hi = lo, hi
        merged.append((cur_lo, cur_hi))

        # Complement: gaps between [1, slen]
        # gap before first
        if merged[0][0] > 1:
            gs, ge = 1, merged[0][0] - 1
            out_rows.append({subject_id_col: sseqid, "gap_start": gs, "gap_end": ge,
                             "gap_len": ge - gs + 1, "n_mapped_intervals": len(merged)})

        # gaps between merged intervals
        for (l1, r1), (l2, r2) in zip(merged, merged[1:]):
            if l2 > r1 + 1:
                gs, ge = r1 + 1, l2 - 1
                out_rows.append({subject_id_col: sseqid, "gap_start": gs, "gap_end": ge,
                                 "gap_len": ge - gs + 1, "n_mapped_intervals": len(merged)})

        # gap after last
        if merged[-1][1] < slen:
            gs, ge = merged[-1][1] + 1, slen
            out_rows.append({subject_id_col: sseqid, "gap_start": gs, "gap_end": ge,
                             "gap_len": ge - gs + 1, "n_mapped_intervals": len(merged)})

    gaps_df = pd.DataFrame(out_rows).sort_values([subject_id_col, "gap_start"], ignore_index=True)
    return gaps_df


#def fuse_contained_intervals_fast(df: pd.DataFrame(), start_col='qstart', end_col='qend'),:
#    """
#    Uses an interval tree to efficiently remove contained alignments.
#    """
#    if df.empty:
#        return df.copy()
#
#    tree = IntervalTree()
#    kept_rows = []
#
#    # Sort rows by length descending (outer intervals come first)
#    df_sorted = df.copy()
#    df_sorted['len'] = df_sorted[end_col] - df_sorted[start_col]
#    df_sorted = df_sorted.sort_values(by='len', ascending=False)
#
#    df_size = df_sorted.size[0]
#    last_printed = -1
#    for i, (_, row) in enumerate(df_sorted.iterrows(), start=1):
#        s, e = row[start_col], row[end_col]
#        
#        percent = int(i / df_size * 100)
#        if percent != last_printed:
#            print(f"{percent}% done")
#            last_printed = percent
#
#        
#        # Check if this interval is contained in any existing one
#        containing = tree.envelop(s, e)
#        if not containing:
#            kept_rows.append(row)
#            tree.addi(s, e)
#
#    return pd.DataFrame(kept_rows).drop(columns='len').reset_index(drop=True)
#