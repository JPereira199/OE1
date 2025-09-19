
from __future__ import annotations
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

################################## Convert list with blocks to ABC format ###################### 

def read_list_to_abc( df: pd.DataFrame, example = False ) -> pd.DataFrame:

    # Create the initial DataFrame
    if example:
        df = pd.DataFrame({'string_code': ['ABBAAA', 'ABSBACA', 'ABBACA', 'ABBACA']})

    # Convert the 'col' column into lists of letters
    df['letters'] = df['string_code'].apply(list)

    # Explode the 'letters' column to create a row for each letter
    df_exploded = df.explode('letters')

    # Reset the index to turn the index into a column
    df_exploded = df_exploded.reset_index()

    # Rename 'index' to 'original_index' for clarity
    df_exploded = df_exploded.rename(columns={'index': 'group'})

    # Select and reorder the necessary columns
    df_exploded = df_exploded[['group', 'letters']]

    df_results_new = pd.DataFrame()
    df_results_new['prev_element'] = df_exploded.groupby(['group'])['letters'].shift(1, fill_value='Start')
    df_results_new['next_element'] = df_exploded.groupby(['group'])['letters'].shift(-1, fill_value='End')

    df_results_new['count'] = 1
    df_results_new = df_results_new.groupby(['prev_element', 'next_element']).sum().reset_index()

    return(df_results_new)

################################### Make Directed graph from ABC input ############################

def abc_to_dir_graph( df: pd.DataFrame, graph_plot_output: str,adjacency_df_output: str, example = False ):

    # Create directed graph with filtered data
    G_filtered = nx.DiGraph()

    # Add edges to the filtered graph
    for _, row in df.iterrows():
        nodeA = row['node']
        nodeB = row['next_node']
        weight = row['count']
        # Add edge with weight to the G_filtered graph
        if weight > 1: 
            G_filtered.add_edge(nodeA, nodeB, weight=weight)

    # Use spring layout for positioning
    pos = nx.spring_layout(G_filtered, k=0.5)  # Adjust k for spacing

    plt.figure(figsize=(12, 10))

    # Compute weighted degree for node size (sum of weights for incoming and outgoing edges)
    node_weighted_degree = {
        node: sum(weight['weight'] for _, _, weight in G_filtered.edges(node, data=True)) +
              sum(weight['weight'] for _, _, weight in G_filtered.in_edges(node, data=True))
        for node in G_filtered.nodes()
    }

    # Set node size based on weighted degree, scale it appropriately
    node_size = [node_weighted_degree[node] * 0.5 for node in G_filtered]

    # Draw the nodes
    nx.draw_networkx_nodes(G_filtered, pos, node_size=node_size, node_color='skyblue', alpha=0.7)

    # Draw the edges
    edge_weights = [G_filtered[u][v]['weight'] for u, v in G_filtered.edges()]
    nx.draw_networkx_edges(G_filtered, pos, width=[weight * 0.05 for weight in edge_weights], alpha=0.6)

    # Draw labels
    nx.draw_networkx_labels(G_filtered, pos, font_size=10, font_color='black', font_weight='bold')

    plt.title("Filtered Directed Graph with Node Size Based on Weighted Degree")
    plt.savefig(graph_plot_output)

    # Get a sorted list of nodes to maintain consistent ordering
    nodes = sorted(G_filtered.nodes())
    
    #################################### Adjacency Matrix ####################################

    # Extract the adjacency matrix as a pandas DataFrame
    adjacency_df = nx.to_pandas_adjacency(G_filtered, nodelist=nodes, weight='weight')

    #print("\nAdjacency matrix as pandas DataFrame:")
    #print(adjacency_df)

    #adjacency_df_output = '/home/jpereira/OEs/OE1/Data_output/adjacency_matrix_new.tsv'
    adjacency_df.to_csv(adjacency_df_output, sep = '\t')
    
    return(adjacency_df)

# Define the blocks_encoding function
def blocks_encoding(group: pd.DataFrame):
    # Concatenate the 'encode' column for the group
    string_code = ''.join(group['encode'])
    return pd.Series({'string_code': string_code})



######################################## Seed MSA and extension ########################################



from pathlib import Path
from io import StringIO
import subprocess
import tempfile
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def align_seed_neighbours(
    reads_fasta: str | Path,
    seed_hits: pd.DataFrame,
    s_extension: int = 100,
    e_extension: int = 50,
    *,
    mafft_bin: str = "mafft",
    max_iter: Optional[int] = None,
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    threads: int = 10,
    alphabet: str = "acgt-",
    print_fragments: bool = False,
    num_prints: int = 10,
    keep_length: bool = False,
    show_cmd: bool = False
) -> Tuple[List[SeqRecord], pd.DataFrame, pd.DataFrame]:
    """
    Extract flanking regions around seed matches, align them with MAFFT,
    and return the alignment plus per-column base frequencies and a consensus
    table.

    Parameters
    ----------
    reads_fasta : str or Path
        FASTA file with raw reads.
    seed_hits : pd.DataFrame
        BLAST-like table with at least the columns
        ['sseqid', 'sstart', 'send', 'slen'].
    s_extension, e_extension : int
        Bases to extend upstream/downstream from the seed match.
    mafft_bin : str, default "mafft"
        MAFFT executable.
    max_iter : int or None, default None
        If given, passed to MAFFT as ``--maxiterate <max_iter>``.
    gap_open : float or None, default None
        Gap-opening penalty (``--op``).  MAFFT default is 1.53.
    gap_extend : float or None, default None
        Gap-extension penalty (``--ep``).  MAFFT default is 0.0.
    alphabet : str, default "acgt-"
        Symbols counted when computing frequencies.

    Returns
    -------
    aligned : list[SeqRecord]
        MAFFT-aligned sequences.
    df_freq : pd.DataFrame
        Per-column symbol frequencies (1-based index).
    df_cons : pd.DataFrame
        Consensus symbol and its frequency per column.
    """
    # ---------- 1. Load relevant reads ----------
    seed_ids = set(seed_hits["sseqid"])
    reads = {
        rec.id: rec
        for rec in SeqIO.parse(str(reads_fasta), "fasta")
        if rec.id in seed_ids
    }
    if not reads:
        raise ValueError("None of the seed IDs were found in the FASTA file.")

    # ---------- 2. Slice neighbouring fragments ----------
    neighbours: list[SeqRecord] = []
    over_extended_start = [] # Used to determine the real aligment size
    over_extended_end = []
    for read_id in seed_ids:
        for _, row in seed_hits.loc[seed_hits["sseqid"] == read_id].iterrows():
            if row["sstart"] < row["send"]:
                
                start    = max(row["sstart"] - s_extension - 1, 0)
                end      = min(row["send"] + e_extension, row["slen"])
                
                fragment = reads[read_id].seq[start:end]
                orientation = "F"
                
                # Add the surplus size if the reads were over extended 
                if start == 0:
                    over_extended_start.append(0 - (row["sstart"] - s_extension - 1) )
                else:
                    over_extended_start.append(0)
                if end == row['slen']:
                    over_extended_end.append( (row["send"] + e_extension) - row['slen'])
                else:
                    over_extended_end.append(0)

            else:
                start    = max(row["send"] - e_extension - 1, 0)
                end      = min(row["sstart"] + s_extension, row["slen"])
                
                fragment = reads[read_id].seq[start:end].reverse_complement()
                orientation = "R"
                
                # Add the surplus size if the reads were over extended 
                if start == 0:
                    over_extended_start.append(0 - (row["send"] - e_extension - 1))
                else:
                    over_extended_start.append(0)
                if end == row['slen']:
                    over_extended_end.append((row["sstart"] + s_extension) - row['slen'])
                else:
                    over_extended_end.append(0)
             
            neighbours.append(
                SeqRecord(
                    Seq(fragment),
                    id=f"{read_id}-{start}:{end}",
                    description=orientation,
                )
            )
    
    real_start_extension = s_extension - min(over_extended_start)
    real_end_extension   = e_extension - min(over_extended_end)
    
    real_start_extension_l = [s_extension - o for o in over_extended_start]
    real_end_extension_l   = [e_extension - o for o in over_extended_end]
    
    if print_fragments:
        for i, record in enumerate(neighbours[0:num_prints]):
            print(f"{i}_{record.description}: {record.seq}")

    if len(neighbours) < 2:
        raise ValueError("Need at least two neighbouring fragments for alignment.")

    # ---------- 3. Build MAFFT command ----------
    mafft_cmd: list[str] = [mafft_bin, "--auto", "--thread", str(threads)]
    if max_iter is not None:
        mafft_cmd += ["--maxiterate", str(max_iter)]
    if gap_open is not None:
        mafft_cmd += ["--op", str(gap_open)]
    if gap_extend is not None:
        mafft_cmd += ["--ep", str(gap_extend)]
    if keep_length:
        mafft_cmd += ["--keeplength", str(gap_extend)]
        
    # ---------- 4. Align ----------
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".fa", delete=False) as tmp:
        SeqIO.write(neighbours, tmp, "fasta")
        tmp.flush()
        try:
            result = subprocess.run(
                mafft_cmd + [tmp.name],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
            )
        finally:
            Path(tmp.name).unlink(missing_ok=True)
    
    if show_cmd:
        print(f"Mafft Commmad {' '.join(mafft_cmd)}")


    aligned = list(SeqIO.parse(StringIO(result.stdout), "fasta"))

    # ---------- 5. Column-frequency table ----------
    arr = np.array([list(str(rec.seq).lower()) for rec in aligned], dtype="U1")
    symbols = np.array(list(alphabet))
    counts = {s: (arr == s).sum(axis=0) for s in symbols}

    df_freq = pd.DataFrame(counts, index=np.arange(1, arr.shape[1] + 1))
    df_freq = df_freq.div(df_freq.sum(axis=1).replace(0, np.nan), axis=0)

    # ---------- 6. Consensus ----------
    df_cons = pd.DataFrame(
        {"base": df_freq.idxmax(axis=1), "freq": df_freq.max(axis=1)}
    )

    return aligned, df_freq, df_cons, real_start_extension_l, real_end_extension_l



############################ Find block boundaries ######################################

# Split the function in two sections (one to find the start block pos, and other to find the end)
# Then, make the function recursive (or include a while loop), extending the the windows size
# Stop the function after X number of recursions or if the windows size exceed the start (or end)
# seed position. 
def find_block_boundaries(df_clean, start_signal_pos, end_signal_pos, start_seed_pos, end_seed_pos,
                          window_radius=15, threshold=0.10, max_windows_times=5):
    """
    Find the start and end positions of a block around signal peaks.

    Parameters:
        df_clean (pd.DataFrame): DataFrame with at least a 'freq' column.
        start_signal_pos (int): Position of the start signal in original coordinates (before offset).
        end_signal_pos (int): Position of the end signal in original coordinates (before offset).
        window_radius (int): Window size used to define local frequency mean.
        threshold (float): Threshold below mean to define block boundaries.
        max_windows_times (int): Maximum number times the windows could be extended when searching for 
                                 the block postion

    Returns:
        tuple: (start_block_pos, end_block_pos)
    """
    
    # --- Start Block ---
    start_block_pos = None
    windows_times = 0
    while not start_block_pos:
        windows_times += 1
        if windows_times > max_windows_times:
            print(f"Warning!: Couldn't find a block START boundary even after extendid the block search to \
                  {max_windows_times} the size of the windows_radious ({window_radius})")
            break
        
        max_window_start_pos = start_signal_pos + (windows_times + 1) * window_radius
        #min_windows_start_pos = start_signal_pos +  window_radius * windows_times
        if max_window_start_pos > end_seed_pos:
            print(f"Warning!: Maximum windows end position ({max_window_start_pos}) exceded the original seed \
                    position ({end_seed_pos}) when searching for the START block boundary")
            break

        start_local_freq_mean = df_clean['freq'][
             max_window_start_pos - window_radius : max_window_start_pos 
        ].mean()

        for i, row in df_clean.iloc[: max_window_start_pos - window_radius][::-1].iterrows():
            if row['freq'] < start_local_freq_mean - threshold:
                start_block_pos = i + 1
                break

    # --- End Block ---
    end_block_pos = None
    windows_times = 0
    while not end_block_pos:
        windows_times += 1
        if windows_times > max_windows_times:
            print(f"Warning!: Couldn't find a block END boundary even after extendid the block search to \
                  {max_windows_times} the size of the windows_radious ({window_radius})")
            break
        
        #max_window_start_pos = end_signal_pos -  window_radius * windows_times
        min_windows_start_pos = end_signal_pos -  window_radius * (windows_times + 1)
        if min_windows_start_pos < start_seed_pos:
            print(f"Warning!: Minimum windows start position ({min_windows_start_pos}) exceded the original seed \
                    position ({start_seed_pos}) when searching for the END block boundary")
            break

        end_local_freq_mean = df_clean['freq'][
             min_windows_start_pos : min_windows_start_pos + window_radius
        ].mean()

        for i, row in df_clean.iloc[min_windows_start_pos + window_radius :].iterrows():
            if row['freq'] < end_local_freq_mean - threshold:
                end_block_pos = i - 1
                break

    return start_block_pos, end_block_pos

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from typing import Iterable, Optional, Union

def plot_windowed_metrics(
    *,
    x_vals,
    plot_means,
    plot_stds,
    plot_signal,
    block_num: int,
    start_seed_pos: int,
    end_seed_pos: int,
    max_indices: Iterable[int],
    # optional signal windows
    start_sigs=None,
    mean_start_sigs=None,
    end_sigs=None,
    mean_end_sigs=None,
    start_signal_pos: Optional[int] = None,
    mean_start_signal_pos: Optional[int] = None,
    end_signal_pos: Optional[int] = None,
    mean_end_signal_pos: Optional[int] = None,
    params_window_radius: int = 0,
    # layout / saving
    n: Optional[int] = None,
    visuals_dir: Union[str, Path, None] = None,
    filename: Optional[str] = None,
    save: bool = True,
    show: bool = True,
) -> Optional[Path]:
    """
    Plot windowed MSA metrics for a given block and (optionally) save the figure.

    Returns the saved Path if `save=True`, otherwise None.
    """

    def _nonempty(x):
        if x is None:
            return False
        if hasattr(x, "empty"):
            return not x.empty  # pandas Series/DataFrame
        try:
            return len(x) > 0   # lists/arrays
        except Exception:
            return True         # truthy fallback

    visuals_dir = Path(visuals_dir) if visuals_dir is not None else Path(".")
    if filename is None:
        filename = f"windowed_msa_metrics_block_{block_num}.svg"
    out_path = visuals_dir / filename

    plt.figure(figsize=(12, 5))
    plt.yticks(np.linspace(0, 1, 21).round(2))
    plt.xlabel("Alignment Position")
    plt.ylabel("Value")
    plt.title(f"Windowed Metrics Block: {block_num}")

    # main series
    plt.plot(x_vals, plot_means,  label="Mean", linewidth=2)
    plt.plot(x_vals, plot_stds,   label="Std Dev", linewidth=2)
    plt.plot(x_vals, plot_signal, label="Signal", linewidth=2)

    # seed window
    plt.axvline(start_seed_pos, color='gray', linestyle='-', linewidth=1.2, label='Seed Start')
    plt.axvline(end_seed_pos,   color='gray', linestyle='-', linewidth=1.2, label='Seed End')

    # highlight max-mean indices
    for idx in max_indices:
        plt.axvline(idx, color='green', linestyle='--', label='Max Means')

    # start thresholds
    if _nonempty(start_sigs):
        plt.axvline(start_signal_pos, color='red', linestyle='-', label='Signal Threshold')
        plt.axvline(start_signal_pos + params_window_radius, color='Blue', linestyle='-', label='Signal + Windows')
        plt.axvline(start_signal_pos + 2 * params_window_radius, color='Skyblue', linestyle='-', label='Signal + (2 x Windows)')
    elif _nonempty(mean_start_sigs):
        plt.axvline(mean_start_signal_pos, color='brown', linestyle='-', label='Mean Signal Threshold')
        plt.axvline(mean_start_signal_pos + params_window_radius, color='darkgoldenrod', linestyle='-', label='Mean Signal + Windows')
        plt.axvline(mean_start_signal_pos + 2 * params_window_radius, color='gold', linestyle='-', label='Mean Signal + (2 x Windows)')

    # end thresholds
    if _nonempty(end_sigs):
        plt.axvline(end_signal_pos, color='red', linestyle='-', label='Signal Threshold')
        plt.axvline(end_signal_pos - params_window_radius, color='Blue', linestyle='-', label='Signal + Windows')
        plt.axvline(end_signal_pos - 2 * params_window_radius, color='Skyblue', linestyle='-', label='Signal + (2 x Windows)')
    elif _nonempty(mean_end_sigs):
        plt.axvline(mean_end_signal_pos, color='brown', linestyle='-', label='Mean Signal Threshold')
        plt.axvline(mean_end_signal_pos - params_window_radius, color='darkgoldenrod', linestyle='-', label='Mean Signal + Windows')
        plt.axvline(mean_end_signal_pos - 2 * params_window_radius, color='gold', linestyle='-', label='Mean Signal + (2 x Windows)')

    # legend (unique)
    handles, labels = plt.gca().get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    plt.legend(uniq.values(), uniq.keys(), loc='center left', bbox_to_anchor=(1, 0.75))

    plt.grid(True)
    if n is not None:
        plt.xlim(left=0, right=n)

    if save:
        visuals_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()

    return out_path if save else None



from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Union

def plot_consensus_frequency(
    *,
    df_clean,                       # DataFrame with 'freq' column OR a 1D array/Series of frequencies
    block_num: int,
    start_seed_pos: int,
    end_seed_pos: int,
    # signal windows
    start_sigs=None,
    end_sigs=None,
    start_signal_pos: Optional[int] = None,
    end_signal_pos: Optional[int] = None,
    params_window_radius: int = 0,
    # block detection flags & positions
    signal_start: bool = False,
    mean_signal_start: bool = False,
    signal_end: bool = False,
    mean_signal_end: bool = False,
    start_block_pos: Optional[int] = None,
    end_block_pos: Optional[int] = None,
    # layout / saving
    n: Optional[int] = None,
    visuals_dir: Union[str, Path, None] = None,
    filename: Optional[str] = None,
    save: bool = True,
    show: bool = True,
):
    """
    Plot consensus alignment base frequency for a block.
    Saves the figure (SVG) if `save=True` and returns the Path; otherwise returns None.
    """

    def _nonempty(x):
        if x is None:
            return False
        if hasattr(x, "empty"):
            return not x.empty  # pandas objects
        try:
            return len(x) > 0
        except Exception:
            return True

    # resolve y-values from input
    try:
        y = df_clean["freq"]
    except Exception:
        y = df_clean  # assume it's already a 1D sequence

    visuals_dir = Path(visuals_dir) if visuals_dir is not None else Path(".")
    if filename is None:
        filename = f"consensus_msa_frequency_block_{block_num}.svg"
    out_path = visuals_dir / filename

    plt.figure(figsize=(12, 5))
    plt.yticks(np.linspace(0, 1, 21).round(2))
    plt.xlabel("Alignment Position")
    plt.ylabel("Frequency")
    plt.title(f"Consensus Alignment Base Frequency Block: {block_num}")

    plt.plot(y, label="Gap-cleaned Alignment", linewidth=2)

    # start side thresholds
    if _nonempty(start_sigs):
        plt.axvline(start_signal_pos, color='red', linestyle='--', label='Signal Threshold')
        plt.axvline(start_signal_pos +  params_window_radius, color='Blue', linestyle='--', label='Signal + Windows')
        plt.axvline(start_signal_pos + 2 * params_window_radius, color='Skyblue', linestyle='--', label='Signal + (2 x Windows)')

    # end side thresholds
    if _nonempty(end_sigs):
        plt.axvline(end_signal_pos, color='red', linestyle='--', label='Signal Threshold')
        plt.axvline(end_signal_pos - params_window_radius, color='Blue', linestyle='--', label='Signal + Windows')
        plt.axvline(end_signal_pos - 2 * params_window_radius, color='Skyblue', linestyle='--', label='Signal + (2 x Windows)')

    # seed boundaries
    plt.axvline(start_seed_pos, color='gray', linestyle='-', linewidth=2, label='Seed Start')
    plt.axvline(end_seed_pos,   color='gray', linestyle='-', linewidth=2, label='Seed End')

    # block boundaries when both ends are confirmed
    if (signal_start or mean_signal_start) and (signal_end or mean_signal_end):
        plt.axvline(start_block_pos,   color="Green",   linestyle='-', label="Block Start")
        plt.axvline((end_block_pos or 0) + 1, color="Fuchsia", linestyle='-', label="Block End")

    # legend (unique)
    handles, labels = plt.gca().get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    plt.legend(unique.values(), unique.keys(), loc='center left', bbox_to_anchor=(1, 0.75))

    plt.grid(visible=True)
    if n is not None:
        plt.xlim(left=0, right=n)
    plt.ylim(bottom=0, top=1.05)

    if save:
        visuals_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, bbox_inches="tight")
    if show:
        plt.show()
    plt.close()

    return out_path if save else None


from pathlib import Path
import re
from typing import Callable, Dict, Iterable, List, Optional, Union

import numpy as np
import pandas as pd


def extend_seeds_to_blocks(
    *,
    output_dir: Union[str, Path],
    seeds_df: pd.DataFrame,
    seeds_hits_df: pd.DataFrame,
    blast_df: pd.DataFrame,
    input_reads_fasta: Union[str, Path],
    # external helpers you already have:
    align_seed_neighbours: Callable,
    find_block_boundaries: Callable,
    blastn_subject: Callable,
    default_blast_columns: Iterable[str],
    # params (kept 1:1 with your variables)
    params_aln_w_step_size: int,
    params_window_radius: int,
    params_signal_threshold: float,
    params_mean_signal_threshold: float,
    params_max_aln_extension: int,
    params_seed_map_coverage_threshold: float,
    # tuning knobs with sensible defaults
    boundary_k: int = 10,
    boundary_threshold: float = 0.10,
    min_aln_len: int = 20,
    min_pident: float = 85.0,
    word_size: int = 10,
    gap_extend_blast: int = 2,
    gap_open_blast: int = 5,
    penalty_blast: int = -2,
    threads: int = 10,
    max_iter: int = 100,
    gap_open: float = 1.0,
    gap_extend: float = 0.1,
    show_cmd: bool = False,
    make_plots: bool = True,
    log: Callable[[str], None] = print,
):
    """
    Extends each seed to detect a block, writes outputs, optionally plots, and maps seeds covered by the block.

    Returns:
        dict with:
          - 'blocks_fasta' (Path): path to the aggregated blocks FASTA
          - 'mapped_seeds' (set[str])
          - 'visuals_dir' (Path)
          - 'blocks' (List[Dict]): per-block summary (seed, positions, lengths, sequence, etc.)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---------- Prepare FASTA with all seeds (subject for later BLAST) ----------
    temp_seeds_fasta = output_dir / "temp.seeds.fasta"
    with open(temp_seeds_fasta, "w") as fasta:
        for _, row in seeds_df.iterrows():
            fasta.write(f">{row['qseqid']}\n{row['seq']}\n")

    mapped_seeds: set = set()

    output_blocks_fasta = output_dir / "blocks.fasta"
    temp_block_fasta = output_dir / "block-temp.fasta"
    visuals_dir = output_dir / "block_extension_plots"
    visuals_dir.mkdir(parents=True, exist_ok=True)

    # init final FASTA
    with open(output_blocks_fasta, "w") as fasta:
        fasta.write("")

    blocks_summary: List[Dict] = []
    block_num = 0

    # ------------------------------ main loop ------------------------------
    for _, hit_row in seeds_hits_df.iterrows():
        seed = hit_row["qseqid"]
        log(f"Working on target seed: {seed}")

        if seed in mapped_seeds:
            continue
        else:
            log(f"Mapping seed: {seed}")

        block_num += 1
        log(f"\nProcessing Block Number: {block_num}\n")

        # reads mapped to this seed
        seed_mapped_reads = blast_df.loc[blast_df["qseqid"] == seed, ["sseqid", "sstart", "send", "slen"]]
        seed_mapped_id = set(seed_mapped_reads["sseqid"])

        # flags
        mean_signal_start = False
        mean_signal_end = False
        signal_start = False
        signal_end = False

        # seed size
        try:
            seed_size = int(blast_df.loc[blast_df["qseqid"] == seed, "qlen"].iloc[0])
        except Exception:
            log(f"Warning: cannot determine seed size for {seed}. Skipping.")
            continue
        log(f"Seed initial size: {seed_size}")

        # state variables
        s_aln_extension = 0
        e_aln_extension = 0
        block_seq: Optional[str] = None

        # placeholders for plotting / summary (defined before loop to avoid UnboundLocal)
        start_seed_pos = None
        end_seed_pos = None
        start_sigs = pd.Series(dtype=int)
        end_sigs = pd.Series(dtype=int)
        mean_start_sigs = pd.Series(dtype=int)
        mean_end_sigs = pd.Series(dtype=int)
        start_signal_pos = None
        end_signal_pos = None
        mean_start_signal_pos = None
        mean_end_signal_pos = None
        start_block_pos = None
        end_block_pos = None
        n = None
        x_vals = np.array([], dtype=int)
        plot_means = np.array([], dtype=float)
        plot_stds = np.array([], dtype=float)
        plot_signal = np.array([], dtype=float)
        df_clean = None

        # ------------------------------ grow windows until we can define a block or give up ------------------------------
        while not block_seq:
            if not (signal_start or mean_signal_start):
                s_aln_extension += params_aln_w_step_size
            if not (signal_end or mean_signal_end):
                e_aln_extension += params_aln_w_step_size

            # ---- MSA ----
            aligned, freq_table, consensus_df, real_s_extension, real_e_extension = align_seed_neighbours(
                reads_fasta=input_reads_fasta,
                seed_hits=seed_mapped_reads,
                s_extension=int(s_aln_extension),
                e_extension=int(e_aln_extension),
                threads=threads,
                max_iter=max_iter,
                gap_open=gap_open,
                gap_extend=gap_extend,
                show_cmd=show_cmd,
            )

            # remove gaps from consensus trace
            df_clean = consensus_df[consensus_df["base"] != "-"].copy()
            df_clean = df_clean.reset_index(names="aln_pos")

            n = len(df_clean)
            if n < 2 * params_window_radius:
                log("Not enough data points to compute windowed metrics. Extending more...")
                # continue to extend more
                if (s_aln_extension + e_aln_extension + seed_size) > params_max_aln_extension:
                    log(f"Max alignment extension reached for {seed}. Aborting this seed.")
                    break
                continue

            # ---- windowed signals ----
            w_means = np.full(n, np.nan)
            w_stds = np.full(n, np.nan)
            w_signal = np.full(n, np.nan)

            freqs = df_clean["freq"].values
            signal_pos_list: List[int] = []
            for pos in range(params_window_radius, n - params_window_radius):
                segment = freqs[pos - params_window_radius : pos + params_window_radius]
                mean_val = segment.mean()
                std_val = segment.std()
                sig_val = std_val * (2 - mean_val)

                w_means[pos] = mean_val
                w_stds[pos] = std_val
                w_signal[pos] = sig_val

                if sig_val > params_signal_threshold:
                    signal_pos_list.append(pos)

            valid_mask = ~np.isnan(w_signal)
            x_vals = np.where(valid_mask)[0]
            plot_means = w_means[valid_mask]
            plot_stds = w_stds[valid_mask]
            plot_signal = w_signal[valid_mask]

            # moving mean of the signal
            w_mean_sig = np.full(len(plot_signal), np.nan)
            mean_signal_pos_list: List[int] = []
            for pos in range(params_window_radius, len(plot_signal) - params_window_radius):
                segment = plot_signal[pos - params_window_radius : pos + params_window_radius]
                mean_sig = segment.mean()
                w_mean_sig[pos] = mean_sig
                if mean_sig > params_mean_signal_threshold:
                    mean_signal_pos_list.append(pos + params_window_radius)

            if plot_signal.size == 0:
                log("Warning: no valid signal values. Skipping plot/seed.")
                break

            # ---- seed position in consensus ----
            total_aln_extension = s_aln_extension + seed_size + e_aln_extension
            if total_aln_extension != n:
                log(
                    f"Warning! Alignment extension ({n}) != expected ({total_aln_extension}). "
                    "Attempting to locate seed in consensus."
                )
                try:
                    seed_seq = str(seeds_df.loc[seeds_df["qseqid"] == seed, "seq"].values[0]).upper()
                except Exception:
                    log("Could not retrieve seed sequence. Skipping seed.")
                    break
                consensus_seq = "".join(df_clean["base"].astype(str)).upper()
                matches = list(re.finditer(seed_seq, consensus_seq))

                if len(matches) > 1:
                    log("Multiple matches detected for the same seed. Skipping seed.")
                    break
                elif len(matches) == 0:
                    log("No exact matches for seed in consensus. Skipping seed.")
                    break
                else:
                    start_seed_pos = matches[0].start()
                    end_seed_pos = matches[0].end()

                if total_aln_extension > params_max_aln_extension:
                    log(
                        f"Warning! Skipping seed {seed}, exceeded max alignment extension "
                        f"({params_max_aln_extension})."
                    )
                    break
            else:
                start_seed_pos = s_aln_extension
                end_seed_pos = s_aln_extension + seed_size

            # ---- find signals relative to the seed window ----
            signal_pos = pd.Series(signal_pos_list, dtype=int)
            mean_signal_pos = pd.Series(mean_signal_pos_list, dtype=int)

            start_sigs = signal_pos[signal_pos < start_seed_pos]
            end_sigs = signal_pos[signal_pos > end_seed_pos]

            mean_start_sigs = mean_signal_pos[mean_signal_pos < start_seed_pos]
            mean_end_sigs = mean_signal_pos[mean_signal_pos > end_seed_pos]

            # positions of max mean (for plotting)
            try:
                max_indices = x_vals[plot_means == plot_means.max()]
            except ValueError:
                max_indices = np.array([], dtype=int)

            if not mean_start_sigs.empty:
                mean_start_signal_pos = int(max(mean_start_sigs))
                mean_signal_start = True
            if not mean_end_sigs.empty:
                mean_end_signal_pos = int(min(mean_end_sigs))
                mean_signal_end = True
            if not start_sigs.empty:
                start_signal_pos = int(max(start_sigs))
                signal_start = True
            if not end_sigs.empty:
                end_signal_pos = int(min(end_sigs))
                signal_end = True

            # ---- define block if both ends are present ----
            if (signal_start or mean_signal_start) and (signal_end or mean_signal_end):
                s_sig_pos = start_signal_pos if signal_start else (mean_start_signal_pos + params_window_radius)
                e_sig_pos = end_signal_pos if signal_end else (mean_end_signal_pos - params_window_radius)

                log(f"\nStart signal Pos: {start_signal_pos}")
                log(f"Mean Start signal Pos: {mean_start_signal_pos}")
                log(f"End signal Pos: {end_signal_pos}")
                log(f"Mean End signal Pos: {mean_end_signal_pos}\n")

                start_block_pos, end_block_pos = find_block_boundaries(
                    df_clean,
                    s_sig_pos,
                    e_sig_pos,
                    start_seed_pos,
                    end_seed_pos,
                    boundary_k,
                    threshold=boundary_threshold,
                )

                if start_block_pos is None or end_block_pos is None:
                    log("Warning: find_block_boundaries couldn't find a result. Skipping seed.")
                    log(f"Start Block Position: {start_block_pos}")
                    log(f"End Block Position: {end_block_pos}")
                    break

                # extract sequence
                block_seq = "".join(df_clean["base"][start_block_pos : end_block_pos + 1])

                # write final and temp block FASTAs
                with open(output_blocks_fasta, "a") as fasta, open(temp_block_fasta, "w") as temp_fasta:
                    fasta.write(f">block_{block_num}\n{block_seq}\n")
                    temp_fasta.write(f">block_{block_num}\n{block_seq}\n")

                # mark seeds included in this block (by BLAST mapping of seeds -> block)
                mapped_seeds.add(seed)
                
                #### - - - - Add to mapped seeds those that are present in the block - - - - - - - - - - - - - 
                from pandas.errors import EmptyDataError
                
                blast_output_table_tsv = output_dir / "temp.blast.seeds_fragment.tsv"
                blastn_subject(
                    blast_input_seqs=temp_seeds_fasta,
                    blast_subject_seqs=temp_block_fasta,
                    blast_output_table_tsv=blast_output_table_tsv,
                    word_size=word_size,
                    gap_extend=gap_extend_blast,
                    gapopen=gap_open_blast,
                    penalty=penalty_blast,
                )
                
                # Try to read BLAST output
                try:
                    seeds_block_df = pd.read_csv(blast_output_table_tsv, sep='\t', header=None)
                except EmptyDataError:
                    # File exists but is empty → skip this iteration
                    print(f"Skipping {blast_output_table_tsv}, file is empty.")
                    continue  

                seeds_block_df = pd.read_csv(blast_output_table_tsv, sep="\t")
                seeds_block_df.columns = list(default_blast_columns)

                seed_coverage_filter = (seeds_block_df["length"] / seeds_block_df["qlen"]) > params_seed_map_coverage_threshold
                aln_len_filter = seeds_block_df["length"] >= min_aln_len
                aln_pident_filter = seeds_block_df["pident"] > min_pident
                seeds_block_df = seeds_block_df[aln_len_filter & aln_pident_filter & seed_coverage_filter]

                mapped_seeds.update(seeds_block_df["qseqid"])

                #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                
                # block summary
                blocks_summary.append(
                    dict(
                        block_num=block_num,
                        seed=seed,
                        seed_size=seed_size,
                        s_aln_extension=s_aln_extension,
                        e_aln_extension=e_aln_extension,
                        start_seed_pos=int(start_seed_pos),
                        end_seed_pos=int(end_seed_pos),
                        start_block_pos=int(start_block_pos),
                        end_block_pos=int(end_block_pos),
                        block_len=int(end_block_pos - start_block_pos + 1),
                        s_sig_pos=int(s_sig_pos),
                        e_sig_pos=int(e_sig_pos),
                        s_sig_source="signal" if signal_start else "mean",
                        e_sig_source="signal" if signal_end else "mean",
                        block_seq=block_seq,
                        n=n,
                    )
                )

            # ---- plotting (if we got a block OR both ends detected) ----
            if make_plots and (
                (block_seq is not None) or ((signal_start or mean_signal_start) and (signal_end or mean_signal_end))
            ):
                try:
                    plot_windowed_metrics(
                        x_vals=x_vals,
                        plot_means=plot_means,
                        plot_stds=plot_stds,
                        plot_signal=plot_signal,
                        block_num=block_num,
                        start_seed_pos=start_seed_pos,
                        end_seed_pos=end_seed_pos,
                        max_indices=(x_vals[plot_means == plot_means.max()] if plot_means.size else []),
                        start_sigs=start_sigs,
                        mean_start_sigs=mean_start_sigs,
                        end_sigs=end_sigs,
                        mean_end_sigs=mean_end_sigs,
                        start_signal_pos=start_signal_pos,
                        mean_start_signal_pos=mean_start_signal_pos,
                        end_signal_pos=end_signal_pos,
                        mean_end_signal_pos=mean_end_signal_pos,
                        params_window_radius=params_window_radius,
                        n=n,
                        visuals_dir=visuals_dir,
                    )

                    plot_consensus_frequency(
                        df_clean=df_clean,
                        block_num=block_num,
                        start_seed_pos=start_seed_pos,
                        end_seed_pos=end_seed_pos,
                        start_sigs=start_sigs,
                        end_sigs=end_sigs,
                        start_signal_pos=start_signal_pos,
                        end_signal_pos=end_signal_pos,
                        params_window_radius=params_window_radius,
                        signal_start=signal_start,
                        mean_signal_start=mean_signal_start,
                        signal_end=signal_end,
                        mean_signal_end=mean_signal_end,
                        start_block_pos=start_block_pos,
                        end_block_pos=end_block_pos,
                        n=n,
                        visuals_dir=visuals_dir,
                    )
                except NameError:
                    # If the user didn't import the plot helpers, just skip plotting silently
                    log("Plot helpers not found; skipping plots for this block.")

            # if we found a block, end the while loop; else continue extending
            if block_seq is not None:
                break

            # safety: stop if we already exceeded maximum extension
            if (s_aln_extension + e_aln_extension + seed_size) > params_max_aln_extension:
                log(f"Reached max extension for seed {seed} without forming a block. Skipping.")
                break

    return dict(
        blocks_fasta=output_blocks_fasta,
        mapped_seeds=mapped_seeds,
        visuals_dir=visuals_dir,
        blocks=blocks_summary,
    )

import pandas as pd
from Bio import SeqIO
from utils.blast_utils import makeblast_db, blastn, default_blast_columns


def run_blast_mapping(
    input_clean_fasta,
    input_auto_blocks_fasta,
    output_map_data_dir,
    logs_dir,
    word_size=20,
    num_threads=10,
    pident_threshold=0.9,
    gapopen_threshold=20,
    query_coverage_threshold=0.9
) -> float:
    """
    Run BLAST mapping of auto blocks against a clean fasta database and compute mapping statistics.

    Parameters
    ----------
    input_clean_fasta : str or Path
        Path to the clean fasta file.
    input_auto_blocks_fasta : str or Path
        Path to the fasta file with auto blocks.
    output_map_data_dir : Path
        Directory where BLAST outputs will be written.
    logs_dir : Path
        Directory where log files will be written.
    word_size : int, optional
        Word size for BLASTN, default = 20.
    num_threads : int, optional
        Number of threads for BLASTN, default = 10.
    pident_threshold : float, optional
        Minimum identity threshold to keep alignments, default = 0.9.
    gapopen_threshold : int, optional
        Maximum allowed gapopen value, default = 20.

    Returns
    -------
    map_proportion: float 
    """

    # Build BLAST database
    clean_db = makeblast_db(
        seqs_path=input_clean_fasta,
        db_out=output_map_data_dir / "clean_fasta_db",
        log_file=logs_dir / "make_db.log"
    )

    # Run BLAST
    blast_tsv = output_map_data_dir / "blastn.auto_bl_clean_db.tsv"
    blastn(
        blast_input_seqs=input_auto_blocks_fasta,
        blast_db_file=clean_db,
        word_size=word_size,
        num_threads=num_threads,
        blast_output_table_tsv=blast_tsv,
        log_file=logs_dir / "blastn.auto_bl.log"
    )

    # Load BLAST results
    blast_df = pd.read_csv(blast_tsv, sep='\t', header=None)
    blast_df.columns = default_blast_columns

    mapped_subjects = set(blast_df['sseqid'])

    # Filter BLAST results
    blast_df = blast_df[blast_df['pident'] > pident_threshold]
    blast_df = blast_df[blast_df['gapopen'] < gapopen_threshold]
    blast_df = blast_df[(blast_df['length'] / blast_df['qlen']) > query_coverage_threshold]

    unmapped_id, unmapped_fasta = [], []
    for record in SeqIO.parse(input_clean_fasta, 'fasta'):
        if record.id not in mapped_subjects:
            unmapped_id.append(record.id)
            unmapped_fasta.append(record.seq)

    unmapped_df = pd.DataFrame({'id': unmapped_id, 'fasta': unmapped_fasta})

    if unmapped_df.empty:
        unmapped_bases = 0
    else:
        unmapped_bases = unmapped_df['fasta'].apply(len).sum()

    total_bases = blast_df.groupby('sseqid')['slen'].first().sum() + unmapped_bases
    mapped_bases = blast_df['length'].sum()

    mapped_prop = mapped_bases / total_bases if total_bases > 0 else 0

    return mapped_prop
