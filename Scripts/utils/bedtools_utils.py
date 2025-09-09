import subprocess
import pandas as pd
from pathlib import Path
import subprocess
import shutil
import sys

def bedtools_getfasta(
    fasta_path: str,
    bed_path: str,
    extracted_regions_path: str,
    show_command: bool = False,
    log_file: str = None
) -> None:
    """
    Extracts regions from `fasta_path` as defined in `bed_path` using bedtools,
    writing to `extracted_regions_path`.

    If `show_command` is True, prints the command before running.
    If `log_file` is provided, appends stdout/stderr to that file.
    """
    cmd = [
        "bedtools", "getfasta",
        "-fi", fasta_path,
        "-bed", bed_path,
        "-fo", extracted_regions_path
    ]

    if show_command:
        print("Running:", " ".join(map(str, cmd)))

    try:
        if log_file:
            with open(log_file, "a") as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        else:
            subprocess.run(cmd, check=True)
        #print("bedtools ran successfully.")
    except subprocess.CalledProcessError as e:
        print("Error running bedtools:", e)
        raise

            
def coords_to_regions(df: pd.DataFrame, id_col: str = 'qseqid', coords_col: str = 'coords') -> pd.DataFrame:
    """
    Generates a DataFrame suitable for FASTA extraction using start and end coordinates.

    Parameters:
        df (pd.DataFrame): A DataFrame with at least two columns: an ID column and a coordinate column.
        id_col (str): Name of the column containing sequence IDs.
        coords_col (str): Name of the column containing coordinates (sorted integers).

    Returns:
        pd.DataFrame: A DataFrame with 'id', 'start', and 'end' columns.
    """
    sorted_coords = sorted(df[coords_col])
    
    # Create pairs of consecutive coordinates: (start, end)
    starts = sorted_coords[:-1]
    ends = sorted_coords[1:]

    # Repeat the same ID for all coordinate pairs
    fasta_df = pd.DataFrame({
        'id': [df[id_col].iloc[0]] * len(starts),
        'start': starts,
        'end': ends
    })
    
    # Ensure start and end columns are integers (in case of type coercion)
    fasta_df['start'] = fasta_df['start'].astype(int)
    fasta_df['end'] = fasta_df['end'].astype(int)
    
    return fasta_df

def blast_to_bed(blast_df: pd.DataFrame, seqid: str = 'qseqid', start_coord: str = 'qstart', end_coord: str = 'qend') -> pd.DataFrame:
    bed_df = blast_df[[seqid, start_coord, end_coord]]
    return(bed_df)


def bbtools_kmercountexact(
    infile: Path,
    outfile: Path,
    k: int = 31,
    min_counts: int = 30,
    khist: Path | None = None,
    make_fasta: bool = False,
    threads: int = 1,
    extra: list[str] | None = None,
    exe: str = "kmercountexact.sh",
    overwrite = True,
    log_file: str | None = None 
) -> None:
    """
    Build and execute the kmercountexact.sh command.

    Parameters
    ----------
    infile : Path
        Input FASTA/FASTQ file (gzipped ok).
    outfile : Path
        Tab-separated file of k-mer counts.
    k : int, default 31
        K-mer length.
    min_counts : int, default 30
        Minimum count threshold for reporting.
    khist : Path, optional
        File to write the k-mer frequency histogram.
    make_fasta : bool, default False
        Whether to dump k-mers in FASTA format.
    threads : int, default 1
        Number of CPU threads.
    extra : list of str, optional
        Additional BBTools parameters (e.g. ["mem=32g"]).
    exe : str, default "kmercountexact.sh"
        Path to the BBTools launcher.
    """
    if shutil.which(exe) is None:
        sys.exit(f"❌ {exe} not found on $PATH – install BBMap/BBTools first.")

    infile = Path(infile).expanduser().resolve()
    if not infile.exists():
        sys.exit(f"❌ Input file {infile} does not exist.")

    fastadump = 't' if make_fasta else 'f'

    cmd = [
        exe,
        f"in={infile}",
        f"k={k}",
        f"threads={threads}",
        f"mincount={min_counts}",
        f"fastadump={fastadump}",
        f"out={outfile}"
    ]
    if khist:
        cmd.append(f"khist={khist}")
    if overwrite:
        cmd.append(f"overwrite=true")
    if extra:
        cmd.extend(extra)

    print("▶ Running:", " ".join(cmd))
    
    try:
        if log_file:
            with open(log_file, "a") as log:
                subprocess.run(cmd, check=True, stdout=log, stderr=log)
        else:
            subprocess.run(cmd, check=True)
        
    except subprocess.CalledProcessError as e:
        sys.exit(f"❌ kmercountexact.sh failed (exit code {e.returncode}).")
        raise
