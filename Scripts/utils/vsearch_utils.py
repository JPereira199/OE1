import subprocess

def vsearch_sort_by_length(fasta_file: str, sorted_output_file: str, 
                           log_file: str = None, reverse: bool = False, threads: int = 40) -> None:
    """
    Sorts sequences by length using VSEARCH.

    Args:
        fasta_file (str): Path to the input FASTA file with sequences to be sorted.
        sorted_output_file (str): Path to output the sorted sequences.
        log_file (str): Path to the log file. If None, logs are printed to the console.
        reverse (bool): If True, sorts sequences from shortest to longest. Default is False (longest to shortest).
    """
    # Construct the vsearch command
    cmd = [
        'vsearch', '--sortbylength', fasta_file,
        '--output', sorted_output_file,
        '--threads', str(threads)
    ]
    
    # Reverse sort if specified
    if reverse:
        cmd.append('--reversesort')
    
    print("Running command:", " ".join(map(str, cmd)))
    
    # Run the command and log output
    log = open(log_file, "a") if log_file else None
    try:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)
        print("vsearch ran successfully.")
    except subprocess.CalledProcessError as e:
        print("Error in running vsearch:", e)
    finally:
        if log:
            log.close()


def vsearch_cluster_unoise(fasta_file: str, fasta_denoised_file: str, maxseqlength: int = 100000, 
                           minseqlength: int = 25, minsize: int = 1, unoise_alpha: float = 2.0, strand_both: bool = True,
                           sizein: bool = False, sizeout: bool = False, log_file: str = None, threads: int = 40) -> None:
    cmd = [
        'vsearch', '--cluster_unoise', fasta_file, '--centroids', fasta_denoised_file,
        '--maxseqlength', str(maxseqlength), '--minseqlength', str(minseqlength),
        '--minsize', str(minsize), '--unoise_alpha', str(unoise_alpha), '--threads', str(threads)
    ]
    
    if sizein:
        cmd.append('--sizein')
    if sizeout:
        cmd.append('--sizeout')
    if strand_both:
        cmd.extend(['--strand', 'both'])

    print("Running command:", " ".join(map(str, cmd)))

    log = open(log_file, "a") if log_file else None
    try:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)
        print("vsearch ran successfully.")
    except subprocess.CalledProcessError as e:
        print("Error in running vsearch:", e)
    finally:
        if log:
            log.close()


def vsearch_cluster_fast(fasta_file: str, fasta_clustered_file: str, identity: float = 0.99, 
                         maxseqlength: int = 100000, minseqlength: int = 25, strand_both: bool = True,
                         sizein: bool = False, sizeout: bool = False, log_file: str = None, threads: int = 40) -> None:
    cmd = [
        'vsearch', '--cluster_fast', fasta_file,
        '--centroids', fasta_clustered_file,
        '--id', str(identity),
        '--maxseqlength', str(maxseqlength),
        '--minseqlength', str(minseqlength),
        '--threads', str(threads)
    ]
    
    if sizein:
        cmd.append('--sizein')
    if sizeout:
        cmd.append('--sizeout')
    if strand_both:
        cmd.extend(['--strand', 'both'])
    
    print("Running command:", " ".join(map(str, cmd)))
    
    log = open(log_file, "a") if log_file else None
    try:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)
        print("vsearch ran successfully.")
    except subprocess.CalledProcessError as e:
        print("Error in running vsearch:", e)
    finally:
        if log:
            log.close()

def vsearch_dereplication(input_fasta: str, derep_fasta: str, min_seq_length: int = 32, log_file: str = None):
    """Run VSEARCH dereplication."""
    cmd = [
        "vsearch", "--derep_fulllength", input_fasta,
        "--output", derep_fasta,
        "--sizeout", "--relabel", "Seq",
        "--minseqlength", str(min_seq_length),
    ]
    print(f"Executing command: {' '.join(map(str, cmd))}")

    with open(log_file, "a") if log_file else subprocess.DEVNULL as log:
        try:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
        except subprocess.CalledProcessError as e:
            print("❌ Dereplication failed.")
            raise

def vsearch_sortbysize(derep_fasta: str, sorted_fasta: str, log_file: str = None):
    """Run VSEARCH sortbysize."""
    cmd = [
        "vsearch", "--sortbysize", derep_fasta,
        "--output", sorted_fasta
    ]
    print(f"Executing command: {' '.join(map(str, cmd))}")

    with open(log_file, "a") if log_file else subprocess.DEVNULL as log:
        try:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
        except subprocess.CalledProcessError as e:
            print("❌ Sorting by size failed.")
            raise

def vsearch_clustersize(sorted_fasta: str, centroids_fasta: str, uc_file: str,
                        id_thresh: float = 0.96, query_cov: float = 0.96, target_cov: float = 0.96,
                        min_seq_length: int = 32, log_file: str = None,
                        use_both_strands: bool = False, use_sizein: bool = False):
    """
    Run VSEARCH clustering by size with coverage filters, optional strand awareness,
    and optional use of abundance information from dereplication (--sizein).
    """
    cmd = [
        "vsearch", "--cluster_size", sorted_fasta,
        "--id", str(id_thresh),
        "--query_cov", str(query_cov),
        "--target_cov", str(target_cov),
        "--minseqlength", str(min_seq_length),
        "--centroids", centroids_fasta,
        "--uc", uc_file
    ]

    if use_both_strands:
        cmd.extend(["--strand", "both"])

    if use_sizein:
        cmd.append("--sizein")

    print(f"Executing command: {' '.join(map(str, cmd))}")

    with open(log_file, "a") if log_file else subprocess.DEVNULL as log:
        try:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
        except subprocess.CalledProcessError as e:
            print("❌ Clustering failed.")
            raise


def vsearch_pipeline_clustersize(input_fasta: str, derep_fasta: str,
                         sorted_fasta: str, centroids_fasta: str, uc_file: str):
    """Main pipeline: dereplicate → sort → cluster."""
    vsearch_dereplication(input_fasta, derep_fasta)
    vsearch_sortbysize(derep_fasta, sorted_fasta)
    vsearch_clustersize(sorted_fasta, centroids_fasta, uc_file)