############################################ Loading libraries and user arguments ###########################################

import os
from utils.blast_utils import unnoise_coords, makeblast_db, blastn, alignment_absolute_start_end, select_internal_aligns, retrieve_single_alignments
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from Bio import SeqIO
from pathlib import Path  

def str2bool(v: str) -> bool:
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ("yes", "y", "true", "t", "1"):
        return True
    if v in ("no", "n", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("boolean value expected")

set_trim_nam = {
    'input_fasta' : '/home/jpereira/OEs/Blocks_Namasivayam.fa',
    'param_min_instances' : 0, # Must be always 0 when you are working with extended seeds
    'param_threads' : 20,
    'params_cluster_identity' : 0.85,
    'params_cluster_coverage' : 0.85,
    'params_select_internals_alns' : False, # Must be always false when you are working with extended seeds
    'params_retrive_sinlge_aln' : False,
    'output_dir' : Path('/home/jpereira/OEs/Results/OE1/NamSeqs/Data/define_blocks/')
    }

set_auto_bl_nam = {
    'input_fasta' : '/home/jpereira/OEs/Results/OE1/NamSeqs/Data/seed_extension/kmer30/blocks.fasta',
    'param_min_instances' : 0, # Must be always 0 when you are working with extended seeds
    'param_threads' : 20,
    'params_cluster_identity' : 0.85,
    'params_cluster_coverage' : 0.85,
    'params_select_internals_alns' : False, # Must be always false when you are working with extended seeds
    'params_retrive_sinlge_aln' : False,
    'output_dir' : Path('/home/jpereira/OEs/Results/OE1/NamSeqs/Data/define_blocks/')
    }

set_auto_bl_toxop1 = {
    'input_fasta' : '/home/jpereira/OEs/Results/OE1/ToxoPasteur1/Data/seed_extension/toxo_pasteur1/blocks.fasta',
    'param_min_instances' : 0, # Must be always 0 when you are working with extended seeds
    'param_threads' : 20,
    'params_cluster_identity' : 0.85,
    'params_cluster_coverage' : 0.85,
    'params_select_internals_alns' : False, # Must be always false when you are working with extended seeds
    'params_retrive_sinlge_aln' : False,
    'output_dir' : Path('/home/jpereira/OEs/Results/OE1/ToxoPasteur1/Data/define_blocks/')
    }

set_auto_bl_toxop2 = {
    'input_fasta' : '/home/jpereira/OEs/Results/OE1/ToxoPasteur1/Data/seed_extension/toxo_pasteur2/blocks.fasta',
    'param_min_instances' : 0, # Must be always 0 when you are working with extended seeds
    'param_threads' : 20,
    'params_cluster_identity' : 0.85,
    'params_cluster_coverage' : 0.85,
    'params_select_internals_alns' : False, # Must be always false when you are working with extended seeds
    'params_retrive_sinlge_aln' : False,
    'output_dir' : Path('/home/jpereira/OEs/Results/OE1/ToxoPasteur1/Data/define_blocks/')
    }

parser = argparse.ArgumentParser()
parser.add_argument('--input-fasta', required=True, type=Path)
parser.add_argument('--param-min-instances',required=True, type=int )
parser.add_argument('--param-threads', type=int, default=10 )
parser.add_argument('--params-cluster-identity',type=float, default=0.85 )
parser.add_argument('--params-cluster-coverage', type=float, default=0.85 )
parser.add_argument('--params-select-internals-alns',  choices=('False', 'false', 'True', 'true' ) )
parser.add_argument('--params-retrive-sinlge-aln', choices=('False', 'false', 'True', 'true' ) )
parser.add_argument('--output-dir', type=Path)

args = parser.parse_args()


set_input_args = {
    'input_fasta' : args.input_fasta,
    'param_min_instances' : args.param_min_instances,
    'param_threads' : args.param_threads,
    'params_cluster_identity' : args.params_cluster_identity,
    'params_cluster_coverage' : args.params_cluster_coverage,
    'params_select_internals_alns' : str2bool(args.params_select_internals_alns),
    'params_retrive_sinlge_aln' : str2bool(args.params_retrive_sinlge_aln),
    'output_dir' : args.output_dir,
    
}

settings_dict = {
    'trim_nam' : set_trim_nam,
    'auto_bl_nam' : set_auto_bl_nam,
    'auto_bl_toxop1': set_auto_bl_toxop1, 
    'auto_bl_toxop2': set_auto_bl_toxop2, 
    None: set_input_args
}

work_name = None

input_fasta = settings_dict[work_name]['input_fasta']
param_min_instances = settings_dict[work_name]['param_min_instances']
param_threads = settings_dict[work_name]['param_threads']
params_cluster_identity = settings_dict[work_name]['params_cluster_identity']
params_cluster_coverage = settings_dict[work_name]['params_cluster_coverage']
params_select_internals_alns = settings_dict[work_name]['params_select_internals_alns']
params_retrive_sinlge_aln = settings_dict[work_name]['params_retrive_sinlge_aln']
output_dir = settings_dict[work_name]['output_dir']


if work_name:
    output_dir = output_dir / work_name    

output_dir.mkdir(exist_ok=True, parents=True)

########################################## Defining Helper Function ##########################################

from Bio import SeqIO

def get_fasta_from_bed(fasta_file, bed_file, output_file):
    # Read the FASTA file into a dictionary for quick access
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    with open(bed_file, 'r') as bed, open(output_file, 'w') as out:
        for line in bed:
            if line.strip() == "" or line.startswith("#"):
                continue  # Skip empty lines and comments

            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
            
            if chrom not in fasta_dict:
                raise ValueError(f"Sequence {chrom} not found in FASTA file.")
            
            seq = fasta_dict[chrom].seq[start-1:end]
            out.write(f">{chrom}:{start}-{end}\n{seq}\n")


##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

output_reciprocal_blast_tsv = output_dir / "blastn.infasta_reciprocal.tsv"
infasta_db_dir = os.path.join(os.path.dirname(output_reciprocal_blast_tsv), 'blast_db', 'infasta_db') #'/home/jpereira/OEs/OE1.v2/Data/blocks/blast_db/infasta_db' 
os.makedirs(infasta_db_dir, exist_ok=True)

## Make a reciprocal blast using blastn function from blast utils
infasta_db = makeblast_db(seqs_path=input_fasta, db_out=infasta_db_dir, remove_old_db=True)
blastn(blast_input_seqs=input_fasta, blast_db_file=infasta_db, blast_output_table_tsv=output_reciprocal_blast_tsv,
             num_threads=100, reward=2, gap_extend=2, gapopen=4, penalty=-2, word_size=15)


############################################ Preparing sequences for iterative block finding ###########################################

blast_df = pd.read_csv(output_reciprocal_blast_tsv, sep='\t')
blast_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']

blast_size = blast_df.shape[0]
print(f'Org_size: {blast_size}')
s0_df = blast_df[(blast_df['pident'] > 90) & (blast_df['length'] > 15)]
print(f'Selected pident 90%: {s0_df.shape[0]/blast_size}')
s1_df = alignment_absolute_start_end(s0_df)
s2_df = select_internal_aligns(df=s1_df,border=10) 

# Filter by query alignment interval
if params_select_internals_alns:
    s3_df = s2_df# blast.fuse_contained_intervals(df=s2_df, start_col='qstart', end_col='qend', group_cols=['qseqid', 'sseqid'], keep_highest_score_col='bitscore')
else:
    s3_df = s1_df

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

############################################ Iterative Block finding ###########################################

import warnings
warnings.filterwarnings('ignore')

import shutil
import pandas as pd
from utils.bedtools_utils import bedtools_getfasta, coords_to_regions, blast_to_bed
from utils.vsearch_utils import vsearch_dereplication, vsearch_clustersize, vsearch_sortbysize


iteration_dir = output_dir /  'block_iterations'
print(f"Iteration dir: {iteration_dir}")
if os.path.isdir(iteration_dir): 
    print("Removieng previous iteration dir")
    shutil.rmtree(iteration_dir)
iteration_dir.mkdir(exist_ok=True)

iteration = 0
num_blocks = []
while iteration < 2:
    
    print(f"\n###### ITERATION {iteration} ######\n")

    in_fasta                         = iteration_dir / f'{iteration}'     / 'centroids.fasta'
    blast_db                         = iteration_dir / f'{iteration + 1}' / 'blast_db'
    reciprocal_blast_tsv             = iteration_dir / f'{iteration + 1}' / 'blastn.reciprocal.tsv'
    coords_bed                       = iteration_dir / f'{iteration + 1}' / 'regions.bed'
    extracted_regions_fasta          = iteration_dir / f'{iteration + 1}' / 'regions.fasta'
    retrieve_coords_bed              = iteration_dir / f'{iteration}'     / 'retrieve' / 'regions.bed'
    retrieve_extracted_regions_fasta = iteration_dir / f'{iteration}'     / 'retrieve' / 'regions.fasta'
    derep_fasta                      = iteration_dir / f'{iteration + 1}' / 'derep.fasta'
    sorted_fasta                     = iteration_dir / f'{iteration + 1}' / 'sorted.fasta'
    centroids_fasta                  = iteration_dir / f'{iteration + 1}' / 'centroids.fasta'
    clusters_uc                      = iteration_dir / f'{iteration + 1}' / 'clusters.uc'
    
    iteration_path = os.path.dirname(clusters_uc)
    iteration_log =  os.path.join(iteration_path, 'logs')
    
    (clusters_uc.parent).mkdir( exist_ok=True)
    (retrieve_coords_bed.parent).mkdir(exist_ok=True, parents=True)
    
    ## Make a reciprocal blast using blastn function from blast utils
    if iteration == 0:
        blast_df = s3_df
        in_fasta = input_fasta
    else:
        blast_db = makeblast_db(seqs_path=in_fasta, db_out=blast_db, remove_old_db=True, log_file=iteration_log)
        
        blastn(blast_input_seqs=in_fasta,
               blast_db_file=blast_db,
               blast_output_table_tsv=reciprocal_blast_tsv,
               num_threads=80, reward=1, gap_extend=5, gapopen=5,
               penalty=-1, word_size=10, log_file=iteration_log)
        
        blast_df = pd.read_csv(reciprocal_blast_tsv, sep='\t', header=None)
        blast_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
        
    if iteration == 0:
        blast_df0 = blast_df.copy()

    blast_df = blast_df[blast_df['pident'] > 90]
    blast_df = blast_df[blast_df['length'] > 20]
    
    if params_retrive_sinlge_aln:
        # Retrive aligmentes that only map with itself
        retrive_df, blast_df = retrieve_single_alignments(blast_df)
        if retrive_df.shape[0] != 0:
            print(f"Retrieving sinlge mapped blocks iteration: {iteration}")
            retrive_bed_df = blast_to_bed(retrive_df)
            retrive_bed_df = retrive_bed_df.drop_duplicates()
            retrive_bed_df.to_csv(retrieve_coords_bed, sep='\t', header=False, index=False)
            #bedtools_getfasta(fasta_path=in_fasta, bed_path=retrieve_coords_bed, extracted_regions_path=retrieve_extracted_regions_fasta, show_command=True, log_file=iteration_log)
            get_fasta_from_bed(in_fasta, retrieve_coords_bed, retrieve_extracted_regions_fasta)
    else:
        print(f"Retrieving all blocks iteration: {iteration}")
        retrive_bed_df = blast_to_bed(blast_df)
        retrive_bed_df = retrive_bed_df.drop_duplicates()
        retrive_bed_df.to_csv(retrieve_coords_bed, sep='\t', header=False, index=False)
        #bedtools_getfasta(fasta_path=in_fasta, bed_path=retrieve_coords_bed, extracted_regions_path=retrieve_extracted_regions_fasta, show_command=True, log_file=iteration_log)
        get_fasta_from_bed(in_fasta, retrieve_coords_bed, retrieve_extracted_regions_fasta)
    
    print(f"Number of Retrived Seeds in iteration {iteration}: {len(retrive_bed_df)}")

        
    if blast_df.empty:
        print("Descomposition finished, all blocks are independent")
        break

    unnoised_df = blast_df.groupby('qseqid').apply(lambda x: unnoise_coords(x, radious=10))
    unnoised_df = unnoised_df.reset_index().drop(columns='level_1')

    # Filter coordinates with low number of instances 
    unnoised_df = unnoised_df[unnoised_df['instances'] > param_min_instances]
    
    if iteration == 0:
        unnoised_df0 = unnoised_df.copy()


    # Apply the function across groups
    get_fasta_df = unnoised_df.groupby('qseqid', group_keys=False).apply(coords_to_regions)
    get_fasta_df.to_csv(coords_bed, sep='\t', header=False, index=False)
    #bedtools_getfasta(fasta_path=in_fasta, bed_path=coords_bed, extracted_regions_path=extracted_regions_fasta, show_command=True, log_file=iteration_log)
    get_fasta_from_bed(in_fasta, coords_bed, extracted_regions_fasta)
    print(f"Extracted regions from bed file: {extracted_regions_fasta}")
    vsearch_dereplication(input_fasta=extracted_regions_fasta, derep_fasta=derep_fasta, min_seq_length=20, log_file=iteration_log)
    vsearch_sortbysize(derep_fasta=derep_fasta, sorted_fasta=sorted_fasta, log_file=iteration_log)
    vsearch_clustersize( sorted_fasta=sorted_fasta, centroids_fasta=centroids_fasta, uc_file=clusters_uc,query_cov=params_cluster_coverage,
                        target_cov=params_cluster_coverage, id_thresh=params_cluster_identity, min_seq_length=20, log_file=iteration_log,
                        use_both_strands=True,use_sizein=True)
    
    
    iteration += 1

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
