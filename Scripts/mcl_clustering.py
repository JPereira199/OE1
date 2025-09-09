#################################################### Loading libraries and script arguments ###################################################

import matplotlib.pyplot as plt 
from Bio import SeqIO
import numpy as np
import seaborn as sns
import subprocess 
import argparse
import os

test=False
if test:
    print('Warning! Developer mode activated! ')
    input_graph_abc = '/home/jpereira/OEs/OE1.v2/Data/mcl_processing/blast_map.abc'
    input_fasta = '/home/jpereira/OEs/OE1.v2/Data/short_tandem_repeats/sequences_wo_str.fasta'
    param_inflation = float(1.4)
    output_names_tab =  '/home/jpereira/OEs/OE1.v2/Data/mcl_processing/blast_map.names.tab'
    output_matrix_mci = '/home/jpereira/OEs/OE1.v2/Data/mcl_processing/blast_map.mci'
    output_cluster_dump = '/home/jpereira/OEs/OE1.v2/Data/mcl_processing/dump.blast_map.mci'
    output_clust_fasta_dir = '/home/jpereira/OEs/OE1.v2/Data/mcl_processing/sequences_cluster'
    output_barplot_cluster_sizes = '/home/jpereira/OEs/OE1.v2/Visuals/mcl_processing/barplot.cluster_sizes.svg'
    
else:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-graph-abc', type=str, help='A table representing a graph in ABC format')
    parser.add_argument('--input-fasta', type=str, help='Sequences used to define the sequences clusters')
    parser.add_argument('--param-inflation', type=float)
    parser.add_argument('--output-blast-names-tab', type=str, help='Internal file created by mcxload. A file with the names of each sequence')
    parser.add_argument('--output-matrix-mci', type=str, help='Internal file created by mcl. A matrix in format mci with the distances between sequences')
    parser.add_argument('--output-cluster-dump', type=str, help='A file with the members of each cluster. Each cluster has its member in the same row.')
    parser.add_argument('--output-clust-fasta-dir', type=str, help='Directory with the clustericed sequences')
    parser.add_argument('--output-barplot-cluster-sizes', type=str)
    args = parser.parse_args()

    input_graph_abc              = args.input_graph_abc
    input_fasta                  = args.input_fasta
    param_inflation              = args.param_inflation
    output_names_tab             = args.output_blast_names_tab
    output_matrix_mci            = args.output_matrix_mci
    output_cluster_dump          = args.output_cluster_dump
    output_clust_fasta_dir       = args.output_clust_fasta_dir
    output_barplot_cluster_sizes = args.output_barplot_cluster_sizes

os.makedirs(os.path.dirname(output_names_tab), exist_ok=True)
os.makedirs(os.path.dirname(output_barplot_cluster_sizes), exist_ok=True)
os.makedirs(output_clust_fasta_dir, exist_ok=True)

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

################################################## Creating distance matrix from an ABC graph #################################################

command = (f"mcxload" 
           f" -abc {input_graph_abc}"
           f" --stream-mirror"
           f" -re max" 
           f" -write-tab {output_names_tab}"
           f" -o {output_matrix_mci}"
           f" ")

try:
    subprocess.run(command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    print('Error in running mcxload: ', e)
    print(f"Command {command}")

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

########################################################## Clustering graph using MCL #########################################################

# Compute the rounded inflation value and construct output paths
rounded_i = round(param_inflation * 10)
mcl_output_name = f"out.{os.path.basename(output_matrix_mci)}.I{rounded_i}"
mcl_output_path = os.path.join(os.path.dirname(output_matrix_mci), mcl_output_name)

command = (f"mcl"
           f" {output_matrix_mci}"
           f" -I {param_inflation}"
           f" -o {mcl_output_path}")

try:
    subprocess.run(command, shell=True, check=True)
except subprocess.CalledProcessError as e:
    print('Error in running mcl: ', e)

# Construct the mcxdump command, ensuring proper spacing
command = (
    f"mcxdump "
    f"-icl {mcl_output_path} "
    f"-tabr {output_names_tab} "
    f"-o {output_cluster_dump}"
)

try:
    subprocess.run(command, shell=True, check=True)
    print(f"mcxdump completed successfully. Output: {output_cluster_dump}")
except subprocess.CalledProcessError as e:
    print('Error in running mcxdump:', e)

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

############################################################ Saving sequences of each cluster ############################################################

# Load all sequences into a dictionary {record.id: record}
seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, format='fasta'))

with open(output_cluster_dump, 'r') as file:
    clust_num = 0

    for line in file:
        line = line.strip()
        clust_ids = line.split('\t')  # IDs in this cluster

        clust_num += 1
        output_path = os.path.join(output_clust_fasta_dir, f'cluster_{clust_num}.fasta')

        # Filter the sequences that belong to this cluster
        cluster_seqs = [seq_dict[seq_id] for seq_id in clust_ids if seq_id in seq_dict]

        # Write them as a FASTA
        with open(output_path, 'w') as out_fasta:
            SeqIO.write(cluster_seqs, out_fasta, 'fasta')


##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

############################################################ Ploting cluster sizes ############################################################

sns.set_theme(style='whitegrid')

# Assuming that output_cluster_dump and param_inflation are defined
clust_sizes = []
clust_names = []

with open(output_cluster_dump, 'r') as file:
    for line in file:
        line = line.strip()  # Remove whitespace and newlines
        split_line = line.split('\t')
        clust_names.append(split_line)
        clust_sizes.append(len(split_line))

# Use a range of integers for x-axis: cluster numbers from 1 to n
x_array = np.arange(1, len(clust_sizes) + 1)
y_array = clust_sizes

plt.figure(figsize=(9, 4))
plt.bar(x_array, y_array, edgecolor='black')

plt.xlabel('Cluster Number', fontsize=12)
plt.ylabel('Cluster Size', fontsize=12)
plt.title(f'Graph Clusterization by MCL (inflation={param_inflation})', fontsize=14)

plt.xticks(ticks=x_array, labels=x_array)

plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.grid(axis='x', visible=False)    
plt.tight_layout()
plt.savefig(output_barplot_cluster_sizes)
if test:
    plt.show()
else:
    plt.close()

##### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -