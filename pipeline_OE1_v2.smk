# snakefile
import os, re

configfile: 'config.yaml'

PATH = os.getcwd()
RESULTS_D = config['Results.Folder'] 
S_PATH = os.path.join(PATH, 'Scripts')
V_PATH = os.path.join(RESULTS_D, 'Visuals') 
D_PATH = os.path.join(RESULTS_D, 'Data')

FILTERS_D    = 'filter_sizes'
JOIN_CREFS_D = 'join_cont_refs'
MAPC_D       = 'map_contamination'
DECONT_D     = 'decontamination'
ANALYZEC_D   = 'annalyze_contamination'
CHECK_GC     = 'check_gc'   
SREPEATS     = 'short_tandem_repeats'
ABC_D        = 'blast_to_abc'
MCL_D        = 'mcl_clustering'
SEEDS_D      = 'make_seeds'

task = {
    'filter_sizes.rule'             : (D_PATH, FILTERS_D, 'toxo1.size_filtered.fasta'),
    'join_cont_refs.rule'           : (D_PATH, JOIN_CREFS_D, 'join_crefs.fastq'),
    'map_contamination.rule'        : (D_PATH, MAPC_D, 'query_x_ref-contamination.paf'),
    'discard_contamination.rule'    : (D_PATH, DECONT_D, 'decont.fasta'),
    'GC_check.rule'                 : (D_PATH, CHECK_GC, 'decont.low_gc.fasta' ),
    'find_str.rule'                 : (D_PATH, SREPEATS, 'short_tandem_repeats.tsv'),
    'discard_str.rule'              : (D_PATH, SREPEATS, 'sequences_wo_str.fasta'),
    'blast_to_abc.rule'             : (D_PATH, ABC_D, 'blast_map.abc'),
    'mcl_clustering.rule'           : (D_PATH, MCL_D, 'seq_distances.mci'),
    'make_seeds.rule'               : (D_PATH, SEEDS_D, '.sentinel.make_seeds')
}

final_files = []
for rule, tuple_path in task.items():
    if config.get(rule):
        path = os.path.join(*tuple_path)
        final_files.append(path)

print('\nFinal Files:\n')
for file in final_files:
    print(file)

rule all:
    input:
        final_files        

rule filter_sizes:
    input: 
        fasta = config['query_sequences.Path'],
        fasta_sizes = os.path.join(S_PATH, 'fasta_sizes.py')
    params:
        min_read_length = config['filter_sizes.min_read_length']
    output: 
        hist_fig = os.path.join(V_PATH, FILTERS_D, 'hist.read_lengths.pdf' ),
        filtered_fasta = os.path.join(*task['filter_sizes.rule'])
    message: 
        """
        Discard any sequences from a query fasta file with a size below {params} nucleotides 
        """
    shell:
        """
        python3 "{input.fasta_sizes}" \
         "{input.fasta}" \
         "{params.min_read_length}" \
         "{output.hist_fig}" \
         "{output.filtered_fasta}"
        """

rule join_cont_refs:
    input:
        cont_ref_dir = config['cont_references.Folder']
    output: 
        join_refs = os.path.join(*task['join_cont_refs.rule'])
    message:
        """
        Joins all fastq files from cont_reference.folder into a single file
        """
    shell:
        """
        cat {input.cont_ref_dir}/* > {output.join_refs}
        """

rule map_contamination:
    input:
        join_refs = rules.join_cont_refs.output.join_refs,
        filtered_fasta = rules.filter_sizes.output.filtered_fasta
    output:
        aligns_paf = os.path.join(*task['map_contamination.rule'])
    threads: config['map_contamination.threads']
    message:
        """
        Maps the query fasta against the reference genomes from different contamintan organisms
        """
    shell:
        """
        minimap2 \
          {input.join_refs} \
          {input.filtered_fasta} \
          -x lr:hq \
          -O 6,36 \
          -E 6,4 \
          --secondary=no \
          -t {threads} > {output.aligns_paf}
        """

rule discard_contamination:
    input:
        discard_contamination_py = os.path.join(S_PATH, 'discard_contamination.py'),
        filtered_fasta = rules.filter_sizes.output.filtered_fasta,
        aligns_paf = rules.map_contamination.output.aligns_paf
    params: 
        min_coverage = config['discard_contamination.min_coverage']
    output:
        decont_fasta = os.path.join(*task['discard_contamination.rule']),
        decont_stats = os.path.join(D_PATH, DECONT_D, 'decont.tsv')
    message:
        """
        Given a PAF alignment table, this rule discards all sequences from the 
        filtered FASTA that have a good mapping against any contamination references.
        """
    shell:
        """
        python {input.discard_contamination_py} \
          --fasta_in {input.filtered_fasta} \
          --paf_in {input.aligns_paf} \
          --fasta_out {output.decont_fasta} \
          --stats_out {output.decont_stats} \
          --min_coverage {params.min_coverage}
        """

rule check_GC:
    input: 
        decont_fasta = rules.discard_contamination.output.decont_fasta
    params: 
        max_gc_content = config['GC_check.max_GC']
    output: 
        gc_table = os.path.join(D_PATH, CHECK_GC, 'input_gc_table.txt'),
        decont_low_gc_fasta = os.path.join(D_PATH, CHECK_GC, 'decont.low_gc.fasta'),
        decont_high_gc_fasta = os.path.join(D_PATH, CHECK_GC, 'decont.high_gc.fasta'),
        gc_stats = os.path.join(D_PATH, CHECK_GC, 'gc_stats.tsv')
    shell:
        """
        # Step 1: Generate GC table for input FASTA
        seqkit fx2tab -n -g {input.decont_fasta} > {output.gc_table}

        # Step 2: Split into low-GC and high-GC outputs
        # Low-GC sequences
        awk -v max_gc={params.max_gc_content} '$2 < max_gc {{print $1}}' {output.gc_table} \
            | seqkit grep -f - {input.decont_fasta} > {output.decont_low_gc_fasta}
        # High-GC sequences
        awk -v max_gc={params.max_gc_content} '$2 >= max_gc {{print $1}}' {output.gc_table} \
            | seqkit grep -f - {input.decont_fasta} > {output.decont_high_gc_fasta}

        # Step 3: Compute input FASTA stats
        input_num=$(wc -l < {output.gc_table})
        input_avg=$(awk '{{total += $2}} END {{printf "%.2f", (NR>0 ? total/NR : 0)}}' {output.gc_table})
        input_std=$(awk '{{sum += $2; sum2 += $2^2}} END {{printf "%.2f", (NR>0 ? sqrt(sum2/NR - (sum/NR)^2) : 0)}}' {output.gc_table})

        # Step 4: Compute low-GC output stats
        tmp_low_gc=$(mktemp)
        seqkit fx2tab -n -g {output.decont_low_gc_fasta} > $tmp_low_gc
        low_num=$(wc -l < $tmp_low_gc)
        low_avg=$(awk '{{total += $2}} END {{printf "%.2f", (NR>0 ? total/NR : 0)}}' $tmp_low_gc)
        low_std=$(awk '{{sum += $2; sum2 += $2^2}} END {{printf "%.2f", (NR>0 ? sqrt(sum2/NR - (sum/NR)^2) : 0)}}' $tmp_low_gc)
        rm $tmp_low_gc

        # Step 5: Compute high-GC output stats
        tmp_high_gc=$(mktemp)
        seqkit fx2tab -n -g {output.decont_high_gc_fasta} > $tmp_high_gc
        high_num=$(wc -l < $tmp_high_gc)
        high_avg=$(awk '{{total += $2}} END {{printf "%.2f", (NR>0 ? total/NR : 0)}}' $tmp_high_gc)
        high_std=$(awk '{{sum += $2; sum2 += $2^2}} END {{printf "%.2f", (NR>0 ? sqrt(sum2/NR - (sum/NR)^2) : 0)}}' $tmp_high_gc)
        rm $tmp_high_gc

        # Step 6: Write stats to file
        echo -e "file_name\tnum_seqs\tavg_gc_content\tmean_gc_content\tstd_deviation" > {output.gc_stats}
        echo -e "$(basename {input.decont_fasta})\t$input_num\t$input_avg\t$input_avg\t$input_std" >> {output.gc_stats}
        echo -e "$(basename {output.decont_low_gc_fasta})\t$low_num\t$low_avg\t$low_avg\t$low_std" >> {output.gc_stats}
        echo -e "$(basename {output.decont_high_gc_fasta})\t$high_num\t$high_avg\t$high_avg\t$high_std" >> {output.gc_stats}
        """

rule find_str:
    input:
        decont_low_gc_fasta=rules.check_GC.output.decont_low_gc_fasta,
        trf_script=os.path.join(S_PATH, "trf")
    params:
        trf_params="2 7 7 80 10 50 2000"
    output:
        str_table=os.path.join(D_PATH, SREPEATS, 'short_tandem_repeats.tsv')
    shell:
        """
        mkdir -p $(dirname {output})
        cd $(dirname {output})
        {input.trf_script} {input.decont_low_gc_fasta} {params.trf_params} -h -d > log.trf 2>&1 || true
        ls

        trf_params_name=$(echo "{params.trf_params}" | sed 's/ /./g')
        input_fasta_name=$(basename {input.decont_low_gc_fasta})

        if [ -s "$input_fasta_name.$trf_params_name.dat" ]; then
            tail -n +6 "$input_fasta_name.$trf_params_name.dat" | \
            grep -vE '^$|Parameters' | \
            awk 'BEGIN{{RS="Sequence: "; FS="\\n"; OFS="\\t"}} ($2 != "" ){{print $0}}' | \
            awk 'BEGIN{{RS="\\n\\n" ; FS="\\n"}} {{for(i=2;i<=NF;i++) printf("%s %s\\n", $1, $i)}}' | \
            awk 'BEGIN{{OFS="\\t" ; print "seqid", "rstart", "rend", "period_size", "copy_number", "consensus_size", "percent_matches", "percent_indels", "score", "A", "C", "G", "T", "entropy(0-2)", "repeat_sequence"}} {{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}}' > {output}

        else
            echo -e "seqid\\trstart\\trend\\tperiod_size\\tcopy_number\\tconsensus_size\\tpercent_matches\\tpercent_indels\\tscore\\tA\\tC\\tG\\tT\\tentropy(0-2)\\trepeat_sequences" > {output}
        fi
        """

rule discard_str:
    input:
        discard_str_py      = os.path.join(S_PATH, 'discard_str.py'),
        decont_low_gc_fasta = rules.check_GC.output.decont_low_gc_fasta,
        str_table           = rules.find_str.output.str_table
    params: 
        period_size = 3.0,
        total_tandem_size = 150,
        tandem_coverage = config['discard_str.tandem_coverage']
    output:
        passed_sequences = os.path.join(D_PATH, SREPEATS, 'sequences_wo_str.fasta'),
        discarded_sequences = os.path.join(D_PATH, SREPEATS, 'short_tandem_repeats.fasta'),
        repeats = os.path.join(D_PATH, SREPEATS, 'repeats.fasta'),
        str_table = os.path.join(D_PATH, SREPEATS, 'filtered_str.tsv')
    shell:
        """
            python3 {input.discard_str_py} \
             --input-fasta {input.decont_low_gc_fasta} \
             --input-str-table {input.str_table} \
             --param-period-size {params.period_size} \
             --param-tandem-coverage {params.tandem_coverage} \
             --param-total-tandem-coverage {params.total_tandem_size} \
             --output-passed-sequences {output.passed_sequences} \
             --output-discarded-sequences {output.discarded_sequences} \
             --output-short-repeats {output.repeats}  \
             --output-str-table {output.str_table}
        """

rule blast_to_abc:
    input:
        blast_to_abc_py     = os.path.join(S_PATH, 'blast_to_abc.py'),
        seqs_wo_str_fasta   = rules.discard_str.output.passed_sequences
    params: 
        min_aln_identity = config['blast_to_abc.min_aln_identity'],
        min_aln_length = config['blast_to_abc.min_aln_length'],
        min_blast_word_size = config['blast_to_abc.min_blast_word_size']
    threads: config['blast_to_abc.threads']
    output:
        blast_tsv                = os.path.join(D_PATH, ABC_D, 'blastn.seqs_wo_str.inner.tsv'),
        coverage_graph_abc       = os.path.join(D_PATH, ABC_D, 'blast_map.abc'),
        uninformative_seqs_fasta = os.path.join(D_PATH, ABC_D, 'uninformative_sequences.fasta')
    shell:
        """
            python3 {input.blast_to_abc_py} \
             --input-fasta {input.seqs_wo_str_fasta} \
             --param-min-aln-identity {params.min_aln_identity} \
             --param-min-aln-length {params.min_aln_length} \
             --param-min-blast-word-size {params.min_blast_word_size} \
             --param-threads {threads} \
             --output-blast-table {output.blast_tsv} \
             --output-coverage-graph-abc {output.coverage_graph_abc} \
             --output-uninformative-seqs-list {output.uninformative_seqs_fasta}  
        """

rule mcl_clustering:
    input:
        mcl_clustering_py  = os.path.join(S_PATH, 'mcl_clustering.py'),
        coverage_graph_abc = rules.blast_to_abc.output.coverage_graph_abc,
        sequences_wo_str   = rules.discard_str.output.passed_sequences
    params: 
        inflation = config['mcl_clustering.inflation']
    output:
        blast_names_tab       = os.path.join(D_PATH, MCL_D, 'seq_names.list'),
        matrix_mci            = os.path.join(D_PATH, MCL_D, 'seq_distances.mci'),
        cluster_dump          = os.path.join(D_PATH, MCL_D, 'dump.mci'),
        barplot_cluster_sizes = os.path.join(V_PATH, MCL_D, 'barplot.cluster_sizes.svg'),
        cluster_seqs_dir      = directory(os.path.join(D_PATH, MCL_D, 'sequneces_clusters'))
    message:
        """
        Clusterize a set sequences with MCL using the total coverage between them.
        """
    shell:
        """
        python3 {input.mcl_clustering_py} \
         --input-graph-abc {input.coverage_graph_abc} \
         --input-fasta {input.sequences_wo_str} \
         --param-inflation {params.inflation} \
         --output-blast-names-tab {output.blast_names_tab} \
         --output-matrix-mci {output.matrix_mci} \
         --output-cluster-dump {output.cluster_dump} \
         --output-clust-fasta-dir {output.cluster_seqs_dir} \
         --output-barplot-cluster-sizes {output.barplot_cluster_sizes}
        """

rule make_seeds:
    input:
        make_seeds_py   = os.path.join(S_PATH, "make_seeds.py"),
        fasta_dir = rules.mcl_clustering.output.cluster_seqs_dir,
    params:
        min_fasta_size  = config['make_seeds.min_fasta_size'],
        kmer_size = config['make_seeds.kmer_size'],
        min_kmer_counts = config['make_seeds.min_kmer_counts']
    threads:
        config['make_seeds.threads']
    output:
        seeds_dir = directory(os.path.join(D_PATH, SEEDS_D))
        sentinel_make_seeds = os.path.join(D_PATH, SEEDS_D, '.sentinel.make_seeds')
    shell:
        """
        python3 {input.make_seeds_py} \
         --input-fasta-dir {input.fasta_dir} \
         --params-min-fasta-size {params.min_fasta_size} \
         --params-kmer-size {params.kmer_size} \
         --params-kmer-counts {params.min_kmer_counts} \
         --params-threads {threads} \
         --output-seeds-dir {output.seeds_dir}

        touch {output.sentinel_make_seeds}
        """
        





#rule validate_sequences_before:
#    input: 
#        continuos_sequence_aligment_py = os.path.join(S_PATH, 'continuos_sequence_aligment.py'),
#        query_sequences = rules.filter_sizes.output.filtered_fasta,
#        sequence_blocks = config['validation.sequence_blocks']
#    output:
#        sequence_stats = os.path.join(D_PATH, SEQ_VALIDATION, 'before.continuos_mit-block_coverage.tsv')
#    message:
#        """
#        Given a set of short sequence blocks, this script produces a TSV table with all adyacent aligments  
#        that can be made when mapping to a set of sequences from a fasta file.   
#        """
#    shell:
#        """
#
#        """



## Maybe the rules map_contamination, analyze contamination, and filter contamination could be generalyzed 
## to map, analyze and filter target and contamined reads (or map, analyze and filter between multiple targets)

# Work in progress
# Pipeline Bottle Neck: Possibly this rule requires the greatest processing power of this pipeline

#print(f'Filtered Fasta: {rules.filter_sizes.output.filtered_fasta}')
#print(f'Cont_refs: {os.path.join(REFS_D)}')
#
#ref_genome_wildcards = glob_wildcards(os.path.join(REFS_D,'{fname}.fna')).fname
#print(f'ref_genome_wildcards: {ref_genome_wildcards}')
#print(expand(os.path.join(D_PATH, MAPC_D, '{fname}.sam'), fname=ref_genome_wildcards))
#print(f'cont_sam = os.path.join(D_PATH, MAPC_D, '{fname}.sam')')


#rule map_refs:
#    input:
#        filtered_fasta = rules.filter_sizes.output.filtered_fasta,
#        ref_genome = os.path.join(REFS_D, '{fname}.fna')
#    params:
#        cont_paf = os.path.join(D_PATH, MAPC_D)
#    output:
#        cont_paf = os.path.join(D_PATH, MAPC_D, '{fname}.paf')
#    shell:
#        """
#        mkdir -p {params.cont_paf}
#        minimap2 {input.ref_genome} {input.filtered_fasta} > {output.cont_paf}
#
#        """
#
#rule map_refs_mito:
#    input:
#        mit_frags = os.path.join(PATH, 'Blocks_Namasivayam.fa'),
#        ref_genome = os.path.join(REFS_D, '{fname}.fna')
#    params:
#        cont_paf = os.path.join(D_PATH, MAPC_D)
#    output:
#        cont_paf = os.path.join(D_PATH, MAPC_D, '{fname}.mito.paf')
#    shell:
#        """
#        mkdir -p {params.cont_paf}
#        minimap2  -k15 -w5 -A2 -B8 -O20,50 -E5,2 -r2k -Y \
#         --secondary=no  \
#         {input.ref_genome} {input.mit_frags} > {output.cont_paf}
#
#        """
#
#rule blast_db:
#    input:
#        ref_genome = os.path.join(REFS_D, '{fname}.fna')
#    params:
#        blast_db = os.path.join(D_PATH, MAPC_D, 'blast_db', '{fname}_db'),
#        blast_db_dir = os.path.join(D_PATH, MAPC_D, 'blast_db')
#    output:
#        blast_db = os.path.join(D_PATH, MAPC_D, 'blast_db', '{fname}_db.nsq')
#    shell:
#        """
#        mkdir -p {params.blast_db_dir}
#        makeblastdb -in {input.ref_genome} -dbtype nucl -out {params.blast_db}
#        """
#
#rule blast_map:
#    input:
#        filtered_fasta = os.path.join(D_PATH, FILTERS_D, 'toxo1.size_filtered.fasta'),
#        blast_db = rules.blast_db.output.blast_db
#    params:
#        blast_db = rules.blast_db.params.blast_db,
#        word_size = 11
#    output:
#        cont_tsv = os.path.join(D_PATH, MAPC_D, '{fname}.tsv')
#    shell:
#        """
#        blastn \
#         -query {input.filtered_fasta} \
#         -db {params.blast_db} \
#         -out {output.cont_tsv} \
#         -outfmt 6 \
#         -word_size {params.word_size}
#         -threads 4
#        """

#rule decont_sequences:
#    input:
#        filtered_fasta =
#
#rule gc_discard:
#    input:  
#        decont_fasta = os.path.join()
        

# Collect all .sam files
#cont_sam_files = glob_wildcards(os.path.join(D_PATH, MAPC_D, '{fname}.sam')).fname

# Rule to analyze contamination
#rule analyze_contamination:
#    input:
#        cont_sam = expand(os.path.join(D_PATH, MAPC_D, '{fname}.sam'), fname=cont_sam_files)
#    output:
#        decont_stats = os.path.join(D_PATH, ANALYZEC_D, 'decont-stats.tsv')
#    run:
#        # Creating a table where each row represents a different sequence, and each column shows 
#        # the mapping stats on each reference.
#        with open(output.decont_stats, 'w') as output_file:
#            output_file.write('Example output')
#
#cont_sam_files = glob_wildcards(os.path.join(D_PATH, '{fname}.sam')).fname
#FILTERC_D = 'filter_contamination'
#rule filter_contamination:
#    input:
#        filtered_fasta = rules.filter_sizes.output.filtered_fasta,
#        cont_sam = expand(os.path.join(D_PATH, MAPC_D, '{fname}.sam'), fname=cont_sam_files)
#    output:
#        decont_fasta = os.path.join(D_PATH, FILTERC_D, 'decont.fasta'),
#        decont_stats = os.path.join(D_PATH, FILTERC_D, 'decont-stats.tsv')
#    run:
#        ...
