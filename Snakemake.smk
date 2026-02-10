# snakefile
import os, re
import shlex

configfile: 'config.yaml'

PATH = os.getcwd()
RESULTS_D = config['Results.Folder'] 
S_PATH = os.path.join(PATH, 'Scripts')
V_PATH = os.path.join(RESULTS_D, 'Visuals') 
D_PATH = os.path.join(RESULTS_D, 'Data')
L_PATH = os.path.join(RESULTS_D, 'Logs')

FILTERS_D    = 'filter_sizes'
CLEAN_N_FILTER_D = 'clean_n_filter'
JOIN_CREFS_D = 'join_cont_refs'
MAPC_D       = 'map_contamination'
DECONT_D     = 'decontamination'
ANALYZEC_D   = 'annalyze_contamination'
CHECK_GC     = 'check_gc'   
SREPEATS     = 'short_tandem_repeats'
ABC_D        = 'blast_to_abc'
MCL_D        = 'mcl_clustering'
SEEDS_D      = 'make_seeds'
BLOCKS_EXT_D = 'block_extension'
DECONT_STATS_D = "decont_stats"
DEFINE_BLOCKS_D = "define_blocks"
DOTPLOT_D = "blastn_dotplot"


task = {
    'filter_sizes.rule'             : (D_PATH, FILTERS_D, 'reads.size_filtered.fasta'),
    'clean_n_filter.rule'           : (D_PATH, CLEAN_N_FILTER_D, 'reads.size_filtered.fasta'),
    'join_cont_refs.rule'           : (D_PATH, JOIN_CREFS_D, 'join_crefs.fastq'),
    'map_contamination.rule'        : (D_PATH, MAPC_D, 'query_x_ref-contamination.paf'),
    'discard_contamination.rule'    : (D_PATH, DECONT_D, 'decont.fasta'),
    'GC_check.rule'                 : (D_PATH, CHECK_GC, 'decont.low_gc.fasta' ),
    'find_str.rule'                 : (D_PATH, SREPEATS, 'short_tandem_repeats.tsv'),
    'discard_str.rule'              : (D_PATH, SREPEATS, 'sequences_wo_str.fasta'),
    'blast_to_abc.rule'             : (D_PATH, ABC_D, 'blast_map.abc'),
    'mcl_clustering.rule'           : (D_PATH, MCL_D, 'seq_distances.mci'),
    'make_seeds.rule'               : (D_PATH, SEEDS_D, '.sentinel.make_seeds'),
    'block_extension.rule'          : (D_PATH, BLOCKS_EXT_D, '.sentinel.block_extension'),
    'decont_stats.rule'             : (D_PATH, DECONT_STATS_D, "summary.tsv"),
    'define_blocks.rule'            : (D_PATH, DEFINE_BLOCKS_D, ".sentinel.define_blocks"),
    'blastn_dotplot.rule'           : (D_PATH, DOTPLOT_D, "dotplot_table.tsv")
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


rule clean_n_filter:
    input:
        reads = config["query_sequences.Path"]
        clean_n_filter_py = os.path.join(S_PATH, "clean_n_filter.py")
    output:
        filtered_fasta = os.path.join(D_PATH, "clean_n_filter", "reads.size_filtered.fasta"),
        fig   = os.path.join(V_PATH, "clean_n_filter", "hist.read_lengths.svg")
    params:
        min_len = config["clean_n_filter.min_read_length"],
        min_q   = config.get("clean_n_filter.min_mean_quality", 10)
    threads: config.get('clean_n_filter.threads', 10)
    log:
        os.path.join(L_PATH, "clean_n_filter.log")
    conda:
        "envs/clean_n_filter.yaml"
    message:
        """
        Filtering reads by:
          - mean base quality (FASTQ only)
          - minimum length
        Generating length histograms before and after filtering.
        """
    shell:
        """
        set -euo pipefail

        python {input.clean_n_filter_py} \
          --input {input.reads} \
          --min_len {params.min_len} \
          --min_q {params.min_q} \
          --out_fasta {output.filtered_fasta} \
          --out_fig {output.fig} \
          > {log} 2>&1
        """
   

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
        filtered_fasta = rules.clean_n_filter.output.filtered_fasta
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
         -x map-ont \
         -c \
         --secondary=no \
         -t {threads} \
         -o {output.aligns_paf} \
         {input.join_refs} \
         {input.filtered_fasta}
        """
    
        # """
        # minimap2 \
        #   {input.join_refs} \
        #   {input.filtered_fasta} \
        #   -x lr:hq \
        #   -O 6,36 \
        #   -E 6,4 \
        #   --secondary=no \
        #   -t {threads} > {output.aligns_paf}
        # """

rule discard_contamination:
    input:
        discard_contamination_py = os.path.join(S_PATH, "discard_contamination2.py"),
        filtered_fasta = rules.filter_sizes.output.filtered_fasta,
        aligns_paf = rules.map_contamination.output.aligns_paf
    params:
        min_alignment_length = config.get("discard_contamination.min_alignment_length", 700),
        max_de = config.get("discard_contamination.max_de", 0.30),
        min_query_coverage = config.get("discard_contamination.min_query_coverage", 0.70)
    output:
        decont_fasta = os.path.join(*task["discard_contamination.rule"]),
        decont_stats = os.path.join(D_PATH, DECONT_D, "decont.tsv"),
        decont_align_stats = os.path.join(D_PATH, DECONT_D, "alignment-decont.tsv")
    message:
        """
        Given a PAF alignment table, discard sequences from the filtered FASTA
        that have a "good" mapping against any contamination references, based on:
          - min alignment length
          - max de
          - min query coverage
        """
    shell:
        """
        python {input.discard_contamination_py} \
          --fasta_in {input.filtered_fasta} \
          --paf_in {input.aligns_paf} \
          --fasta_out {output.decont_fasta} \
          --stats_out {output.decont_stats} \
          --align_stats_out {output.decont_align_stats} \
          --min_alignment_length {params.min_alignment_length} \
          --max_de {params.max_de} \
          --min_query_coverage {params.min_query_coverage}
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
            echo -e "seqid\\trstart\\trend\\tperiod_size\\tcopy_number\\tconsensus_size\\tpercent_matches\\tpercent_indels\\tscore\\tA\\tC\\tG\\tT\\tentropy(0-2)\\trepeat_sequence" > {output}
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
        cluster_seqs_dir      = directory(os.path.join(D_PATH, MCL_D, 'sequneces_clusters')),
        cluster_1_fasta       = os.path.join(D_PATH, MCL_D, 'sequneces_clusters', 'cluster_1.fasta')
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
# Defining subjects to make the barplots of decont stats


NAM_SUBJECTS = []
NAM_LABELS = []
NAM_LABELS_QUOTED = ""

if config.get("decont_stats.rule"):

    if config.get("filter_sizes.rule"):
        NAM_LABELS.append("Input")
        NAM_SUBJECTS.append(os.path.join(*task["filter_sizes.rule"]))

    if config.get("discard_contamination.rule"):
        NAM_LABELS.append("Decontamination")
        NAM_SUBJECTS.append(os.path.join(*task["discard_contamination.rule"]))

    if config.get("GC_check.rule"):
        NAM_LABELS.append("GC filter")
        NAM_SUBJECTS.append(os.path.join(D_PATH, CHECK_GC, "decont.low_gc.fasta"))

    if config.get("discard_str.rule"):
        NAM_LABELS.append("STR filter")
        NAM_SUBJECTS.append(os.path.join(D_PATH, SREPEATS, "sequences_wo_str.fasta"))

    # Optional: add a cluster fasta (because your mcl step outputs a directory, not a tracked fasta)
    if config.get("mcl_clustering.rule"):
        NAM_LABELS.append("MCL filter")
        NAM_SUBJECTS.append(os.path.join(D_PATH, MCL_D, 'sequneces_clusters', 'cluster_1.fasta'))

    NAM_LABELS_QUOTED = " ".join(shlex.quote(x) for x in NAM_LABELS)


rule decont_stats:
    input:
        script = os.path.join(S_PATH, "decont_stats.py"),
        #blast_utils = os.path.join(PATH, "utils", "blast_utils.py"),
        validation_seqs = config["validation.sequence_blocks"],
        subjects = NAM_SUBJECTS
    output:
        summary_tsv = os.path.join(D_PATH, DECONT_STATS_D, "summary.tsv"),
        barplots_svg = os.path.join(D_PATH, DECONT_STATS_D, "summary_barplots.svg")
    threads:
        config.get("decont_stats.threads", 8)
    params:
        out_dir = os.path.join(D_PATH, DECONT_STATS_D),
        word_size = config.get("decont_stats.word_size", 11),
        min_pident = config.get("decont_stats.min_pident", 90.0),
        nam_threshold = config.get("decont_stats.nam_threshold", 0.50),
        labels = NAM_LABELS_QUOTED
    log:
        os.path.join(L_PATH, "decont_stats.log")
    message:
        """
        Computing decontamination statistics and NAM-block coverage across pipeline steps.
        Generating per-step coverage plots and a global summary barplot.
        """
    conda:
        "envs/decont_stats.yaml"
    shell:
        """
        set -euo pipefail

        MPLBACKEND=Agg \
        python {input.script} \
          --validation_seqs {input.validation_seqs} \
          --subjects {input.subjects} \
          --subject_labels {params.labels} \
          --out_dir {params.out_dir} \
          --threads {threads} \
          --word_size {params.word_size} \
          --min_pident {params.min_pident} \
          --nam_threshold {params.nam_threshold} \
          #--blast_utils_path {{input.blast_utils}} \
          > {log} 2>&1
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
        seeds_dir = directory(os.path.join(D_PATH, SEEDS_D)),
        seeds_cluster_1_tsv = os.path.join(D_PATH, SEEDS_D, "cluster_1.tsv"),
        sentinel_make_seeds = os.path.join(D_PATH, SEEDS_D, '.sentinel.make_seeds')
    shell:
        """
        python3 {input.make_seeds_py} \
         --input-fasta-dir {input.fasta_dir} \
         --params-min-fasta-size {params.min_fasta_size} \
         --params-kmer-size {params.kmer_size} \
         --params-min-kmer-counts {params.min_kmer_counts} \
         --params-threads {threads} \
         --output-seeds-dir {output.seeds_dir}

        touch {output.sentinel_make_seeds}
        """
        
rule block_extension:
    input:
        block_extension_py   = os.path.join(S_PATH, "block_extension.py"),
        reads_fasta = rules.mcl_clustering.output.cluster_1_fasta,
        seeds_tsv = rules.make_seeds.output.seeds_cluster_1_tsv,
    params:
        min_hits =                    config['block_extension.min_hits'],
        max_aln_extension =           config['block_extension.max_aln_extension'],
        seed_map_coverage_threshold = config['block_extension.seed_map_coverage_threshold'],
        window_radius =               config['block_extension.window_radius'],
        signal_threshold =            config['block_extension.signal_threshold'],
        aln_w_step_size =             config['block_extension.aln_w_step_size'],
        try_mean_signal =             config['block_extension.try_mean_signal'],
        mean_signal_threshold =       config['block_extension.mean_signal_threshold'],
        plot =                        config['block_extension.plot']
    threads:
        config['block_extension.threads']
    output:
        block_extension_data_dir = directory(os.path.join(D_PATH, BLOCKS_EXT_D)),
        block_extension_vis_dir = directory(os.path.join(V_PATH, BLOCKS_EXT_D)),
        blocks_fasta = os.path.join(D_PATH, BLOCKS_EXT_D, "blocks.fasta"),
        sentinel_block_extension = os.path.join(D_PATH, BLOCKS_EXT_D, '.sentinel.block_extension')
    shell:
        """
        python3 {input.block_extension_py} \
         --input-reads-fasta {input.reads_fasta} \
         --input-seeds-tsv {input.seeds_tsv} \
         --params-threads {threads} \
         --params-min-hits {params.min_hits} \
         --params-max-aln-extension {params.max_aln_extension} \
         --params-seed-map-coverage-threshold {params.seed_map_coverage_threshold} \
         --params-window-radius {params.window_radius} \
         --params-signal-threshold {params.signal_threshold} \
         --params-aln-w-step-size {params.aln_w_step_size} \
         --params-try-mean-signal {params.try_mean_signal} \
         --params-mean-signal-threshold {params.mean_signal_threshold} \
         --params-plot {params.plot} \
         --output-data-dir {output.block_extension_data_dir} \
         --output-vis-dir {output.block_extension_vis_dir} 
         
        touch {output.sentinel_block_extension}
        """






rule define_blocks:
    input:
        script = os.path.join(S_PATH, "define_blocks.py"),
        blocks_fasta = rules.block_extension.output.blocks_fasta,
        sentinel_block_extension = rules.block_extension.output.sentinel_block_extension
    output:
        out_dir = directory(os.path.join(D_PATH, DEFINE_BLOCKS_D)),
        reciprocal_tsv = os.path.join(D_PATH, DEFINE_BLOCKS_D, "blastn.infasta_reciprocal.tsv"),
        blocks_fasta = os.path.join(D_PATH, DEFINE_BLOCKS_D, "block_iterations", "1", "retrieve", "regions.fasta"),
        sentinel = os.path.join(D_PATH, DEFINE_BLOCKS_D, ".sentinel.define_blocks")
    threads:
        config.get("define_blocks.threads", 20)
    params:
        min_instances = config.get("define_blocks.min_instances", 0),
        cluster_identity = config.get("define_blocks.cluster_identity", 0.85),
        cluster_coverage = config.get("define_blocks.cluster_coverage", 0.85),
        select_internals = config.get("define_blocks.select_internals_alns", False),
        retrieve_single = config.get("define_blocks.retrieve_single_aln", False),
    log:
        os.path.join(L_PATH, DEFINE_BLOCKS_D, "define_blocks.log")
    message:
        """
        Defining blocks by reciprocal BLAST + iterative clustering (define_blocks).
        Input: {input.blocks_fasta}
        """
    #conda:
    #    "envs/define_blocks.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})

        python3 {input.script} \
          --input-fasta {input.blocks_fasta} \
          --param-min-instances {params.min_instances} \
          --param-threads {threads} \
          --params-cluster-identity {params.cluster_identity} \
          --params-cluster-coverage {params.cluster_coverage} \
          --params-select-internals-alns {params.select_internals} \
          --params-retrive-sinlge-aln {params.retrieve_single} \
          --output-dir {output.out_dir} \
          > {log} 2>&1

        test -s {output.reciprocal_tsv}
        touch {output.sentinel}
        """

rule blastn_dotplot_unmapped:
    input:
        script = os.path.join(S_PATH, "blastn_dotplot.py"),
        regions_fasta = rules.define_blocks.output.blocks_fasta,
        fasta_queries = config["validation.sequence_blocks"],  # e.g. /home/.../Blocks_Namasivayam.fa
    output:
        dotplot_svg = os.path.join(D_PATH, DOTPLOT_D, "dotplot.svg"),
        blast_tsv = os.path.join(D_PATH, DOTPLOT_D, "blastn.fasta1_vs_fasta2.tsv"),
        dotplot_table = os.path.join(D_PATH, DOTPLOT_D, "dotplot_table.tsv"),
        unmapped_queries = os.path.join(D_PATH, DOTPLOT_D, "unmapped_queries.txt"),
        unmapped_subjects = os.path.join(D_PATH, DOTPLOT_D, "unmapped_subjects.txt"),
        unmapped_summary = os.path.join(D_PATH, DOTPLOT_D, "unmapped_summary.tsv"),
        output_base_dir = directory(os.path.join(D_PATH, DOTPLOT_D)),
    threads:
        config.get("blastn_dotplot.threads", 8)
    params:
        title = config["blastn_dotplot.title"],
        word_size = config.get("blastn_dotplot.word_size", 15),
        qnames = config.get("blastn_dotplot.qnames", True),
        snames = config.get("blastn_dotplot.snames", True),
        annot_size = config.get("blastn_dotplot.annot_size", 10),
    log:
        os.path.join(L_PATH, "blastn_dotplot.log")
    message:
        "Running BLASTN dotplot + unmapped ID report for {params.title}"
    shell:
        """
        mkdir -p $(dirname {log})

        python3 {input.script} \
          --fasta1 {input.regions_fasta} \
          --fasta2 {input.fasta_queries} \
          --output_base_dir {output.output_base_dir} \
          --word_size {params.word_size} \
          --qnames {params.qnames} \
          --snames {params.snames} \
          --annot_size {params.annot_size} \
          --title {params.title} \
          > {log} 2>&1
        """
