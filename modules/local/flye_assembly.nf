import groovy.json.JsonBuilder

// processes required for assembly using flye


process assembleCore_flye {
    errorStrategy = {task.attempt <= 4 ? 'retry' : 'ignore'}
    maxRetries 4
    label "wfplasmid"
    cpus params.threads
    memory "4GB"
    input:
        tuple val(sample_id), path(fastq), val(approx_size)
    output:
        tuple val(sample_id), path("${sample_id}.reconciled.fasta"), optional: true, emit: assembly
        tuple val(sample_id), path("${sample_id}.downsampled.fastq"), optional: true, emit: downsampled
        tuple val(sample_id), env(STATUS), emit: status
    script:
        cluster_dir = "trycycler/cluster_001"
        int coverage_target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 100
        int max_len = (approx_size as Integer) * 1.2
        int min_q = 7
        int exit_number = task.attempt <= 4 ? 1 : 0
        // min_overlap normally auto calculated but with a lower limit of 3000
        // assembly with same size as overlap will likely fail
        def min_overlap = approx_size.toInteger() <= 3000 ? '--min-overlap 1000' : ''
        def meta = params.non_uniform_coverage ? '--meta' : ''
        def seqkit_threads = params.threads >= 6 ? 2 : 1
    """
    # STATUS is changed after a command succeeds and is related to the following section
    # This is because we can't put comments in a multiline cmd
    ############################################################
    # Trimming
    ############################################################
    STATUS="Failed to trim reads"
    (seqkit subseq -j $seqkit_threads -r $params.trim_length:-$params.trim_length $fastq | \
        seqkit subseq -j $seqkit_threads -r 1:$max_len | \
        seqkit seq -j $seqkit_threads -m $min_len -Q $min_q -g > "${sample_id}.trimmed.fastq") && 
        

    ############################################################
    # Downsampling
    ############################################################

    STATUS="Failed to downsample reads" &&
    (rasusa \
        --coverage $coverage_target \
        --genome-size $approx_size \
        --input "${sample_id}.trimmed.fastq" > "${sample_id}.downsampled.fastq") &&

    ############################################################
    # Subsetting
    ############################################################
    STATUS="Failed to Subset reads" &&
    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads "${sample_id}.downsampled.fastq" \
        --out_dir sets \
        --genome_size $approx_size) && 

    ############################################################
    # Assembly
    ############################################################
    STATUS="Failed to assemble using Flye" &&
    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        flye \
            --${params.flye_quality} \${SUBSET} \
            --deterministic \
            --threads $task.cpus \
            --genome-size $approx_size \
            --out-dir "assm_\${SUBSET_NAME}" \
            ${meta} \
            $min_overlap 
             
        mv assm_sample_0*/assembly.fasta "assm_\${SUBSET_NAME}/\${SUBSET_NAME}_assembly.fasta" 
    done) && 

    ############################################################
    # Trim assemblies
    ############################################################
    STATUS="Failed to trim Assembly" &&
    (for assembly in \$(ls assm_sample_0*/*assembly.fasta)
    do  
        echo \$assembly
        assembly_name=\$(basename -s .fasta \$assembly)
        ass_stats=\$(dirname \$assembly)/assembly_info.txt
        deconcatenate.py \
            \$assembly \
            -o \${assembly_name}.deconcat.fasta \
            --approx_size $approx_size
    done
    ls *.deconcat.fasta > /dev/null 2>&1) &&


    ############################################################
    # Reconciliation
    ############################################################
    STATUS="Failed to reconcile assemblies" &&
    (trycycler cluster \
        --assemblies *.deconcat.fasta \
        --reads "${sample_id}.downsampled.fastq" \
        --out_dir trycycler) &&
    (trycycler reconcile \
        --reads "${sample_id}.downsampled.fastq" \
        --cluster_dir $cluster_dir \
        --max_trim_seq_percent 20 \
        --max_add_seq_percent 10) &&
    (trycycler msa --cluster_dir $cluster_dir) &&
    (trycycler partition --reads "${sample_id}.downsampled.fastq" --cluster_dirs $cluster_dir) &&
    (trycycler consensus --cluster_dir $cluster_dir)

    ############################################################
    # Exit handling
    ############################################################

    if [ ! -f "${cluster_dir}/7_final_consensus.fasta" ]; then
        if ls ${cluster_dir}/1_contigs/*.fasta 1> /dev/null 2>&1; then
            STATUS="Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > "${sample_id}.reconciled.fasta") \
                && echo "Trycycler failed, outputting un-reconciled assembly"
        elif [ "$exit_number" == "1" ]; then
            echo \$STATUS
            echo "Assembly failed, retrying process"
            exit 1
        elif [ "$exit_number" == "0" ]; then
            echo \$STATUS
            echo "Failed final attempt"
        fi
    else
        mv "${cluster_dir}/7_final_consensus.fasta" "${sample_id}.reconciled.fasta"
        STATUS="Completed successfully"
    fi
    """
}