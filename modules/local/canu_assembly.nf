import groovy.json.JsonBuilder

// processes required for assembly using canu


process assembleCore_canu {
    errorStrategy = {task.attempt <= 4 ? 'retry' : 'ignore'}
    maxRetries 4
    label "canu"
    cpus params.threads
    memory "7GB"
    input:
        tuple val(meta), path(fastq)
    output:
        tuple val(meta), path("${meta.alias}.reconciled.fasta"), optional: true, emit: assembly
        tuple val(meta), path("${meta.alias}.downsampled.fastq"), optional: true, emit: downsampled
        tuple val(meta.alias), env(STATUS), emit: status
    script:
        String cluster_dir = "trycycler/cluster_001"
        int coverage_target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 100
        int max_len = (meta.approx_size as Integer) * 1.2
        int exit_number = task.attempt <= 4 ? 1 : 0
        def fast = params.canu_fast == true ? '-fast' : ''
        // WSL does not support named pipes used by Canu, setting these parameters avoids their use
        def windows_params = System.properties['os.version'].toLowerCase().contains("wsl") ? """\
        -mhapPipe=false \
        -purgeOverlaps=false \
        -saveOverlaps=true """ : ""
        def seqkit_threads = params.threads >= 6 ? 2 : 1

    """
    ############################################################
    # Trimming
    ############################################################
    STATUS="Failed to trim reads"
    (
        if [[ $params.trim_length -gt 0 ]]; then
            seqkit subseq -j $seqkit_threads -r $params.trim_length:-$params.trim_length $fastq
        else
            cat $fastq
        fi \
        | seqkit subseq -j $seqkit_threads -r 1:$max_len \
        | seqkit seq -j $seqkit_threads -m $min_len -Q $params.min_quality -g > "${meta.alias}.trimmed.fastq"
    ) &&

    ############################################################
    # Downsampling
    ############################################################

    STATUS="Failed to downsample reads" &&
    (rasusa \
        --coverage $coverage_target \
        --genome-size ${meta.approx_size} \
        --input "${meta.alias}.trimmed.fastq" > "${meta.alias}.downsampled.fastq") &&
        

    ############################################################
    # Subsetting
    ############################################################
    STATUS="Failed to Subset reads" &&
    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads "${meta.alias}.downsampled.fastq" \
        --out_dir sets \
        --genome_size ${meta.approx_size}) && 

    ############################################################
    # Assembly
    ############################################################
    STATUS="Failed to assemble using Canu" &&
    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        canu \
            -p \$SUBSET_NAME \
            -d assm_\${SUBSET_NAME} \
            -maxThreads=$task.cpus \
            genomeSize=${meta.approx_size} \
            $fast \
            -nanopore \$SUBSET \
            $windows_params
    done) && 

    ############################################################
    # Trim assemblies
    ############################################################
    STATUS="Failed to trim Assembly" &&
    (for assembly in \$(ls assm_sample_0*/*.contigs.fasta)
    do  
        echo \$assembly
        assembly_name=\$(basename -s .fasta \$assembly)
        trim.py \
            \$assembly \
            -o \${assembly_name}.trimmed.fasta
        deconcatenate.py \
            \${assembly_name}.trimmed.fasta \
            -o \${assembly_name}.deconcat.fasta \
            --approx_size ${meta.approx_size}
    done
    ls *.deconcat.fasta > /dev/null 2>&1) &&

    ############################################################
    # Reconciliation
    ############################################################
    STATUS="Failed to reconcile assemblies" &&
    (trycycler cluster \
        --assemblies *.deconcat.fasta \
        --reads "${meta.alias}.downsampled.fastq" \
        --out_dir trycycler) &&
    (trycycler reconcile \
        --reads "${meta.alias}.downsampled.fastq" \
        --cluster_dir $cluster_dir \
        --max_trim_seq_percent 20 \
        --max_add_seq_percent 10) &&
    (trycycler msa --cluster_dir $cluster_dir) &&
    (trycycler partition --reads "${meta.alias}.downsampled.fastq" --cluster_dirs $cluster_dir) &&
    (trycycler consensus --cluster_dir $cluster_dir)

    ############################################################
    # Exit handling
    ############################################################

    if [ ! -f "${cluster_dir}/7_final_consensus.fasta" ]; then
        if ls ${cluster_dir}/1_contigs/*.fasta 1> /dev/null 2>&1; then
            STATUS="Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > "${meta.alias}.reconciled.fasta") \
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
        mv "${cluster_dir}/7_final_consensus.fasta" "${meta.alias}.reconciled.fasta"
        STATUS="Completed successfully"
    fi
    """
}