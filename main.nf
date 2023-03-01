#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'


process combineFastq {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    input:
        tuple val(meta), path(directory), val(approx_size)
    output:
        tuple val(meta.sample_id), path("${meta.sample_id}.fastq.gz"), val(approx_size), optional: true, emit: sample
        path "${meta.sample_id}.stats", optional: true, emit: stats
        tuple val(meta.sample_id), env(STATUS), emit: status
    script:
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8
        def expected_length_max = approx_size.toInteger()
        def expected_length_min = approx_size.toInteger()
        int max = (expected_length_max.toInteger()) * 1.5
        int min = (expected_length_min.toInteger()) * 0.5
    """
    STATUS="Failed due to insufficient reads"
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} > /dev/null
    fastcat -a "$min" -b "$max" -s ${meta.sample_id} -r ${meta.sample_id}.interim -x ${directory} > ${meta.sample_id}.fastq
    if [[ "\$(wc -l <"${meta.sample_id}.interim")" -ge "$value" ]];  then
        gzip ${meta.sample_id}.fastq
        STATUS="Completed successfully"
    fi
    """
}


process filterHostReads {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_id), path(fastq), val(approx_size)
        path reference
        path regions_bedfile
    output:
        tuple val(sample_id), path("*.filtered.fastq"), val(approx_size), optional: true, emit: unmapped
        path "*.stats", optional: true, emit: host_filter_stats
        tuple val(sample_id), env(STATUS), emit: status
    script:
        def name = sample_id
        def regs = regions_bedfile.name != 'NO_REG_BED' ? regs : 'none'
    """
    STATUS="Failed due to filtered host reads"
    (minimap2 -t $task.cpus -y -ax map-ont $reference $fastq \
        | samtools sort -o ${name}.sorted.aligned.bam -
    samtools index ${name}.sorted.aligned.bam 
    samtools view -b -f 4  ${name}.sorted.aligned.bam > unmapped.bam
    samtools view -b -F 4  ${name}.sorted.aligned.bam > mapped.bam
    samtools fastq unmapped.bam > ${name}.filtered.fastq
    fastcat -s ${name} -r ${name}.interim ${name}.filtered.fastq > /dev/null
    if [[ "\$(wc -l <"${name}.stats")" -ge "1" ]];  then
        mv ${name}.interim ${name}.stats
    fi
    if [[ -f "$regs" ]]; then
        bedtools intersect -a mapped.bam -b $regs -wa \
            | samtools view -bh - > retained.bam
        samtools fastq retained.bam >> ${name}.filtered.fastq
    fi ) && STATUS="Completed successfully"
    """
}


process assembleCore {
    errorStrategy = {task.attempt <= 4 ? 'retry' : 'ignore'}
    maxRetries 4
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_id), path(fastq), val(approx_size)
    output:
        tuple val(sample_id), path("*.reconciled.fasta"), optional: true, emit: assembly
        tuple val(sample_id), path("*.downsampled.fastq"), optional: true, emit: downsampled
        tuple val(sample_id), env(STATUS), emit: status
    script:
        name = sample_id
        cluster_dir = "trycycler/cluster_001"
        int target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 100
        int max_len = approx_size.toInteger() * 1.2
        int min_q = 7
        int exit_number = task.attempt <= 4 ? 1 : 0
        def fast = params.fast == true ? '-fast' : ''
        def cluster_option = (params.canu_useGrid == false) ? """\
        -useGrid=false \
        -obtovlThreads=$task.cpus \
        -utgovlThreads=$task.cpus \
        -corThreads=$task.cpus \
        -redThreads=$task.cpus \
        -batThreads=$task.cpus """ : ""
        def windows_params = System.properties['os.version'].toLowerCase().contains("wsl") ? """\
        -mhapPipe=false \
        -purgeOverlaps=false \
        -saveOverlaps=true """ : ""

    """

    ############################################################
    # Trimming
    ############################################################
    STATUS="Failed to trim reads"
    (seqkit subseq -j $task.cpus -r $params.trim_length:-$params.trim_length $fastq | \
        seqkit subseq -j $task.cpus -r 1:$max_len | \
        seqkit seq -j $task.cpus -m $min_len -Q $min_q -g > ${name}.trimmed.fastq) \
        && STATUS="Failed to downsample reads" &&

    ############################################################
    # Downsampling
    ############################################################


    (rasusa \
        --coverage $target \
        --genome-size $approx_size \
        --input ${name}.trimmed.fastq > ${name}.downsampled.fastq) \
        && STATUS="Failed to Subset reads" &&

    ############################################################
    # Subsetting
    ############################################################

    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads ${name}.downsampled.fastq \
        --out_dir sets \
        --genome_size $approx_size) \
        && STATUS="Failed to assemble using Canu" &&

    ############################################################
    # Assembly
    ############################################################

    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        canu \
            -p \$SUBSET_NAME \
            -d assm_\${SUBSET_NAME} \
            -maxThreads=$task.cpus \
            genomeSize=$approx_size \
            $fast \
            -nanopore \$SUBSET \
            $cluster_option \
            $windows_params
    done) && STATUS="Failed to trim Assembly" &&

    ############################################################
    # Trim assemblies
    ############################################################

    (for ASSEMBLY in \$(ls assm_*/*.contigs.fasta)
    do
        ASSEMBLY_NAME=\$(basename -s .fasta \$ASSEMBLY)
        workflow-glue trim \
            \$ASSEMBLY \
            -o \${ASSEMBLY_NAME}.trimmed.fasta
        workflow-glue deconcatenate \
            \${ASSEMBLY_NAME}.trimmed.fasta \
            -o \${ASSEMBLY_NAME}.deconcat.fasta
    done
    ls *.deconcat.fasta 1> /dev/null 2>&1) \
    && STATUS="Failed to reconcile assemblies" &&

    ############################################################
    # Reconciliation
    ############################################################

    (trycycler cluster \
        --assemblies *.deconcat.fasta \
        --reads ${name}.downsampled.fastq \
        --out_dir trycycler) &&
    (trycycler reconcile \
        --reads ${name}.downsampled.fastq \
        --cluster_dir $cluster_dir \
        --max_trim_seq_percent 20 \
        --max_add_seq_percent 10) &&
    (trycycler msa --cluster_dir $cluster_dir) &&
    (trycycler partition --reads ${name}.downsampled.fastq --cluster_dirs $cluster_dir) &&
    (trycycler consensus --cluster_dir $cluster_dir)

    ############################################################
    # Exit handling
    ############################################################

    if [ ! -f "${cluster_dir}/7_final_consensus.fasta" ]; then
        if ls ${cluster_dir}/1_contigs/*.fasta 1> /dev/null 2>&1; then
            STATUS="Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > ${name}.reconciled.fasta) \
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
        mv ${cluster_dir}/7_final_consensus.fasta ${name}.reconciled.fasta
        STATUS="Completed successfully"
    fi
    """
}

process lookup_medaka_model {
    label "wfplasmid"
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model lookup_table '!{basecall_model}')
    echo $medaka_model
    '''
}

process medakaPolishAssembly {
    label "medaka"
    cpus params.threads
    input:
        tuple val(sample_id), path(draft), path(fastq), val(medaka_model)
    output:
        tuple val(sample_id), path("*.final.fasta"), emit: polished
        tuple val(sample_id), env(STATUS), emit: status
    script:
        def model = medaka_model
    """
    STATUS="Failed to polish assembly with Medaka"
    medaka_consensus -i "${fastq}" -d "${draft}" -m "${model}" -o . -t $task.cpus -f
    echo ">${sample_id}" >> "${sample_id}.final.fasta"
    sed "2q;d" consensus.fasta >> "${sample_id}.final.fasta"
    STATUS="Completed successfully"
    """
}


process downsampledStats {
    label "wfplasmid"
    cpus 1
    input:
        tuple val(sample_id), path(sample)
    output:
        path "*.stats", optional: true
    """
    fastcat -s ${sample_id} -r ${sample_id}.downsampled $sample > /dev/null
    if [[ "\$(wc -l <"${sample_id}.downsampled")" -ge "2" ]];  then
        mv ${sample_id}.downsampled ${sample_id}.stats
    fi
    """
}


process findPrimers {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    input:
        path primers
        tuple val(sample_id), path(sequence)
    output:
        path "*.bed", optional: true
    shell:
    '''
    cat !{sequence} | seqkit amplicon -p !{primers} -m 3 -j !{task.cpus} --bed >> !{sample_id}.interim
    if [[ "$(wc -l <"!{sample_id}.interim")" -ge "1" ]];  then
        mv !{sample_id}.interim !{sample_id}.bed
    fi
    '''
}

process medakaVersion {
    label "medaka"
    output:
        path "medaka_version.txt"
    """
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}

process getVersions {
    label "wfplasmid"
    cpus 1
    input:
        path "input_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    cat "input_versions.txt" >> "versions.txt"
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    trycycler --version | sed 's/ /,/' >> versions.txt
    porechop --version | sed 's/^/porechop,/'  >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    canu -version | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    last --version | sed 's/ /,/' >> versions.txt
    rasusa --version | sed 's/ /,/' >> versions.txt
    python -c "import spoa; print(spoa.__version__)" | sed 's/^/spoa,/'  >> versions.txt
    """
}


process getParams {
    label "wfplasmid"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process runPlannotate {
    label "wfplasmid"
    cpus 1
    input:
        path annotation_database
        path "assemblies/*"
        path final_status
    output:
        path "feature_table.txt", emit: feature_table
        path "plannotate.json", emit: json
        path "*annotations.bed", optional: true, emit: annotations
        path "plannotate_report.json", emit: report
    script:
        def database =  annotation_database.name.startsWith('OPTIONAL_FILE') ? "Default" : "${annotation_database}"
    """
    if [ -e "assemblies/OPTIONAL_FILE" ]; then
        assemblies=""
    else
        assemblies="--sequences assemblies/"
    fi
    workflow-glue run_plannotate \$assemblies --database $database
    """
}


process inserts {
    label "wfplasmid"
    cpus 1
    input:
         path "primer_beds/*"
         path "assemblies/*"
         path align_ref
    output:
        path "inserts/*", optional: true, emit: inserts
        path "*.json", emit: json
    script:
        def ref =  align_ref.name.startsWith('OPTIONAL_FILE') ? '' : "--reference ${align_ref}"
    """
    if [ -e "primer_beds/OPTIONAL_FILE" ]; then
        inserts=""
    else
        inserts="--primer_beds primer_beds/*"
    fi
    workflow-glue find_inserts \$inserts $ref  
    """
}


process report {
    label "wfplasmid"
    cpus 1
    input:
        path "downsampled_stats/*"
        path final_status
        path "per_barcode_stats/*"
        path "host_filter_stats/*"
        path "versions/*"
        path "params.json"
        path plannotate_json
        path inserts_json
        path lengths
    output:
        path "wf-clone-validation-*.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "inserts/*", optional: true, emit: inserts
    script:
        report_name = "wf-clone-validation-report.html"
    """
    workflow-glue report \
    --downsampled_stats downsampled_stats/* \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --status $final_status \
    --per_barcode_stats per_barcode_stats/* \
    --host_filter_stats host_filter_stats/* \
    --params params.json \
    --versions versions \
    --report_name $report_name \
    --plannotate_json $plannotate_json \
    --lengths $lengths \
    --inserts_json $inserts_json
    """
}

workflow pipeline {
    take:
        samples
        host_reference
        regions_bedfile
        database
        primers
        align_ref

    main:
        // Combine fastq from each of the sample directories into
        // a single per-sample fastq file
        named_samples = samples.map { it -> return tuple(it[1],it[0])}
        if(params.approx_size_sheet != null) {
            approx_size = Channel.fromPath(params.approx_size_sheet) \
            | splitCsv(header:true) \
            | map { row-> tuple(row.sample_id, row.approx_size) }
            final_samples = named_samples.join(approx_size)}
        else {
            final_samples = samples.map  { it -> return tuple(it[1],it[0], params.approx_size)}
        }
        sample_fastqs = combineFastq(final_samples)
        // Optionally filter the data, removing reads mapping to
        // the host or background genome
        if (host_reference.name != "NO_HOST_REF") {
            filtered = filterHostReads(
                    sample_fastqs.sample, host_reference, regions_bedfile)
            samples_filtered = filtered.unmapped
            updated_status = filtered.status
            filtered_stats = filtered.host_filter_stats.collect()
                             .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        }
        else {
            samples_filtered = sample_fastqs.sample
            updated_status = sample_fastqs.status
            filtered_stats = file("$projectDir/data/OPTIONAL_FILE")
        }
       
        // Core assembly and reconciliation
        assemblies = assembleCore(samples_filtered)
        
        named_drafts = assemblies.assembly.groupTuple()
        named_samples = assemblies.downsampled.groupTuple()
        named_drafts_samples = named_drafts.join(named_samples)

        if(params.medaka_model) {
            log.warn "Overriding Medaka model with ${params.medaka_model}."
            medaka_model = Channel.fromList([params.medaka_model])
        }
        else {
            // map basecalling model to medaka model
            lookup_table = Channel.fromPath("${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_model = lookup_medaka_model(lookup_table, params.basecaller_cfg)
        }
        // Polish draft assembly
        polished = medakaPolishAssembly(named_drafts_samples.combine(medaka_model))
       
        // Concat statuses and keep the last of each
        final_status = sample_fastqs.status.concat(updated_status)
        .concat(assemblies.status).concat(polished.status).groupTuple()
        .map { it -> it[0].toString() + ',' + it[1][-1].toString() }
        final_status = final_status.collectFile(name: 'final_status.csv', newLine: true)
    
        downsampled_stats = downsampledStats(assemblies.downsampled)

        primer_beds = findPrimers(primers, polished.polished)
        medaka_version = medakaVersion()
        software_versions = getVersions(medaka_version)
        workflow_params = getParams()

        annotation = runPlannotate(
            database, polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status)

        insert = inserts(primer_beds.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            align_ref)
        report = report(
            downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status,
            sample_fastqs.stats.collect(),
            filtered_stats,
            software_versions.collect(),
            workflow_params,
            annotation.report,
            insert.json,
            annotation.json)

        results = polished.polished.map { it -> it[1] }.concat(
            report.html,
            report.sample_stat,
            annotation.feature_table,
            insert.inserts,
            annotation.json,
            annotation.annotations,
            workflow_params)
    emit:
        results
        telemetry = workflow_params
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfplasmid"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: {
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "min_barcode":params.min_barcode,
        "max_barcode":params.max_barcode])
    host_reference = params.host_reference ?: 'NO_HOST_REF'
    host_reference = file(host_reference, checkIfExists: host_reference == 'NO_HOST_REF' ? false : true)
    regions_bedfile = params.regions_bedfile ?: 'NO_REG_BED'
    regions_bedfile = file(regions_bedfile, checkIfExists: regions_bedfile == 'NO_REG_BED' ? false : true)
    primer_file = file("$projectDir/data/OPTIONAL_FILE")
    if (params.primers != null){
        primer_file = file(params.primers, type: "file")
    }
    align_ref = file("$projectDir/data/OPTIONAL_FILE")
    if (params.reference != null){
        align_ref = file(params.reference, type: "file")
    }
    database = file("$projectDir/data/OPTIONAL_FILE")
    if (params.db_directory != null){
         database = file(params.db_directory, type: "dir")

    }

    // Run pipeline
    results = pipeline(samples, host_reference, regions_bedfile, database, primer_file, align_ref)

    output(results[0])
   
}


if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
