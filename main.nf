#!/usr/bin/env nextflow

import groovy.json.JsonBuilder 
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
include { start_ping; end_ping } from './lib/ping'


def groupIt(ch) {
    return ch.map { it -> return tuple(it.simpleName, it) }.groupTuple()
}

def nameIt(ch) {
    return ch.map { it -> return tuple(it.simpleName, it) }
}

process combineFastq {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(approx_size), val(type)
    output:
        path "${sample_id}.fastq.gz", optional: true, emit: sample
        path "${sample_id}.stats", optional: true, emit: stats
        env STATUS, emit: status
        tuple val(sample_id), val(approx_size), emit: approx_size 
    script:
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8
        def expected_length_max = approx_size.toInteger()
        def expected_length_min = approx_size.toInteger()
        int max = (expected_length_max.toInteger()) * 1.5
        int min = (expected_length_min.toInteger()) * 0.5
    """
    STATUS=${sample_id}",Failed due to insufficient reads"
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} > /dev/null
    fastcat -a "$min" -b "$max" -s ${sample_id} -r ${sample_id}.interim -x ${directory} > ${sample_id}.fastq
    if [[ "\$(wc -l <"${sample_id}.interim")" -ge "$value" ]];  then 
        gzip ${sample_id}.fastq
        STATUS=${sample_id}",Completed successfully"
    fi
    """
}


process filterHostReads {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    input:
        file fastq
        file reference
        file regions_bedfile
    output:
        path "*.filtered.fastq", optional: true, emit: unmapped
        path "*.stats", optional: true, emit: host_filter_stats
        env STATUS, emit: status
    script:
        def name = fastq.simpleName
        def regs = regions_bedfile.name != 'NO_REG_BED' ? regs : 'none'
    """
    STATUS="${name},Failed due to filtered host reads"
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
    fi ) && STATUS="${name},Completed successfully"
    """
}


process assembleCore {
    errorStrategy = {task.attempt <= 5 ? 'retry' : 'ignore'}
    maxRetries 5
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_name), file(fastq), val(approx_size)
    output:
        path "*.reconciled.fasta", optional: true, emit: assembly
        path "*.downsampled.fastq", optional: true, emit: downsampled
        env STATUS, emit: status
    script:
        name = fastq.simpleName
        cluster_dir = "trycycler/cluster_001"
        int target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 1000
        int max_len = approx_size.toInteger() * 1.2
        int min_q = 7
        int exit_number = task.attempt <= 5 ? 1 : 0
        def fast = params.fast == true ? '-fast' : ''
        def cluster_option = (params.canu_useGrid == false) ? """\
        -useGrid=false \
        -obtovlThreads=$params.threads \
        -utgovlThreads=$params.threads \
        -corThreads=$params.threads \
        -redThreads=$params.threads \
        -batThreads=$params.threads """ : ""

    """
    
    ############################################################
    # Trimming
    ############################################################
    STATUS="${name},Failed to trim reads"
    (seqkit subseq -j $params.threads -r $params.trim_length:-$params.trim_length $fastq | \
        seqkit subseq -j $params.threads -r 1:$max_len | \
        seqkit seq -j $params.threads -m $min_len -Q $min_q -g > ${name}.trimmed.fastq) \
        && STATUS="${name},Failed to downsample reads" &&

    ############################################################
    # Downsampling
    ############################################################
    
    
    (rasusa \
        --coverage $target \
        --genome-size $approx_size \
        --input ${name}.trimmed.fastq > ${name}.downsampled.fastq) \
        && STATUS="${name},Failed to Subset reads" &&

    ############################################################
    # Subsetting
    ############################################################
    
    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads ${name}.downsampled.fastq \
        --out_dir sets \
        --genome_size $approx_size) \
        && STATUS="${name},Failed to assemble using Canu" &&
    
    ############################################################
    # Assembly
    ############################################################

    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        canu \
            -p \$SUBSET_NAME \
            -d assm_\${SUBSET_NAME} \
            -maxThreads=$params.threads \
            genomeSize=$approx_size \
            $fast \
            -nanopore \$SUBSET \
            $cluster_option
    done) && STATUS="${name},Failed to trim Assembly" &&

    ############################################################
    # Trim assemblies
    ############################################################

    (for ASSEMBLY in \$(ls assm_*/*.contigs.fasta)
    do
        ASSEMBLY_NAME=\$(basename -s .fasta \$ASSEMBLY)
        trim.py \
            \$ASSEMBLY \
            -o \${ASSEMBLY_NAME}.trimmed.fasta
        deconcatenate.py \
            \${ASSEMBLY_NAME}.trimmed.fasta \
            -o \${ASSEMBLY_NAME}.deconcat.fasta
    done
    ls *.deconcat.fasta 1> /dev/null 2>&1) \
    && STATUS="${name},Failed to reconcile assemblies" &&

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
            STATUS=${name}",Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > ${name}.reconciled.fasta) \
                && echo "Trycycler failed, outputting un-reconciled assembly"
        elif [ "$exit_number" == "1" ]; then
            echo "Assembly failed, retrying process"
            exit 1
        elif [ "$exit_number" == "0" ]; then
            echo "Failed final attempt"
        fi
    else
        mv ${cluster_dir}/7_final_consensus.fasta ${name}.reconciled.fasta
        STATUS=${name}",Completed successfully"
    fi
    """
}


process medakaPolishAssembly {
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample), file(draft), file(fastq)
    output:
        path "*.final.fasta", emit: polished
        env STATUS, emit: status
    """
    STATUS=${fastq.simpleName}",Failed to polish assembly with Medaka"
    medaka_consensus -i $fastq -d $draft -m r941_min_high_g360 -o . -t $task.cpus -f
    echo ">${fastq.simpleName}" >> ${fastq.simpleName}.final.fasta
    sed "2q;d" consensus.fasta >> ${fastq.simpleName}.final.fasta
    STATUS=${fastq.simpleName}",Completed successfully"
    """
}


process downsampledStats {
    label "wfplasmid"
    cpus 1 
    input:
        file sample
    output:
        path "*.stats", optional: true
    """
    fastcat -s ${sample.simpleName} -r ${sample.simpleName}.downsampled $sample > /dev/null
    if [[ "\$(wc -l <"${sample.simpleName}.downsampled")" -ge "2" ]];  then 
        mv ${sample.simpleName}.downsampled ${sample.simpleName}.stats
    fi
    """
}


process assemblyMafs {
    label "wfplasmid"
    cpus 1
    input:
        file assemblies
    output:
        path "*.maf", emit: assembly_maf
    shell:
    '''
    # Get maf files for dotplots
    for assm in !{assemblies}
    do
        lastdb ${assm}.lastdb $assm
        lastal ${assm}.lastdb $assm > ${assm}.maf
    done
    '''
}

process findPrimers {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    input:
        file primers
        file sequence 
    output:
        path "*.bed", optional: true 
    shell:
    '''
    cat !{sequence} | seqkit amplicon -p !{primers} -m 3 -j !{params.threads} --bed >> !{sequence.simpleName}.interim
    if [[ "$(wc -l <"!{sequence.simpleName}.interim")" -ge "1" ]];  then 
        mv !{sequence.simpleName}.interim !{sequence.simpleName}.bed
    fi
    '''
}

process getVersions {
    label "wfplasmid"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    medaka --version | sed 's/ /,/' >> versions.txt
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
    cpus params.threads
    input:
        path annotation_database
        path "assemblies/*"
    output:    
        path "feature_table.txt", emit: feature_table 
        path "plannotate.json", emit: json
        path "*annotations.bed", emit: annotations
        path "plannotate_report.json", emit: report
    """
    if [ -e "assemblies/OPTIONAL_FILE" ]; then
        assemblies=""
    else 
        assemblies="--sequences assemblies/"
    fi 
    run_plannotate.py \$assemblies --database $annotation_database
    """
}


process inserts {
    label "wfplasmid"
    cpus params.threads
    input:
         path "primer_beds/*"
         path "assemblies/*"
         file align_ref
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
    find_inserts.py \$inserts $ref
    """
}


process report {
    label "wfplasmid"
    cpus 1
    input:
        path "downsampled_stats/*"
        file final_status
        path "per_barcode_stats/*"
        path "host_filter_stats/*"
        path "versions/*"
        path "params.json"
        file plannotate_json
        file inserts_json
    output:
        path "wf-clone-validation-*.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "inserts/*", optional: true, emit: inserts
    script:
        report_name = "wf-clone-validation-" + params.report_name + '.html'
    """
    report.py \
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
      
        sample_fastqs = combineFastq(samples)
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
        sample_with_size = nameIt(samples_filtered).join(sample_fastqs.approx_size)
        // Core assembly and reconciliation
        assemblies = assembleCore(sample_with_size)

        named_drafts = groupIt(assemblies.assembly)
        named_samples = groupIt(assemblies.downsampled)
        named_drafts_samples = named_drafts.join(named_samples)


        // Polish draft assembly
        polished = medakaPolishAssembly(named_drafts_samples)

        final_status = sample_fastqs.status
            .join(assemblies.status,remainder: true)
            .join(polished.status,remainder: true)
            .join(updated_status,remainder: true)
            .collectFile(name: 'final_status.csv', newLine: true)

        downsampled_stats = downsampledStats(
            assemblies.downsampled)
        
        assembly_mafs = assemblyMafs(
            polished.polished.collect())
        
        primer_beds = findPrimers(primers, polished.polished)
        software_versions = getVersions()
        workflow_params = getParams()

        annotation = runPlannotate(
            database, polished.polished.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")))
        
        insert = inserts(primer_beds.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            polished.polished.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            align_ref)
        report = report(
            downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status,
            sample_fastqs.stats.collect(),
            filtered_stats,
            software_versions.collect(),
            workflow_params,
            annotation.report,
            insert.json)
        
        results = polished.polished.concat(
            report.html,
            report.sample_stat,
            annotation.feature_table,
            insert.inserts,
            annotation.json,
            annotation.annotations)
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
    // Start ping
    start_ping()
    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq,
        params.min_barcode, params.max_barcode, params.approx_size)
   
    host_reference = file(params.host_reference, type: "file")
    regions_bedfile = file(params.regions_bedfile, type: "file")
    primer_file = file("$projectDir/data/OPTIONAL_FILE")
    if (params.primers != null){
        primer_file = file(params.primers, type: "file")
    }
    align_ref = file("$projectDir/data/OPTIONAL_FILE")
    if (params.reference != null){
        align_ref = file(params.reference, type: "file")
    }
    database = file(params.db_directory, type: "dir", checkIfExists:true)
    
    // Run pipeline
    results = pipeline(samples, host_reference, regions_bedfile, database, primer_file, align_ref)

    output(results[0])
    // End ping
    end_ping(pipeline.out.telemetry)
}
