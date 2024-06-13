#!/usr/bin/env nextflow

import groovy.json.JsonBuilder

nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'
if (params.assembly_tool == 'canu') {
    include { assembleCore_canu as assembleCore } from './modules/local/canu_assembly.nf'
} else {
    include { assembleCore_flye as assembleCore } from './modules/local/flye_assembly.nf'
}

process checkIfEnoughReads {
    label 'wfplasmid'
    cpus params.threads
    memory '2GB'
    input:
        tuple val(meta),
            path('input.fastq.gz'),
            path('per-read-stats.tsv.gz'),
            val(approx_size)
        val extra_args
    output:
        tuple val(meta.alias), path("${meta.alias}.fastq.gz"), val(approx_size),
            optional: true, emit: sample
        path "${meta.alias}.stats.gz", emit: stats
        tuple val(meta.alias), env(STATUS), emit: status
    script:
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8
        int bgzip_threads = task.cpus == 1 ? 1 : task.cpus - 1
    """
    STATUS="Failed due to insufficient reads"
    mv per-read-stats.tsv.gz ${meta.alias}.stats.gz
    fastcat -s ${meta.alias} -r ${meta.alias}.interim $extra_args input.fastq.gz \
    | bgzip -@ $bgzip_threads > interim.fastq.gz
    if [[ "\$(wc -l < "${meta.alias}.interim")" -ge "$value" ]]; then
        mv interim.fastq.gz ${meta.alias}.fastq.gz
        STATUS="Completed successfully"
    fi
    """
}

process filterHostReads {
    errorStrategy 'ignore'
    label 'wfplasmid'
    cpus params.threads
    memory '4GB'
    input:
        tuple val(sample_id), path(fastq), val(approx_size)
        path reference
        path regions_bedfile
    output:
        tuple val(sample_id), path('*.filtered.fastq'), val(approx_size), optional: true, emit: unmapped
        path '*.stats', optional: true, emit: host_filter_stats
        tuple val(sample_id), env(STATUS), emit: status
    script:
        def name = sample_id
        def regs = regions_bedfile.name != 'NO_REG_BED' ? regions_bedfile : 'none'
    """
    STATUS="Failed due to filtered host reads"
    (minimap2 -t $task.cpus -y -ax map-ont $reference $fastq \
        | samtools sort -o ${name}.sorted.aligned.bam -
    samtools index ${name}.sorted.aligned.bam
    samtools view -b -f 4  ${name}.sorted.aligned.bam > unmapped.bam
    samtools view -b -F 4  ${name}.sorted.aligned.bam > mapped.bam
    samtools fastq unmapped.bam > ${name}.filtered.fastq
    fastcat -s ${name} -r ${name}.interim ${name}.filtered.fastq > /dev/null
    if [[ "\$(wc -l <"${name}.interim")" -ge "1" ]];  then
        mv ${name}.interim ${name}.stats
    fi
    if [[ -f "$regs" ]]; then
        bedtools intersect -a mapped.bam -b $regs -wa \
            | samtools view -bh - > retained.bam
        samtools fastq retained.bam >> ${name}.filtered.fastq
    fi ) && STATUS="Completed successfully"
    """
}

process lookup_medaka_model {
    label 'wfplasmid'
    cpus 1
    memory '1GB'
    input:
        path('lookup_table')
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
    label 'medaka'
    cpus params.threads
    memory '4GB'
    input:
        tuple val(sample_id), path(draft), path(fastq), val(medaka_model)
    output:
        tuple val(sample_id), path('*.final.fasta'), emit: polished
        tuple val(sample_id), env(STATUS), emit: status
        tuple val(sample_id), path("${sample_id}.final.fastq"), emit: assembly_qc
    script:
        def model = medaka_model

    """
    STATUS="Failed to polish assembly with Medaka"
    medaka_consensus -i "${fastq}" -d "${draft}" -m "${model}" -o . -t $task.cpus -f -q
    echo ">${sample_id}" >> "${sample_id}.final.fasta"
    sed "2q;d" consensus.fastq >> "${sample_id}.final.fasta"
    mv consensus.fastq "${sample_id}.final.fastq"
    STATUS="Completed successfully"
    """
}

process downsampledStats {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        tuple val(sample_id), path(sample)
    output:
        path '*.stats', optional: true
    """
    fastcat -s ${sample_id} -r ${sample_id}.downsampled $sample > /dev/null
    if [[ "\$(wc -l <"${sample_id}.downsampled")" -ge "2" ]];  then
        mv ${sample_id}.downsampled ${sample_id}.stats
    fi
    """
}

process findPrimers {
    errorStrategy 'ignore'
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        path primers
        tuple val(sample_id), path(sequence)
    output:
        path '*.bed', optional: true
    shell:
    '''
    cat !{sequence} | seqkit amplicon -p !{primers} -m 3 -j !{task.cpus} --bed >> !{sample_id}.interim
    if [[ "$(wc -l <"!{sample_id}.interim")" -ge "1" ]];  then
        mv !{sample_id}.interim !{sample_id}.bed
    fi
    '''
}

process medakaVersion {
    label 'medaka'
    cpus 1
    memory '2GB'
    output:
        path 'medaka_version.txt'
    """
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}

process flyeVersion {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        path 'versions.txt'
    output:
        path 'assembly_version.txt'
    """
    cat "versions.txt" >> "assembly_version.txt"
    flye --version |  sed 's/^/flye,/' >> "assembly_version.txt"
    """
}

process canuVersion {
    label 'canu'
    cpus 1
    memory '2GB'
    input:
        path 'versions.txt'
    output:
        path 'assembly_version.txt'
    """
    cat "versions.txt" >> "assembly_version.txt"
    canu -version | sed 's/ /,/' >>"assembly_version.txt"
    """
}

process getVersions {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        path 'input_versions.txt'
    output:
        path 'versions.txt'
    script:
    """
    cat "input_versions.txt" >> "versions.txt"
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    trycycler --version | sed 's/ /,/' >> versions.txt
    bedtools --version | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    rasusa --version | sed 's/ /,/' >> versions.txt
    python -c "import spoa; print(spoa.__version__)" | sed 's/^/spoa,/'  >> versions.txt
    python -c "import pandas; print(pandas.__version__)" | sed 's/^/pandas,/'  >> versions.txt
    python -c "import plannotate; print(plannotate.__version__)" | sed 's/^/plannotate,/'  >> versions.txt
    python -c "import bokeh; print(bokeh.__version__)" | sed 's/^/bokeh,/'  >> versions.txt
    """
}

process getParams {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    output:
        path 'params.json'
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process assemblyMafs {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        tuple val(sample_id), path('assembly.fasta')
    output:
        tuple val(sample_id), path("${sample_id}.assembly.maf"), emit: assembly_maf
    // set -m(multiplicity) to 10000 to increase sensitivity from default of 10
    // for assemblies this small computational cost is low
    // reduce offset distance for suppressing repeats inside exact matches -w
    // from default of 1000 to 10.
    """
    lastdb db.lastdb "assembly.fasta"
    lastal -m 10000 -w 10 db.lastdb "assembly.fasta" > "${sample_id}.assembly.maf"

    """
}

process runPlannotate {
    label 'wfplasmid'
    cpus 1
    memory params.large_construct ? '8GB' : '2GB'
    input:
        path annotation_database
        path 'assemblies/*'
        path final_status
    output:
        path 'feature_table.txt', emit: feature_table
        path 'plannotate.json', emit: json
        path '*annotations.bed', optional: true, emit: annotations
        path 'plannotate_report.json', emit: report
        path '*annotations.gbk', optional: true, emit: gbk
    script:
        def database =  annotation_database.name.startsWith('OPTIONAL_FILE') ? 'Default' : "${annotation_database}"
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
    label 'wfplasmid'
    cpus 1
    memory '1GB'
    input:
    path 'primer_beds/*'
    path 'assemblies/*'
    path align_ref
    output:
        path 'inserts/*', optional: true, emit: inserts
        path '*.json', emit: json
    script:
        def ref =  align_ref.name.startsWith('OPTIONAL_FILE') ? '' : "--reference ${align_ref}"
        def large_construct = params.large_construct ? '--large_construct' : ''
    """
    if [ -e "primer_beds/OPTIONAL_FILE" ]; then
        inserts=""
    else
        inserts="--primer_beds primer_beds/*"
    fi
    workflow-glue find_inserts \$inserts $ref $large_construct
    """
}

process insert_qc {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
    tuple val(sample_id), path('insert_assembly.fasta')
    path 'reference_assembly.fasta'
    output:
    tuple val(sample_id), path("${sample_id}.calls.bcf"), path("${sample_id}.stats"), optional: true, emit: insert_stats
    tuple val(sample_id), env(STATUS), emit: status
    script:
    """
    STATUS="Insert found but does not align with provided reference"
    minimap2 -t $task.cpus -y -ax map-ont "reference_assembly.fasta" "insert_assembly.fasta" | samtools sort -o output.bam -
    mapped=\$(samtools view -F 4 -c  output.bam)
    if [ \$mapped != 0 ]; then
        bcftools mpileup -Ou -f "reference_assembly.fasta" output.bam | bcftools call -mv -Ob -o "${sample_id}.calls.bcf"
        bcftools stats "${sample_id}.calls.bcf" > ${sample_id}.stats
        STATUS="Completed successfully"
    fi
    """
}

process assembly_qc {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        tuple val(sample_id), path('assembly.fastq')
    output:
        path "${sample_id}.assembly_stats.tsv"
    script:
    """
    fastcat -s "${sample_id}" -r "${sample_id}.assembly_stats.tsv" assembly.fastq
    """
}

// downsampled, per barcode and host filtered stats files are handled earlier in the workflow and need to be named with the sample alias
process report {
    label 'wfplasmid'
    cpus 1
    memory '2GB'
    input:
        path 'downsampled_stats/*'
        path final_status
        path 'per_barcode_stats/*'
        path 'host_filter_stats/*'
        path 'versions/*'
        path 'params.json'
        path plannotate_json
        path inserts_json
        path lengths
        path 'qc_inserts/*'
        path 'assembly_quality/*'
        path 'mafs/*'
        path client_fields
    output:
        path 'wf-clone-validation-*.html', emit: html
        path 'sample_status.txt', emit: sample_stat
        path 'inserts/*', optional: true, emit: inserts
    script:
        report_name = 'wf-clone-validation-report.html'
        String client_fields_args = client_fields.name != 'OPTIONAL_FILE' ? "--client_fields ${client_fields}" : ''

    """
    workflow-glue report \
     $report_name \
    --downsampled_stats downsampled_stats/* \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --status $final_status \
    --per_barcode_stats per_barcode_stats/* \
    --host_filter_stats host_filter_stats/* \
    --params params.json \
    --versions versions \
    --plannotate_json $plannotate_json \
    --lengths $lengths \
    --inserts_json $inserts_json \
    --qc_inserts qc_inserts \
    --assembly_quality assembly_quality/* \
    --mafs mafs \
    --assembly_tool ${params.assembly_tool} \
    $client_fields_args
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
        min_read_length
        max_read_length

    main:
        // remove samples that didn't have sequences (i.e. metamap entries without
        // corresponding barcode sub-directories) and get the per-read stats file from
        // the fastcat stats dir
        samples = samples
        | filter { it[1] }
        | map { it[2] = it[2].resolve('per-read-stats.tsv.gz'); it }

        // Min/max filter reads
        fastcat_extra_args = "-a $min_read_length -b $max_read_length"

        // drop samples with too low coverage
        sample_fastqs = checkIfEnoughReads(samples, fastcat_extra_args)

        // Optionally filter the data, removing reads mapping to
        // the host or background genome
        if (host_reference.name != 'NO_HOST_REF') {
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
        named_drafts_samples = named_drafts
        .join(named_samples, failOnMismatch: false, remainder: true)
        .filter { alias, assembly, fastq -> (assembly != null) }

        if (params.medaka_model) {
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

        downsampled_stats = downsampledStats(assemblies.downsampled)

        primer_beds = findPrimers(primers, polished.polished)
        medaka_version = medakaVersion()
        if (params.assembly_tool == 'flye') {
            assembly_version = flyeVersion(medaka_version)
        }else {
            assembly_version = canuVersion(medaka_version)
        }
        software_versions = getVersions(assembly_version)
        workflow_params = getParams()

        insert = inserts(primer_beds.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            align_ref)

        // assembly QC

        insert_tuple = insert.inserts.flatten().map { assembly -> tuple(assembly.simpleName, assembly) }
        assembly_quality = assembly_qc(polished.assembly_qc)
        if (params.insert_reference) {
            insert_qc_tuple = insert_qc(insert_tuple, align_ref)
            qc_insert = insert_qc_tuple.insert_stats.map { sample_id, bcf, stats -> stats }
            bcf_insert = insert_qc_tuple.insert_stats.map { sample_id, bcf, stats -> bcf }
            insert_status = insert_qc_tuple.status
        }
        else {
            qc_insert = Channel.empty()
            bcf_insert = Channel.empty()
            insert_status = Channel.empty()
        }

        // Concat statuses and keep the last of each
        final_status = sample_fastqs.status.concat(updated_status)
        .concat(assemblies.status).concat(polished.status).concat(insert_status).groupTuple()
        .map { it -> it[0].toString() + ',' + it[1][-1].toString() }
        final_status = final_status.collectFile(name: 'final_status.csv', newLine: true)

        annotation = runPlannotate(
            database, polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status)

        mafs = assemblyMafs(polished.polished)

        client_fields = params.client_fields && file(params.client_fields).exists() ? file(params.client_fields) : file("$projectDir/data/OPTIONAL_FILE")

        report = report(
            downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status,
            sample_fastqs.stats.collect(),
            filtered_stats,
            software_versions.collect(),
            workflow_params,
            annotation.report,
            insert.json,
            annotation.json,
            qc_insert.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            assembly_quality.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            mafs.map { meta, maf -> maf }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            client_fields
            )

        results = polished.polished.map { meta, polished -> polished }.concat(
            report.html,
            report.sample_stat,
            annotation.feature_table,
            insert.inserts,
            annotation.json,
            annotation.annotations,
            annotation.gbk,
            workflow_params,
            bcf_insert,
            qc_insert,
            polished.assembly_qc.map { meta, assembly_qc -> assembly_qc })

    emit:
        results
        telemetry = workflow_params
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish {
    // publish inputs to output directory
    label 'wfplasmid'
    publishDir "${params.out_dir}", mode: 'copy', pattern: '*', saveAs: {
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        file fname
    output:
        file fname
    '''
    echo "Writing output files"
    '''
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.containsKey('reference')) {
        throw new Exception('--reference is deprecated, use --insert_reference instead.')
    }

    // calculate min and max read length for filtering with `fastcat`
    int min_read_length, max_read_length
    if (params.sample_sheet) {
        sample_sheet_csv = file(params.sample_sheet).splitCsv(header:true)
    }
    if (params.sample_sheet && sample_sheet_csv[0].containsKey('approx_size')) {
        approx_sizes = sample_sheet_csv['approx_size']
        if (approx_sizes.contains(null)) {
            throw new Exception('Either use the `--approx_size` parameter or include an `approx_size` column in the sample sheet for all samples.')
        }
        log.warn 'Overriding the approx size parameter with per sample approx sizes provided by the sample_sheet.'
        // the file provided with `--sample_sheet` contains a size estimate for
        // each sample, but we will filter all samples with the same parameters)
        min_read_length = approx_sizes.collect { it.toInteger() }.min()
        max_read_length = approx_sizes.collect { it.toInteger() }.max()
    } else {
        // we only got a single size estimate --> take as max and min
        min_read_length = max_read_length = params.approx_size
    }
    // +/- 50% margins for read length thresholds
    min_read_length *= 0.5
    // if canu is requested, log a warning
    if (params.assembly_tool == 'canu') {
        log.warn 'Assembly tool Canu. This may result in suboptimal performance on ARM devices.'
    }
    // if large construct don't filter out shorter reads as approx size no longer equal to read length
    if (params.large_construct) {
        min_read_length = 200
    }
    max_read_length *= 1.5

    samples = fastq_ingress([
        'input':params.fastq,
        'sample':params.sample,
        'sample_sheet':params.sample_sheet,
        'stats': true
        ])

    // add the size estimates to the channel with the samples
    // by joining the samples on the alias key with approx size
    if (params.sample_sheet && sample_sheet_csv[0].containsKey('approx_size')) {
        sample_alias = samples.map { [it[0]['alias'], *it] }
        approx_size_alias = Channel.of(*sample_sheet_csv).map { [it['alias'], it['approx_size']] }
        samples = sample_alias.join(approx_size_alias).map {  it[1..-1]  }
    } else {
        samples = samples.map { [*it, params.approx_size] }
    }

    host_reference = params.host_reference ?: 'NO_HOST_REF'
    host_reference = file(host_reference, checkIfExists: host_reference == 'NO_HOST_REF' ? false : true)
    regions_bedfile = params.regions_bedfile ?: 'NO_REG_BED'
    regions_bedfile = file(regions_bedfile, checkIfExists: regions_bedfile == 'NO_REG_BED' ? false : true)
    primer_file = file("$projectDir/data/OPTIONAL_FILE")
    if (params.primers != null) {
        primer_file = file(params.primers, type: 'file')
    }
    align_ref = file("$projectDir/data/OPTIONAL_FILE")
    if (params.insert_reference != null) {
        align_ref = file(params.insert_reference, type: 'file')
    }
    database = file("$projectDir/data/OPTIONAL_FILE")
    if (params.db_directory != null) {
        database = file(params.db_directory, type: 'dir')
    }

    // Run pipeline
    results = pipeline(
        samples,
        host_reference,
        regions_bedfile,
        database,
        primer_file,
        align_ref,
        min_read_length,
        max_read_length)

    publish(results[0])
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
