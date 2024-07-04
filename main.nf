#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
if (params.assembly_tool == 'canu'){
    include {assembleCore_canu as assembleCore} from "./modules/local/canu_assembly.nf"
} else {
    include {assembleCore_flye as assembleCore} from "./modules/local/flye_assembly.nf"
}

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process checkIfEnoughReads {
    label "wfplasmid"
    cpus params.threads
    memory "2GB"
    input:
        tuple val(meta), path("input.fastq.gz"), path("fastcat_stats")
        val min_read_length
        val max_read_length
    output:
        tuple val(meta), path("${meta.alias}.fastq.gz"),
            optional: true, emit: sample
        path("${meta.alias}.fastcat_stats"), emit: stats
        tuple val(meta.alias), env(STATUS), emit: status
    script:
        def extra_args = "-a $min_read_length -b $max_read_length"
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        
        int value = (expected_depth.toInteger()) * 0.8
        int bgzip_threads = task.cpus == 1 ? 1 : task.cpus - 1
    """
    STATUS="Failed due to insufficient reads"
    mv fastcat_stats "${meta.alias}.fastcat_stats"
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
    label "wfplasmid"
    cpus params.threads
    memory "4GB"
    input:
        tuple val(meta), path(fastq)
        path reference
        path regions_bedfile
    output:
        tuple val(meta), path("*.filtered.fastq"), optional: true, emit: unmapped
        path "*.stats", optional: true, emit: host_filter_stats
        tuple val(meta.alias), env(STATUS), emit: status
    script:
        def name = meta.alias
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


process medakaPolishAssembly {
    label "medaka"
    cpus params.threads
    memory "4GB"
    input:
        tuple val(meta), path(draft), path(fastq)
    output:
        tuple val(meta), path("*.final.fasta"), emit: polished
        tuple val(meta.alias), env(STATUS), emit: status
        tuple val(meta), path("${meta.alias}.final.fastq"), emit: assembly_qc
    script:
    String model_args = params.basecaller_cfg ? "-m ${params.basecaller_cfg}:consensus": ""
    """
    STATUS="Failed to polish assembly with Medaka"
    medaka_consensus -i "${fastq}" -d "${draft}" ${model_args} -o . -t $task.cpus -f -q
    echo ">${meta.alias}" >> "${meta.alias}.final.fasta"
    sed "2q;d" consensus.fastq >> "${meta.alias}.final.fasta"
    mv consensus.fastq "${meta.alias}.final.fastq"
    STATUS="Completed successfully"
    """
}


process downsampledStats {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path(sample)
    output:
        path "*.stats", optional: true
    """
    fastcat -s ${meta.alias} -r ${meta.alias}.downsampled $sample > /dev/null
    if [[ "\$(wc -l <"${meta.alias}.downsampled")" -ge "2" ]];  then
        mv ${meta.alias}.downsampled ${meta.alias}.stats
    fi
    """
}


process findPrimers {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
        path primers
        tuple val(meta), path(sequence)
    output:
        tuple val(meta), path(sequence), path("*.bed"), optional: true
    shell:
    '''
    cat !{sequence} | seqkit amplicon -p !{primers} -m 3 -j !{task.cpus} --bed >> !{meta.alias}.interim
    if [[ "$(wc -l <"!{meta.alias}.interim")" -ge "1" ]];  then
        mv !{meta.alias}.interim !{meta.alias}.bed
    fi
    '''
}

process medakaVersion {
    label "medaka"
    cpus 1 
    memory "2GB"
    output:
        path "medaka_version.txt"
    """
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}

process flyeVersion {
    label "wfplasmid"
    cpus 1 
    memory "2GB"
    input:
        path "versions.txt"
    output:
        path "assembly_version.txt"
    """
    cat "versions.txt" >> "assembly_version.txt"
    flye --version |  sed 's/^/flye,/' >> "assembly_version.txt"
    """
}

process canuVersion {
    label "canu"
    cpus 1 
    memory "2GB"
    input:
        path "versions.txt"
    output:
        path "assembly_version.txt"
    """
    cat "versions.txt" >> "assembly_version.txt"
    canu -version | sed 's/ /,/' >>"assembly_version.txt"
    """
}

process getVersions {
    label "wfplasmid"
    cpus 1
    memory "2GB"
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
    label "wfplasmid"
    cpus 1
    memory "2GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process assemblyMafs {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path("assembly.fasta")
    output:
        tuple val(meta), path("${meta.alias}.assembly.maf"), emit: assembly_maf
    // set -m(multiplicity) to 10000 to increase sensitivity from default of 10
    // for assemblies this small computational cost is low
    // reduce offset distance for suppressing repeats inside exact matches -w 
    // from default of 1000 to 10.
    """
    lastdb db.lastdb "assembly.fasta"
    lastal -m 10000 -w 10 db.lastdb "assembly.fasta" > "${meta.alias}.assembly.maf"
  
    """
}

// this process uses plannotate version v1.2.2 which fixes bugs when finding features from the database
// bokeh is pinned to version 2
process runPlannotate {
    label "plannotate"
    cpus 1
    memory params.large_construct ? "8GB" : "2GB"
    input:
        path annotation_database
        path "assemblies/*"
        path final_status
    output:
        path "feature_table.txt", emit: feature_table
        path "plannotate.json", emit: json
        path "*annotations.bed", optional: true, emit: annotations
        path "plannotate_report.json", emit: report
        path "*annotations.gbk", optional: true, emit: gbk
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
    memory "1GB"
    input:
         tuple path(align_ref), path("assemblies/*"), path("primer_beds/*")
    output:
        path "inserts/*", optional: true, emit: inserts
        path "insert_data_${task.index}.json", emit: json
    script:
        def output = "insert_data_${task.index}.json"
        def ref =  align_ref.name.startsWith('OPTIONAL_FILE') ? '' : "--reference ${align_ref}"
        def large_construct = params.large_construct ? "--large_construct" : ""
    // As primer_beds may have *.bed and OPTIONAL_FILE in there we need a second check
    // that OPTIONAL_FILE is the only file there
    """
    if [ -e "primer_beds/OPTIONAL_FILE" ] && [ \$(ls -1 primer_beds | wc -l ) -eq 1 ]; then
        inserts=""
    else
        inserts="--primer_beds primer_beds/*.bed"
    fi
    workflow-glue find_inserts --output $output \$inserts $ref $large_construct
    """
}


process insert_qc {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
         tuple val(meta), path("insert_assembly.fasta"), path("insert_ref.fasta")
    output:
         tuple val(meta), path("${meta.alias}.insert.calls.bcf"), path("${meta.alias}.insert.stats"), optional: true, emit: insert_stats
         tuple val(meta.alias), env(STATUS), emit: status
    script:
    """
    STATUS="Insert found but does not align with provided reference"
    minimap2 -t $task.cpus -y -ax map-ont insert_ref.fasta "insert_assembly.fasta" | samtools sort -o output.bam -
    mapped=\$(samtools view -F 4 -c  output.bam)
    if [ \$mapped != 0 ]; then
        # -m1 reduces evidence for indels to 1 read instead of the default 2
        bcftools mpileup -m1 -Ou -f insert_ref.fasta output.bam | bcftools call -mv -Ob -o "${meta.alias}.insert.calls.bcf"
        bcftools stats "${meta.alias}.insert.calls.bcf" > "${meta.alias}.insert.stats"
        STATUS="Completed successfully"
    fi
    """
}


process align_assembly {
    label "wfplasmid"
    cpus 3
    memory "8GB"
    input:
        tuple val(meta), path("full_assembly.fasta"), path("full_reference.fasta")
    output:
        tuple val(meta), path("${meta.alias}.bam"), path("${meta.alias}.bam.bai"), path("full_reference.fasta"), optional: true
    script:
    // Align assembly with the reference, which results in 2 aligned segments a (primary and supplementary alignment)
    // Create a consensus to get the two sections of the assembled sequence in one Fasta and the same order as the reference
    // Align this with the reference to get final alignment to output to user
    """
    minimap2 -t ${task.cpus - 1} -y -ax asm5 full_reference.fasta "full_assembly.fasta" | samtools sort -o "initial_alignment.bam" -
    mapped=\$(samtools view -F 4 -c  "initial_alignment.bam")
    if [ \$mapped != 0 ]; then
        samtools consensus --mode simple "initial_alignment.bam" -o "consensus.bam"
        samtools fasta consensus.bam > consensus.fasta
        minimap2 -t ${task.cpus - 1} -y -ax asm5 full_reference.fasta "consensus.fasta" \
            | samtools sort -@ ${task.cpus} --write-index -o ${meta.alias}.bam##idx##${meta.alias}.bam.bai -
    fi
    """
}


process assembly_comparison {
    label "wfplasmid"
    cpus 2
    memory "2GB"
    input:
        tuple val(meta), path("${meta.alias}.bam"), path("${meta.alias}.bam.bai"), path("full_reference.fasta")
    output: 
        path("${meta.alias}.bam.stats"), emit: bamstats
        tuple val(meta), path("${meta.alias}.full_construct.calls.bcf"), path("${meta.alias}.full_construct.stats"), emit: full_assembly_stats
    script:
    // use bamstats to get percent identity and reference coverage
    // Also get variants report
    """
    bamstats -t ${task.cpus} -s "${meta.alias}" "${meta.alias}.bam" > "${meta.alias}.bam.stats"
    # -m1 reduces evidence for indels to 1 read instead of the default 2
    bcftools mpileup -m1 -Ou -f full_reference.fasta "${meta.alias}.bam" | bcftools call -mv -Ob -o "${meta.alias}.full_construct.calls.bcf"
    bcftools stats "${meta.alias}.full_construct.calls.bcf" > ${meta.alias}.full_construct.stats
    """
}



process assembly_qc {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path("assembly.fastq")
    output:
        path "${meta.alias}.assembly_stats.tsv"
    script:
    """
    fastcat -s "${meta.alias}" -r "${meta.alias}.assembly_stats.tsv" assembly.fastq
    """
}


// to predict linearisation efficiency
// find cut site location in reference
// align reads with provided reference
// find intersection of alignments and cutsite
// if reads overlap cut site -> circular, not linearized
// linearisation efficiency = total reads - intersection (Calculated in report)
process cutsite_qc {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path("reads.fastq.gz"), path("full_reference.fasta")
    output:
        tuple val(meta.alias), val(meta.n_seqs), env(CUT_COUNT), emit: cut_counts
    script:
    """
    cat full_reference.fasta | seqkit locate -p ${meta.cut_site} -m ${params.cutsite_mismatch} --bed  --circular > locate.bed
    if [[ \$(<locate.bed wc -l) -ne 1 ]]; then
        echo "Found unexpected number of cut sites in reference. Check there is only one cut site within the reference."
        exit 1
    fi
    minimap2 -y -ax map-ont full_reference.fasta "reads.fastq.gz" | samtools view -F 2304  -Sb - | samtools sort -O bam - | bedtools bamtobed -i - > aligned.bed
    bedtools intersect -a aligned.bed -b locate.bed > intersection.txt
    CUT_COUNT=\$(<intersection.txt wc -l)
    """
}

// downsampled, per barcode and host filtered stats files are handled earlier in the workflow and need to be named with the sample alias
// uses plannotate version v1.2.0 where bokeh is not pinned and therefore compatible with 
// bokeh version 3 which is required for ezcharts.
process report {
    label "wfplasmid"
    cpus 1
    memory "4GB"
    input:
        val metadata
        path "downsampled_stats/*"
        path final_status
        tuple path("per_barcode_stats/*"), val(no_stats)
        path "host_filter_stats/*"
        path "versions/*"
        path "params.json"
        path plannotate_json
        path ("inserts_json/*")
        path lengths
        path "qc_inserts/*"
        path "assembly_quality/*"
        path "mafs/*"
        path client_fields
        val wf_version
        path "full_assembly_variants/*"
        path "assembly_bamstats/*"
        path cutsite_csv, stageAs: "cutsite_csv/*"
    output:
        path "wf-clone-validation-*.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "inserts/*", optional: true, emit: inserts
    script:
        report_name = "wf-clone-validation-report.html"
        String stats_args = no_stats ? "" : "--per_barcode_stats per_barcode_stats/*"
        def metadata_obj = new JsonBuilder(metadata)
        String metadata = metadata_obj.toPrettyString()
        String client_fields_args = client_fields.name != "OPTIONAL_FILE" ? "--client_fields ${client_fields}" : ""
        String cutsite = cutsite_csv.fileName.name == "OPTIONAL_FILE" ? "" : "--cutsite_csv cutsite_csv/*"
    """
    echo '${metadata}' > metadata.json
    if [ -f "full_assembly_variants/OPTIONAL_FILE" ]; then
        full_assembly_variants=""
    else
        full_assembly_variants="--full_assembly_variants full_assembly_variants"
    fi
    if [ -f "assembly_bamstats/OPTIONAL_FILE" ]; then
        assembly_bamstats=""
    else
        assembly_bamstats="--reference_alignment_bamstats assembly_bamstats/*"
    fi
    workflow-glue report \
     $report_name \
    --metadata metadata.json \
    --downsampled_stats downsampled_stats/* \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --status $final_status \
    $stats_args \
    --host_filter_stats host_filter_stats/* \
    --params params.json \
    --versions versions \
    --plannotate_json $plannotate_json \
    --lengths $lengths \
    --inserts_json inserts_json/* \
    --qc_inserts qc_inserts \
    --assembly_quality assembly_quality/* \
    --mafs mafs \
    --assembly_tool ${params.assembly_tool} \
    $client_fields_args \
    --wf_version $wf_version \
    \$assembly_bamstats \
    \$full_assembly_variants \
    $cutsite
    """
}

workflow pipeline {
    take:
        samples
        host_reference
        regions_bedfile
        database
        primers
        min_read_length
        max_read_length
        cutsite_csv // alias, read_counts, cutsite_count

    main:
        // Min/max filter reads
        
    
        // drop samples with too low coverage
        sample_fastqs = checkIfEnoughReads(samples, min_read_length, max_read_length)

        // Optionally filter the data, removing reads mapping to
        // the host or background genome
        if (host_reference.name != "NO_HOST_REF") {
            filtered = filterHostReads(
                    sample_fastqs.sample, host_reference, regions_bedfile)
            samples_filtered = filtered.unmapped
            updated_status = filtered.status
            filtered_stats = filtered.host_filter_stats.collect()
                             .ifEmpty(OPTIONAL_FILE)
        }
        else {
            samples_filtered = sample_fastqs.sample
            updated_status = sample_fastqs.status
            filtered_stats = OPTIONAL_FILE
        }

        // Core assembly and reconciliation
        assemblies = assembleCore(samples_filtered)
        named_drafts = assemblies.assembly.groupTuple()
        named_samples = assemblies.downsampled.groupTuple()
        named_drafts_samples = named_drafts
        .join(named_samples, failOnMismatch: false, remainder: true)
        .filter{alias, assembly, fastq -> (assembly != null)}

        // Polish draft assembly
        polished = medakaPolishAssembly(named_drafts_samples)
    
        downsampled_stats = downsampledStats(assemblies.downsampled)

        primer_beds = findPrimers(primers, polished.polished)
        medaka_version = medakaVersion()
        if (params.assembly_tool == "flye"){
            assembly_version = flyeVersion(medaka_version)
        }else{
            assembly_version = canuVersion(medaka_version)
        }
        software_versions = getVersions(assembly_version)
        workflow_params = getParams()




        // Group samples by insert_reference for insert analysis
        // the resulting output should be [ insert_reference, [ assemblies ], [ primers ]]
        ref_groups = primer_beds
        | map { meta, assembly, primers -> [ meta.insert_reference ? meta.insert_reference : OPTIONAL_FILE, assembly, primers ] }
        | groupTuple()
        | ifEmpty([OPTIONAL_FILE, OPTIONAL_FILE, OPTIONAL_FILE]) // to match cardinatlity expected by process insert()


        insert = inserts(ref_groups)
        
        // assembly QC
        // Insert
        assembly_quality = assembly_qc(polished.assembly_qc)
        insert_qc_tuple = insert_qc(polished.polished
        | filter { meta, assembly -> meta.insert_reference}
        | map{meta, assembly -> [meta, assembly, meta.insert_reference]})

        qc_insert = insert_qc_tuple.insert_stats.map{meta, bcf, stats -> stats}
        bcf_insert = insert_qc_tuple.insert_stats.map{meta, bcf, stats -> bcf}
        insert_status = insert_qc_tuple.status
        
        // Full reference
        qc_full = align_assembly( polished.polished
        | filter {meta, assembly -> meta.full_reference}
        | map {meta, assembly -> [meta, assembly, meta.full_reference ]})
        assembly_comparison_output = assembly_comparison( qc_full )
        ref_bamstats = assembly_comparison_output.bamstats
        full_assembly_stats = assembly_comparison_output.full_assembly_stats.map{meta, bcf, stats -> stats}
        bcf = assembly_comparison_output.full_assembly_stats.map{meta, bcf, stats -> bcf}
        bam = qc_full.map{meta, bam, bai, ref -> [bam, bai]}


        // Concat statuses and keep the last of each
        final_status = sample_fastqs.status.concat(updated_status)
        .concat(assemblies.status).concat(polished.status).concat(insert_status).groupTuple()
        .map { it -> it[0].toString() + ',' + it[1][-1].toString() }
        final_status = final_status.collectFile(name: 'final_status.csv', newLine: true)

        annotation = runPlannotate(
            database, polished.polished.map { it -> it[1] }.collect().ifEmpty(OPTIONAL_FILE),
            final_status)

        mafs = assemblyMafs(polished.polished)

        client_fields = params.client_fields && file(params.client_fields).exists() ? file(params.client_fields) : OPTIONAL_FILE
        stats = sample_fastqs.stats | ifEmpty(OPTIONAL_FILE) | collect | map { [it, it[0] == OPTIONAL_FILE] }
        // Can't collect whole meta likely due to full_ref and insert_ref
        // Causes java.lang.StackOverflowError in JsonBuilder so keeping only those used in report
        meta = samples.map{ meta, reads, stats -> ["alias": meta.alias, "n_seqs": meta.n_seqs] }.collect()

        report = report(
            meta,
            downsampled_stats.collect().ifEmpty(OPTIONAL_FILE),
            final_status,
            stats,
            filtered_stats,
            software_versions.collect(),
            workflow_params,
            annotation.report,
            insert.json.collect(),
            annotation.json,
            qc_insert.collect().ifEmpty(OPTIONAL_FILE),
            assembly_quality.collect().ifEmpty(OPTIONAL_FILE),
            mafs.map{ meta, maf -> maf}.collect().ifEmpty(OPTIONAL_FILE),
            client_fields,
            workflow.manifest.version,
            full_assembly_stats.collect().ifEmpty(OPTIONAL_FILE),
            ref_bamstats.collect().ifEmpty(OPTIONAL_FILE),
            cutsite_csv
            )

        results = polished.polished.map { meta, polished -> polished }.concat(
            report.html,
            report.sample_stat,
            annotation.feature_table,
            insert.inserts.collect(),
            annotation.json,
            annotation.annotations,
            annotation.gbk,
            workflow_params,
            bcf_insert,
            qc_insert,
            polished.assembly_qc.map { meta, assembly_qc -> assembly_qc },
            bam.flatten(),
            full_assembly_stats,
            bcf,
            ref_bamstats)
        
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
    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.containsKey("reference")) {
        throw new Exception("--reference is deprecated, use --insert_reference instead.")
    }

   
    if (params.assembly_tool == 'canu'){
        log.warn "Assembly tool Canu. This may result in suboptimal performance on ARM devices."
    }

     Map ingress_args = [
        "sample": params.sample,
        "sample_sheet": params.sample_sheet,
        "analyse_unclassified": params.analyse_unclassified,
        "stats": true,
    ]


    // get input data
    if (params.fastq) {
        samples = fastq_ingress(ingress_args + ["input": params.fastq])
    } else {
        samples = xam_ingress(
            ingress_args + ["input": params.bam, "keep_unaligned": true, "return_fastq": true]
        )
    }

    samples = samples
    | map{ meta, read, stats -> 
        // Check sample_sheet data in metamap
        // If column present in sample sheet but missing item in row value will be "" 
        // Full refernce checks
        if (meta.full_reference && params.full_reference){
            throw new Exception("Parameter --full_reference cannot be used in conjunction with a sample sheet with column 'full_reference'")
        }
        if (meta.full_reference == ""){
            throw new Exception("When using a sample sheet with 'full_reference' column please provide a reference for each sample.")
        }

        // Insert reference checks
        if (meta.insert_reference && params.insert_reference){
            throw new Exception("Parameter --full_reference cannot be used in conjunction with a sample sheet with column 'full_reference'")
        }
        if (meta.insert_reference == ""){
            throw new Exception("When using a sample sheet with 'full_reference' column please provide a reference for each sample.")
        }
        // Approx size checks
        if (meta.approx_size && params.approx_size){
            log.warn("Parameter --approx_size will be overwritten by column 'approx_size'")
        }
        if (meta.approx_size == ""){
            throw new Exception("When using a sample sheet with `approx_size` column please provide  for all samples.")
        }

        // Cut site
        if (meta.cut_site == ""){
            throw new Exception("One or more cut_sites from the sample sheet were missing. Add one cut site per sample.")
        }
        if (meta.cut_site){
            if (!(params.full_reference || meta.full_reference )){
                throw new Exception("As cut site was supplied in sample sheet the full reference must be provided.")
            }
            if (!meta.cut_site.matches('^[ACGT]+$')){
                throw new Exception("Cut site column must not contain base pairs other than AGCT.")
            }
            if (meta.cut_site.size() < 5 || meta.cut_site.size() > 30 ){
                throw new Exception("Cut sites must all be between 5-30 bp long.")
            }
        }
        // Check that files in sample_sheets exist
        // Turn file objects toString() to keep in meta and stop StackOverflowError
	    full_reference = (meta.full_reference ?: params.full_reference) 
        insert_reference = (meta.insert_reference ?: params.insert_reference)
        approx_size = (meta.approx_size ?: params.approx_size) as int
        // Check if null before adding file around it ()
        new_keys = ["full_reference": full_reference ? file(full_reference, checkIfExists: true) : null,
            "insert_reference": insert_reference ? file(insert_reference, checkIfExists: true) : null,
            "approx_size": approx_size
        ]
        return [ meta + new_keys, read, stats ]
    }

    // if large construct don't filter out shorter reads as approx size no longer equal to read length
    if (params.large_construct){
        min_read_length = 200
    } else {
        min_read_length = samples
        | map{ it[0].approx_size }
        | min()
        | map{ it * 0.5 }

    }
    max_read_length = samples
        | map{ it -> it[0].approx_size }
        | max
        | map{ it * 1.5}
    

    host_reference = params.host_reference ?: 'NO_HOST_REF'
    host_reference = file(host_reference, checkIfExists: host_reference == 'NO_HOST_REF' ? false : true)
    regions_bedfile = params.regions_bedfile ?: 'NO_REG_BED'
    regions_bedfile = file(regions_bedfile, checkIfExists: regions_bedfile == 'NO_REG_BED' ? false : true)
    primer_file = OPTIONAL_FILE
    if (params.primers != null){
        primer_file = file(params.primers, type: "file")
    }


    database = OPTIONAL_FILE
    if (params.db_directory != null){
         database = file(params.db_directory, type: "dir")

    }

    // remove samples that didn't have sequences (i.e. metamap entries without
    // corresponding barcode sub-directories)
    samples = samples | filter { it[1] }

    // If cut site in sample sheet
    // Create cutsite csv from read count in meta and cut site alignment counts

    cut = cutsite_qc(samples
    | filter{ meta, reads, stats -> meta.cut_site}
    | map {meta, reads, stats -> [meta, reads, meta.full_reference] }
    )
    cutsite_csv = cut
    | map { alias, read_count, cut_site_count -> "$alias,$read_count,$cut_site_count" }
    | collectFile(name: "cut_sites.csv", newLine: true)
    | ifEmpty(OPTIONAL_FILE)

    // Run pipeline
    results = pipeline(
        samples,
        host_reference,
        regions_bedfile,
        database,
        primer_file,
        min_read_length,
        max_read_length,
        cutsite_csv
        )

    results[0]
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
