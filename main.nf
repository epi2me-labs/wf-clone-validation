#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 

def helpMessage(){
    log.info """
wf-clone-validation

Usage:
    nextflow run epi2melabs/wf-clone-validation [options]

Script Options:
    --fastq             DIR     Path to directory containing FASTQ files (required)
    --db_directory      DIR     Location of annotation database. (required) 
    --samples           FILE    CSV file with columns named `barcode` and `sample_name`
                                (or simply a sample name for non-multiplexed data).
    --out_dir           DIR     Path to output directory (default: output)
    --host_reference    FILE    FASTA file, reads which map to it are discarded.
    --regions_bedfile   FILE    BED file, mask regions within host_reference from filtering.
    --approx_size       INT     Approximate size of the plasmid in base pairs (default: 7000).
    --assm_coverage     INT     Try to use this many fold coverage per assembly (default: 60).
    --flye_overlap      INT     Sets the min overlap that flye requires of reads (default: 2000).
    --no_reconcile      BOOL    If enabled, only a single assembly will be made and polished.
    --prefix            STR     The prefix attached to each of the output filenames.
    --min_barcode       INT     Minimum number in barcode range.
    --max_barcode       INT     Maximmum number in barcode range.
    --threads           INT     Number of threads per process where applicable (default: 4)
    --help

Notes:
    If directories named "barcode*" are found under the `--fastq` directory the
    data is assumed to be multiplex and each barcode directory will be processed
    independently. If `.fastq(.gz)` files are found under the `--fastq` directory
    the sample is assumed to not be multiplexed. In this second case `--samples`
    should be a simple name rather than a CSV file.

    An generic annotation database can be obtained with the following commands:
        wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-clone-validation/wf-clone-validation-db.tar.gz
        tar -xzvf wf-clone-validation-db.tar.gz
    This will create a directory `wf-clone-validation` in the current working directory
    which can be provided as the `--db_directory` parameter.
"""
}


def nameIt(ch) {
    return ch.map { it -> return tuple(it.simpleName, it) }.groupTuple()
}


process combineFastq {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    input:
        tuple file(directory), val(sample_name) 
    output:
        path "${sample_name}.fastq.gz", optional: true, emit: sample
        path "${sample_name}.stats", emit: stats
        env STATUS, emit: status
    script:
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8
        def expected_length_max = "$params.approx_size" 
        def expected_length_min = "$params.approx_size" 
        int max = (expected_length_max.toInteger()) * 1.5
        int min = (expected_length_min.toInteger()) * 0.5
    """
    STATUS=${sample_name}",Insufficient reads"
    fastcat -a "$min" -b "$max" -s ${sample_name} -r ${sample_name}.stats -x ${directory} > ${sample_name}.fastq
    if [[ "\$(wc -l <"${sample_name}.stats")" -ge "$value" ]];  then 
        gzip ${sample_name}.fastq  
        STATUS=${sample_name}",Pass"
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
    STATUS="${name},Filter host read Fail"
    (minimap2 -t $task.cpus -y -ax map-ont $reference $fastq \
        | samtools sort -o ${name}.sorted.aligned.bam -
    samtools index ${name}.sorted.aligned.bam
    samtools view -b -f 4  ${name}.sorted.aligned.bam > unmapped.bam
    samtools view -b -F 4  ${name}.sorted.aligned.bam > mapped.bam
    samtools fastq unmapped.bam > ${name}.filtered.fastq
    fastcat -s ${name} -r ${name}.stats ${name}.filtered.fastq > /dev/null
    if [[ -f "$regs" ]]; then
        bedtools intersect -a mapped.bam -b $regs -wa \
            | samtools view -bh - > retained.bam
        samtools fastq retained.bam >> ${name}.filtered.fastq
    fi ) && STATUS="${name},Pass"
    """
}


process downsampleReads {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus 1
    input:
        path fastq
    output:
        path "*.downsampled.fastq", optional: true, emit: downsampled
        env STATUS, emit: status
    script:
        def target = params.approx_size * params.assm_coverage * 3
    """
    STATUS=${fastq.simpleName}",Downsample reads"
    (filtlong \
        --min_length 1000 \
        --min_mean_q 12 \
        --keep_percent 99 \
        -t $target \
        $fastq > ${fastq.simpleName}.downsampled.fastq ) && STATUS=${fastq.simpleName}",Pass"
    """
}


process subsetReads {
    errorStrategy 'ignore'
    label "wfplasmid"
    input:
        file fastq
    output:
        path "*.fastq", optional: true, emit: subset
        env STATUS, emit: status
    shell:
    '''
    STATUS="!{fastq.simpleName},Subset reads"
    (trycycler subsample \
        --count 3 \
        --min_read_depth !{(params.assm_coverage / 3) * 2} \
        --reads !{fastq} \
        --out_dir sets \
        --genome_size !{params.approx_size} \
    && for sub in $(ls sets/sample_*.fastq)
    do
        mv $sub ./!{fastq.simpleName}.sub$(basename $sub)
    done) && STATUS="!{fastq.simpleName},Pass"
    '''
}


process assembleFlye {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    input:
        file subset_files  
    output:
        path "*.fasta", optional: true, emit: assembly
        env STATUS, emit: status
    script:
        name = subset_files[0]
    """
    STATUS="${name.simpleName},Assemble Flye"
    (mkdir assm
    for FASTQ in $subset_files
    do
    (flye \
        --nano-raw \$FASTQ \
        --meta --plasmids \
        --out-dir assm \
        --min-overlap $params.flye_overlap \
        --threads $task.cpus && mv assm/assembly.fasta \$FASTQ.fasta)
    done) && STATUS="${name.simpleName},Pass"
    """
}


process deconcatenateAssembly {
    errorStrategy 'ignore'
    label "wfplasmid"
    input:
        file assemblies
    output:
        tuple val(sample_name), path("*.deconcatenated.fasta"), optional: true, emit: deconcatenated
        env STATUS, emit: status
    script:
        name = assemblies[0]
        sample_name = name.toString().split("\\.")[0]
    """
    STATUS="${name.simpleName},Deconcatenate Assembly"
    (for ASSEMBLY in $assemblies
    do
    (deconcatenate.py \$ASSEMBLY -o \$ASSEMBLY.deconcatenated.fasta)
    done) && STATUS="${name.simpleName},Pass"
    """
}


process reconcileAssemblies {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample), file(fastq), file(assemblies)
    output:
        path "*.reconciled.fasta", optional: true, emit: reconciled
        env STATUS, emit: status
    script:
        def cluster_dir = "trycycler/cluster_001" 
    """
    STATUS=${fastq.simpleName}",Reconcile Assemblies"
    (trycycler cluster --assemblies $assemblies --reads $fastq --out_dir trycycler \
    && trycycler reconcile --reads $fastq --cluster_dir ${cluster_dir} \
    && trycycler msa --cluster_dir ${cluster_dir} \
    && trycycler partition --reads $fastq --cluster_dirs ${cluster_dir} \
    && trycycler consensus --cluster_dir ${cluster_dir} \
    && mv ${cluster_dir}/7_final_consensus.fasta ${fastq.simpleName}.reconciled.fasta) \
    && STATUS=${fastq.simpleName}",Pass"
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
    STATUS=${fastq.simpleName}",Medaka Polish Assembly"
    medaka_consensus -i $fastq -d $draft -m r941_min_high_g360 -o . -t $task.cpus -f
    mv consensus.fasta ${fastq.simpleName}.final.fasta
    STATUS=${fastq.simpleName}",Pass"
    """
}


process assemblyStats {
    label "wfplasmid"
    cpus 1
    input:
        file assemblies
    output:
        path "assemblies.tsv", emit: assembly_stat
    """
    seqkit stats -T $assemblies > assemblies.tsv
    """
}


process downsampledStats {
    label "wfplasmid"
    cpus 1 
    input:
        file sample
    output:
        path "*.stats"
    """
    fastcat -s ${sample.simpleName} -r ${sample.simpleName}.stats $sample > /dev/null
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


process sampleStatus {
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


process report {
    label "wfplasmid"
    cpus 1
    input:
        path annotation_database
        path "assemblies/*"
        file "assembly_maf/*"
        path "assembly_stat/*"
        path "downsampled_stats/*"
        file final_status
        path "per_barcode_stats/*"
        path "host_filter_stats/*"
    output:
        path "*report.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "feature_table.txt", emit: feature_table
        path "plannotate.json", emit: plannotate_json
    """ 
    report.py \
    --assembly_summary assembly_stat/* \
    --assembly_mafs assembly_maf/* \
    --downsampled_stats downsampled_stats/* \
    --consensus assemblies/* \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --database $annotation_database \
    --status $final_status \
    --per_barcode_stats per_barcode_stats/* \
    --host_filter_stats host_filter_stats/*
    """
}


workflow pipeline {
    take:
        samples
        host_reference
        regions_bedfile
        database
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
        // After host filtering, we reduce our overall read depth
        // to the desired level
        downsampled_fastqs = downsampleReads(samples_filtered)
        all_status = downsampled_fastqs.status.join(sample_fastqs.status)
        // Now we branch depeneding on whether reconciliaton is
        // enabled, if so we will subset the data and create an
        // assembly for each subset, then use trycyler to reconcile
        // them, this can produce a better, circularised assembly
        // more frequently than not.
        if (!params.no_reconcile) {
            // For each sample split the reads into subsets
            subsets = subsetReads(downsampled_fastqs.downsampled)
            // Assemble each subset independently
            assemblies = assembleFlye(subsets.subset)
            // Deconcatenate assemblies
            deconcatenated = deconcatenateAssembly(assemblies.assembly)
            // Group assemblies back together for reconciliation
            named_samples = nameIt(downsampled_fastqs.downsampled)
            named_deconcatenated = named_samples.join(deconcatenated.deconcatenated)
            reconciled = reconcileAssemblies(named_deconcatenated)
            // Re-group reconciled assemblies together for final polish
            named_reconciled = nameIt(reconciled.reconciled).join(named_samples)
            polished = medakaPolishAssembly(named_reconciled)
            final_status = sample_fastqs.status.join(downsampled_fastqs.status,remainder: true)
                           .join(updated_status,remainder: true)
                           .join(assemblies.status,remainder: true)
                           .join(subsets.status,remainder: true)
                           .join(assemblies.status,remainder: true)
                           .join(deconcatenated.status,remainder: true)
                           .join(reconciled.status, remainder: true)
                           .join(polished.status,remainder: true)
                           .collectFile(name: 'final_status.csv', newLine: true)

            
        } else {
            // Given reconciliation is not enabled
            // Assemble the downsampled dataset in one go
            assemblies = assembleFlye(downsampled_fastqs.downsampled)
            // Deconcatenate assemblies
            deconcatenated = deconcatenateAssembly(assemblies.assembly)
            // Final polish
            named_samples = nameIt(downsampled_fastqs.downsampled)
      
            named_deconcatenated = (deconcatenated.deconcatenated).join(named_samples)
         
            polished = medakaPolishAssembly(named_deconcatenated)
            final_status = sample_fastqs.status.join(downsampled_fastqs.status,remainder: true)
                           .join(updated_status,remainder: true)
                           .join(assemblies.status, remainder: true)
                           .join(deconcatenated.status,remainder: true)
                           .join(polished.status,remainder: true)
                           .collectFile(name: 'final_status.csv', newLine: true)
            
        }
       
        downsampled_stats = downsampledStats(downsampled_fastqs.downsampled)
        assembly_stats = assemblyStats(polished.polished.collect())
        
        assembly_mafs = assemblyMafs(polished.polished.collect())

        report = report(database,
                 polished.polished.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                 assembly_mafs.assembly_maf.ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                 assembly_stats.assembly_stat.ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                 downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
                 final_status,
                 sample_fastqs.stats.collect(),
                 filtered_stats)
        results = polished.polished.concat(
                  report.html,
                  report.sample_stat,
                  report.feature_table,
                  report.plannotate_json)

    emit:
        results
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
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
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq || !params.db_directory) {
        helpMessage()
        println("")
        println("Error: `--fastq` and `--db_directory` are required. A suitable")
        println("database can be obtained using the instructions provided by running `--help`")
        exit 1
    }

    if (params.regions_bedfile != "NO_REG_BED" 
        && params.host_reference == "NO_HOST_REF") {
        println("")
        println("Error: `--regions_bedfile` requires `--host_reference` to be set")
        exit 1
    }

    samples = fastq_ingress(
        params.fastq, workDir, params.samples, params.sanitize_fastq,
        params.min_barcode, params.max_barcode)

    database = file(params.db_directory, type: "dir")
    host_reference = file(params.host_reference, type: "file")
    regions_bedfile = file(params.regions_bedfile, type: "file")
    // Run pipeline
    results = pipeline(samples, host_reference, regions_bedfile, database)

    output(results)
}
