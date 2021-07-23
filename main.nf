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
    --host_reference    FILE    FASTA file, reads which map to it are discarded.
    --regions_bedfile   FILE    BED file, mask regions within host_reference from filtering..
    --approx_size       INT     Approximate size of the plasmid in base pairs (default: 10000).
    --assm_coverage     INT     Try to use this many fold coverage per assembly (default: 150).
    --flye_overlap      INT     Sets the min overlap that flye requires of reads (default: 2000).
    --no_reconcile      BOOL    If enabled, only a single assembly will be made and polished.
    --prefix            STR     The prefix attached to each of the output filenames.
    --min_barcode       INT     Minimum number in barcode range.
    --max_barcode       INT     Maximmum number in barcode range.
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
    label "wfplasmid"
    cpus 1
    input:
        tuple file(directory), val(sample_name) 
    output:
        path "${sample_name}.fastq.gz", optional: true, emit: sample
    script:
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8 * 4
    """
    fastcat -x ${directory} > interim.fastq
    if [[ "\$(wc -l <"interim.fastq")" -ge "$value" ]];  then 
        mv interim.fastq ${sample_name}.fastq
        gzip ${sample_name}.fastq  
    fi
    """
}


process filterHostReads {
    label "wfplasmid"
    cpus params.threads
    input:
        file fastq
        file reference
        file regions_bedfile
    output:
        path "*.filtered.fastq", emit: unmapped
    script:
        def name = fastq.simpleName
        def regs = regions_bedfile.name != 'NO_REG_BED' ? regs : 'none'
    """
    minimap2 -t $task.cpus -y -ax map-ont $reference $fastq \
        | samtools sort -o ${name}.sorted.aligned.bam -
    samtools index ${name}.sorted.aligned.bam
    samtools view -b -f 4  ${name}.sorted.aligned.bam > unmapped.bam
    samtools view -b -F 4  ${name}.sorted.aligned.bam > mapped.bam
    samtools fastq unmapped.bam > ${name}.filtered.fastq

    if [[ -f "$regs" ]]; then
        bedtools intersect -a mapped.bam -b $regs -wa \
            | samtools view -bh - > retained.bam
        samtools fastq retained.bam >> ${name}.filtered.fastq
    fi
    """
}


process downsampleReads {
    label "wfplasmid"
    cpus 1
    input:
        path fastq
    output:
        path "*.downsampled.fastq", emit: downsampled
    script:
        def target = params.approx_size * params.assm_coverage * 3
    """
    filtlong \
        --min_length 1000 \
        --min_mean_q 12 \
        --keep_percent 99 \
        -t $target \
        $fastq > ${fastq.simpleName}.downsampled.fastq
    """
}


process subsetReads {
    label "wfplasmid"
    
    input:
        file fastq
    output:
        path "sets/*.fastq", optional: true, emit: subsets 
    shell:
    '''
    (trycycler subsample \
        --count 3 \
        --min_read_depth !{(params.assm_coverage / 3) * 2} \
        --reads !{fastq} \
        --out_dir sets \
        --genome_size !{params.approx_size} || (echo "failed" && exit 1) \
    && for sub in $(ls sets/sample_*.fastq)
    do
        mv $sub sets/!{fastq.simpleName}.sub$(basename $sub)
    done)|| echo "failed something in subsample ${fastq.simpleName}"
    '''
}


process assembleFlye {
    label "wfplasmid"
    cpus params.threads
    input:
        file fastq
    output:
        path "*.fasta", optional: true, emit: assembly
    """
    mkdir assm
    flye \
        --nano-raw $fastq \
        --meta --plasmids \
        --out-dir assm \
        --min-overlap $params.flye_overlap \
        --threads $task.cpus  && mv assm/assembly.fasta ${fastq.baseName}.fasta || echo "failed flye sample ${fastq.simpleName}"
    
    """
}


process deconcatenateAssembly {
    label "wfplasmid"
    input:
        file assembly
    output:
        path "*.deconcatenated.fasta", emit: deconcatenated
    """
    deconcatenate.py $assembly -o ${assembly.baseName}.deconcatenated.fasta
    """
}


process reconcileAssemblies {
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample), file(assemblies), file(fastq)
    output:
        path "*.reconciled.fasta", optional: true, emit: reconciled
    script:
        def cluster_dir = "trycycler/cluster_001" 
    """
    
    (trycycler cluster --assemblies $assemblies --reads $fastq --out_dir trycycler  || (echo "failed" && exit 1)) \
    && (trycycler reconcile --reads $fastq --cluster_dir ${cluster_dir} || (echo "failed" && exit 1)) \
    && (trycycler msa --cluster_dir ${cluster_dir} || (echo "failed" && exit 1 )) \
    && (trycycler partition --reads $fastq --cluster_dirs ${cluster_dir} || (echo "failed" && exit 1 )) \
    && (trycycler consensus --cluster_dir ${cluster_dir} && mv ${cluster_dir}/7_final_consensus.fasta ${fastq.simpleName}.reconciled.fasta) \
    || echo "failed something in Trycycler ${fastq.simpleName}"
 
    """
}


process medakaPolishAssembly {
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample), file(draft), file(fastq)
    output:
        path "*.final.fasta", emit: polished

    """
    medaka_consensus -i $fastq -d $draft -m r941_min_high_g360 -o . -t $task.cpus -f && mv consensus.fasta ${fastq.simpleName}.final.fasta
    """
}


process prokkaAnnotateAssembly {
    label "wfplasmid"
    cpus params.threads
    input:
        file assembly
    output:
        path "*_prokka", type: "dir", emit: annotations
    """
    prokka --outdir ${assembly.simpleName}_prokka --prefix ${assembly.simpleName} $assembly
    """
}


process assemblyStats {
    label "wfplasmid"
    cpus 1
    input:
        file samples
        file assemblies

    output:
        path "assemblies.tsv", emit: assembly_stat
        path "samples_reads.txt", emit: samples_reads
        path "samples_summary.txt", emit: sample_summary
    """
    seqkit stats -T $assemblies > assemblies.tsv
    fastcat --read samples_reads.txt --file samples_summary.txt $samples
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


process report {
    label "wfplasmid"
    cpus 1
    input:
        path annotation_database
        file assemblies
        path assembly_maf
        path assembly_stat
        path samples_reads
        path sample_summary
        file sample_status
        
    output:
        path "*report.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "feature_table.txt", emit: feature_table
    script:
        def sample_sheet = projectDir + "/${params.samples}"
    """ 
    report.py \
    --assembly_summary $assembly_stat \
    --assembly_mafs $assembly_maf \
    --reads_summary $samples_reads \
    --fastq_summary $sample_summary \
    --consensus $assemblies \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --database $annotation_database \
    --status_sheet $sample_status \
    --sample_sheet $sample_sheet
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
                sample_fastqs, host_reference, regions_bedfile)
            sample_fastqs = filtered.unmapped
        }
        // After host filtering, we reduce our overall read depth
        // to the desired level
        downsampled_fastqs = downsampleReads(sample_fastqs)

        // Now we branch depeneding on whether reconciliaton is
        // enabled, if so we will subset the data and create an
        // assembly for each subset, then use trycyler to reconcile
        // them, this can produce a better, circularised assembly
        // more frequently than not.
        if (!params.no_reconcile) {
            // For each sample split the reads into subsets
            subsets = subsetReads(downsampled_fastqs)
            // Assemble each subset independently
            assemblies = assembleFlye(subsets.flatten())
            // Deconcatenate assemblies
            deconcatenated = deconcatenateAssembly(assemblies)
            // Group assemblies back together for reconciliation
            named_samples = nameIt(downsampled_fastqs)
            named_deconcatenated = nameIt(deconcatenated).join(named_samples)

            reconciled = reconcileAssemblies(named_deconcatenated)
            // Re-group reconciled assemblies together for final polish
            named_reconciled = nameIt(reconciled).join(named_samples)
            polished = medakaPolishAssembly(named_reconciled)
        } else {
            // Given reconciliation is not enabled
            // Assemble the downsampled dataset in one go
            assemblies = assembleFlye(downsampled_fastqs)
            // Deconcatenate assemblies
            deconcatenated = deconcatenateAssembly(assemblies)
            // Final polish
            named_samples = nameIt(downsampled_fastqs)
            named_deconcatenated = nameIt(deconcatenated).join(named_samples)
            polished = medakaPolishAssembly(named_deconcatenated)
        }
        // And finally grab the annotations and report
        annotations = prokkaAnnotateAssembly(polished)

        assembly_stats = assemblyStats(downsampled_fastqs.collect(), 
            polished.collect())
        assembly_mafs = assemblyMafs(polished.collect())

        //plannotate
        sample_status = workDir + '/sample_status.txt'

        report = report(database,polished.collect(),assembly_mafs.assembly_maf,
                        assembly_stats.assembly_stat,assembly_stats.samples_reads,
                        assembly_stats.sample_summary, sample_status)
        
        results = polished.concat(
            polished,
            annotations,
            report.html,
            report.sample_stat,
            report.feature_table)

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
        println("database can be obtained using `--download_db`.")
        exit 1
    }

    if (params.regions_bedfile != "NO_REG_BED" 
        && params.host_reference == "NO_HOST_REF") {
        println("")
        println("Error: `--regions_bedfile` requires `--host_reference` to be set")
        exit 1
    }
    // create status file
    barcode_dirs = file("$params.fastq/barcode*", type: 'dir', maxdepth: 1)
    sample_status = workDir + '/sample_status.txt'
    sample_status.write "Sample\tPass/fail\n"
    if(barcode_dirs) {
        for (d in barcode_dirs) {
        sample_status << "${d}" + "\n" 
    }}


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
