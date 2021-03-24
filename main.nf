#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage(){
    log.info """
wf-clone-validation

Usage:
    nextflow run epi2melabs/wf-clone-validation [options]

Script Options:
    --fastq             DIR     Path to directory containing FASTQ files (required)
    --host_reference    FILE    FASTA file, reads which map to it are discarded (optional)
    --regions_bedfile   FILE    BED file, mask regions within host_reference from filtering (optional)
    --approx_size       INT     Approximate size of the plasmid in base pairs (default: 10000)
    --assm_coverage     INT     Try to use this many fold coverage per assembly (default: 150)
    --flye_overlap      INT     Sets the min overlap that flye requires of reads (default: 2000)
    --help

Note:
    --fastq should point to either a directory, or a directory of directories. If the former,
    any fastq files in that single location will be processed together, as a single sample.
    If the latter, i.e. there are sub-directories, each sub-directory will be expected to
    contain .fastq files, and will be treated as a separate sample.
"""
}

if (params.help) {
    helpMessage()
    exit 1
}
if (!params.fastq) {
    helpMessage()
    println("")
    println("Error: `--fastq` is required")
    exit 1
}
if (params.regions_bedfile != "NO_REG_BED" 
    && params.host_reference == "NO_HOST_REF") {
    println("")
    println("Error: `--regions_bedfile` requires `--host_reference` to be set")
    exit 1
}


def nameIt(ch) {
    return ch.map { it -> return tuple(it.simpleName, it) }.groupTuple()
}


process combineFastq {
    cpus params.threads
    input:
        file fastq_dir
    output:
        path "*.fastq", emit: combined
    script:
        catter = fastq_dir.baseName + '.fastq'
    """
    find "$fastq_dir/" -type f -name "*.fastq" -exec cat {} >> $catter \\;
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
    cpus params.threads
    input:
        file fastq
    output:
        path "*.downsampled.fastq", emit: downsampled
    script:
        def target = params.approx_size * params.assm_coverage * 3
    """
    filtlong \
        --min_length 1000 \
        --min_mean_q 12 \
        --keep_percent 90 \
        -t $target \
        $fastq > filtered.temp
    porechop \
        --threads $task.cpus \
        -i filtered.temp \
        -o ${fastq.simpleName}.downsampled.fastq
    rm filtered.temp
    """
}


process subsetReads {
    label "wfplasmid"
    input:
        file fastq
    output:
        path "sets/*.fastq", emit: subsets
    script:
        def name = fastq.simpleName
        def min_depth = (params.assm_coverage / 3) * 2
    """
    trycycler subsample \
        --count 3 \
        --min_read_depth $min_depth \
        --reads $fastq \
        --out_dir sets \
        --genome_size $params.approx_size
    mv sets/sample_01.fastq sets/${name}.subsample01.fastq
    mv sets/sample_02.fastq sets/${name}.subsample02.fastq
    mv sets/sample_03.fastq sets/${name}.subsample03.fastq
    """
}


process assembleFlye {
    label "wfplasmid"
    cpus params.threads
    input:
        file fastq
    output:
        path "*.fasta", emit: assembly
    """
    flye \
        --nano-raw $fastq \
        --meta --plasmids \
        --out-dir assm \
        --min-overlap $params.flye_overlap \
        --threads $task.cpus
    mv assm/assembly.fasta ${fastq.baseName}.fasta
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
        path "*.reconciled.fasta", emit: reconciled
    script:
        def cluster_dir = "trycycler/cluster_001"
    """
    trycycler cluster --assemblies $assemblies --reads $fastq --out_dir trycycler
    trycycler reconcile --reads $fastq --cluster_dir ${cluster_dir}
    trycycler msa --cluster_dir ${cluster_dir}
    trycycler partition --reads $fastq --cluster_dirs ${cluster_dir}
    trycycler consensus --cluster_dir ${cluster_dir}
    mv ${cluster_dir}/7_final_consensus.fasta ${fastq.simpleName}.reconciled.fasta
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
    medaka_consensus -i $fastq -d $draft -m r941_min_high_g360 -o . -t $task.cpus -f
    mv consensus.fasta ${fastq.simpleName}.final.fasta
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


process buildQCReport {
    label "wfplasmid"
    cpus 1
    input:
        file samples
        file assemblies
    output:
        file "report.html"
    shell:
    '''
    # Get maf files for dotplots
    for assm in !{assemblies}
    do
        lastdb ${assm}.lastdb $assm
        lastal ${assm}.lastdb $assm > ${assm}.maf
    done

    # Get summary tables
    seqkit stats -T !{assemblies} > assemblies.tsv
    fastcat --read samples_reads.txt --file samples_summary.txt !{samples} > /dev/null

    aplanat clonevalidation \
        --assembly_summary assemblies.tsv \
        --assembly_mafs *.maf \
        --reads_summary samples_reads.txt \
        --fastq_summary samples_summary.txt
    '''
}


workflow pipeline {
    take:
        fastq_dir
        host_reference
        regions_bedfile
    main:
        // Acquire either top level directory or subdirectories
        // and produce a list of directories each of which will 
        // form a sample
        sample_dirs = channel.fromPath( "$fastq_dir/*", type: 'dir' )
        // Combine fastq from each of the sample directories into 
        // a single per-sample fastq file
        sample_fastqs = combineFastq(sample_dirs.ifEmpty([fastq_dir]))
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
        // And then for each sample split the reads into subsets
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
        // And finally grab the annotations and report
        annotations = prokkaAnnotateAssembly(polished)
        report = buildQCReport(downsampled_fastqs.collect(), 
            polished.collect())
    emit:
        samples = sample_fastqs
        subsets = subsets
        assemblies = assemblies
        deconcatenated = deconcatenated
        reconciled = reconciled
        polished = polished
        annotations = annotations
        report = report
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
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
    // Acquire fastq directory
    fastq_dir = file(params.fastq, type: "dir", checkIfExists: true)
    host_reference = file(params.host_reference, type: "file")
    regions_bedfile = file(params.regions_bedfile, type: "file")
    // Run pipeline
    results = pipeline(fastq_dir, host_reference, regions_bedfile)
    output(results.polished.concat( results.subsets, results.assemblies,
        results.deconcatenated, results.reconciled, results.samples,
        results.report, results.annotations
    ))
}