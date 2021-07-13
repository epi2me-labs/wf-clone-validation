#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage(){
    log.info """
wf-clone-validation

Usage:
    nextflow run epi2melabs/wf-clone-validation [options]

Script Options:
    --fastq             DIR     Path to directory containing FASTQ files (required)
    --samples           FILE    CSV file with columns named `barcode` and `sample_name`
                                (or simply a sample name for non-multiplexed data).
    --host_reference    FILE    FASTA file, reads which map to it are discarded (optional)
    --regions_bedfile   FILE    BED file, mask regions within host_reference from filtering (optional)
    --approx_size       INT     Approximate size of the plasmid in base pairs (default: 10000)
    --assm_coverage     INT     Try to use this many fold coverage per assembly (default: 150)
    --flye_overlap      INT     Sets the min overlap that flye requires of reads (default: 2000)
    --no_reconcile      BOOL    If enabled, only a single assembly will be made and polished (optional)
    --prefix            STR     The prefix attached to each of the output filenames (optional)
    --help

Notes:
    If directories named "barcode*" are found under the `--fastq` directory the
    data is assumed to be multiplex and each barcode directory will be processed
    independently. If `.fastq(.gz)` files are found under the `--fastq` directory
    the sample is assumed to not be multiplexed. In this second case `--samples`
    should be a simple name rather than a CSV file.
"""
}


def nameIt(ch) {
    return ch.map { it -> return tuple(it.simpleName, it) }.groupTuple()
}


process checkSampleSheet {
    label "artic"
    cpus 1
    input:
        file "sample_sheet.txt"
    output:
        file "samples.txt"
    """
    check_sample_sheet.py sample_sheet.txt samples.txt
    """
}


process combineFastq {
    label "wfplasmid"
    cpus 1
    input:
        tuple file(directory), val(sample_name) 
    output:
        file "${sample_name}.fastq.gz"
    """
    fastcat -x ${directory} | gzip > ${sample_name}.fastq.gz
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
        --keep_percent 90 \
        -t $target \
        $fastq > ${fastq.simpleName}.downsampled.fastq
    """
}


process subsetReads {
    label "wfplasmid"
    input:
        file fastq
    output:
        path "sets/*.fastq", emit: subsets
    shell:
    '''
    trycycler subsample \
        --count 3 \
        --min_read_depth !{(params.assm_coverage / 3) * 2} \
        --reads !{fastq} \
        --out_dir sets \
        --genome_size !{params.approx_size}
    for sub in $(ls sets/sample_*.fastq)
    do
        mv $sub sets/!{fastq.simpleName}.sub$(basename $sub)
    done
    '''
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
        samples
        host_reference
        regions_bedfile
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
        report = buildQCReport(downsampled_fastqs.collect(), 
            polished.collect())
    emit:
        polished = polished
        annotations = annotations
        report = report
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


/**
 * Load a sample sheet into a Nextflow channel to map barcodes
 * to sample names.
 *
 * @param samples CSV file according to MinKNOW sample sheet specification
 * @return A Nextflow Channel of tuples (barcode, sample name)
 */ 
def check_sample_sheet(samples)
{
    println("Checking sample sheet.")
    sample_sheet = Channel.fromPath(samples, checkIfExists: true)
    sample_sheet = checkSampleSheet(sample_sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.barcode, row.sample_name) }
    return sample_sheet
}


/**
 * Find fastq data using various globs. Wrapper around Nextflow `file`
 * method.
 *
 * @param patten glob pattern for top level input folder.
 * @param maxdepth maximum depth to traverse
 * @return list of files.
 */ 
def find_fastq(pattern, maxdepth)
{
    files = []
    extensions = ["fastq", "fastq.gz", "fq", "fq.gz"]
    for (ext in extensions) {
        files += file("${pattern}/*.${ext}", type: 'file', maxdepth: maxdepth)
    }
    return files
}


/**
 * Rework EPI2ME flattened directory structure into standard form
 * files are matched on barcode\d+ and moved into corresponding
 * subdirectories ready for processing.
 *
 * @param input_folder Top-level input directory.
 * @param output_folder Top-level output_directory.
 * @return A File object representating the staging directory created
 *     under output_folder
 */ 
def sanitize_fastq(input_folder, output_folder)
{
    println("Running sanitization.")
    println(" - Moving files: ${input_folder} -> ${output_folder}")
    staging = new File(output_folder)
    staging.mkdirs()
    files = find_fastq("${input_folder}/**/", 1)
    for (fastq in files) {
        fname = fastq.getFileName()
        // find barcode
        pattern = ~/barcode\d+/
        matcher = fname =~ pattern
        if (!matcher.find()) {
            // not barcoded - leave alone
            fastq.renameTo("${staging}/${fname}")
        } else {
            bc_dir = new File("${staging}/${matcher[0]}")
            bc_dir.mkdirs()
            fastq.renameTo("${staging}/${matcher[0]}/${fname}")
        }
    }
    println(" - Finished sanitization.")
    return staging
}


/**
 * Resolves input folder containing barcode subdirectories
 * or a flat set of fastq data to a Nextflow Channel. Removes barcode
 * directories with no fastq files.
 *
 * @param input_folder Top level input folder to locate fastq data
 * @param sample_sheet List of tuples mapping barcode to sample name
 *     or a simple string for non-multiplexed data.
 * @return Channel of tuples (path, sample_name)
 */
def resolve_barcode_structure(input_folder, sample_sheet)
{
    println("Checking input directory structure.")
    barcode_dirs = file("$input_folder/barcode*", type: 'dir', maxdepth: 1)
    not_barcoded = find_fastq("$input_folder/", 1)
    samples = null
    if (barcode_dirs) {
        println(" - Found barcode directories")
        // remove empty barcode_dirs
        valid_barcode_dirs = []
        invalid_barcode_dirs = []
        for (d in barcode_dirs) {
            if(!find_fastq(d, 1)) {
                invalid_barcode_dirs << d
            } else {
                valid_barcode_dirs << d
            }
        }
        if (invalid_barcode_dirs.size() > 0) {
            println(" - Some barcode directories did not contain .fastq(.gz) files:")
            for (d in invalid_barcode_dirs) {
                println("   - ${d}")
            }
        }
        // link sample names to barcode through sample sheet
        if (!sample_sheet) {
            sample_sheet = Channel
                .fromPath(valid_barcode_dirs)
                .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
                .map { path -> tuple(path.baseName, path.baseName) }
        }
        samples = Channel
            .fromPath(valid_barcode_dirs)
            .filter(~/.*barcode[0-9]{1,3}$/)  // up to 192
            .map { path -> tuple(path.baseName, path) }
            .join(sample_sheet)
            .map { barcode, path, sample -> tuple(path, sample) }
    } else if (not_barcoded) {
        println(" - Found fastq files, assuming single sample")
        sample = (sample_sheet == null) ? "unknown" : sample_sheet
        samples = Channel
            .fromPath(input_folder, type: 'dir', maxDepth:1)
            .map { path -> tuple(path, sample) }
    }
    return samples
}


// entrypoint workflow
workflow {

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

    // shadow input, since it may get changed
    input_folder = params.fastq
    // EPI2ME harness 
    if (params.sanitize_fastq) {
        staging = "${params.out_dir}/staging"
        input_folder = sanitize_fastq(input_folder, staging)
    }
    // check sample sheet
    sample_sheet = null
    if (params.samples) {
        sample_sheet = check_sample_sheet(params.samples)
    }
    // resolve whether we have demultiplexed data or single sample
    samples = resolve_barcode_structure(input_folder, sample_sheet)


    host_reference = file(params.host_reference, type: "file")
    regions_bedfile = file(params.regions_bedfile, type: "file")
    // Run pipeline
    results = pipeline(samples, host_reference, regions_bedfile)
    output(results.polished.concat( 
        results.polished, results.report, 
        results.annotations 
    ))
}
