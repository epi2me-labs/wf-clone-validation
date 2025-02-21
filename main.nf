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


process filterHostReads {
    errorStrategy 'ignore'
    label "wfplasmid"
    cpus params.threads
    memory "4GB"
    input:
        tuple val(meta), path(fastq), path(fastcat_stats), path(reference), path(regions_bedfile)
    output:
        tuple val(meta), path("*.filtered.fastq"), path(fastcat_stats), optional: true, emit: unmapped
        tuple val(meta), path("${meta.alias}.host.bam"), path("${meta.alias}.host.bam.bai"), optional: true, emit: host_bam
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
    samtools view --bam --require-flags 4 --write-index -o "unmapped.bam##idx##unmapped.bam.bai" --output-unselected  "${name}.host.bam##idx##${name}.host.bam.bai" "${name}.sorted.aligned.bam"
    samtools fastq unmapped.bam > ${name}.filtered.fastq

    # Run fastcat on filtered reads
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


process checkIfEnoughReads {
    label "wfplasmid"
    cpus params.threads
    memory "2GB"
    input:
        tuple val(meta), path("input.fastq.gz"), path("fastcat_stats")
    output:
        tuple val(meta), path("${meta.alias}.fastq.gz"), optional: true, emit: sample
        path("${meta.alias}.fastcat_stats"), emit: stats
        tuple val(meta.alias), env(STATUS), emit: status
    script:
        int min_read_length = params.large_construct ? 200 : meta.approx_size * 0.5
        int max_read_length = meta.approx_size * 1.5
        def extra_args = "-a $min_read_length -b $max_read_length"
        def expected_depth = "$params.assm_coverage"
        // a little heuristic to decide if we have enough data
        int value = (expected_depth.toInteger()) * 0.8
        int bgzip_threads = task.cpus == 1 ? 1 : task.cpus - 1
    """
    STATUS="Failed due to insufficient reads"
    # Raw fastcat stats - carried forward from filterHost (if run)
    mv fastcat_stats "${meta.alias}.fastcat_stats"

    fastcat -s ${meta.alias} -r ${meta.alias}.interim $extra_args input.fastq.gz \
    | bgzip -@ $bgzip_threads > interim.fastq.gz
    if [[ "\$(wc -l < "${meta.alias}.interim")" -ge "$value" ]]; then
        mv interim.fastq.gz ${meta.alias}.fastq.gz
        STATUS="Completed successfully"
    fi
    """
}


process medakaPolishAssembly {
    label "medaka"
    cpus params.threads
    memory "4GB"
    input:
        tuple val(meta), path(draft), path(fastq), path(medaka_model)
    output:
        tuple val(meta.alias), env(STATUS), emit: status
        tuple val(meta), path("${meta.alias}.final.fastq"), emit: fastq
    script:
    // we use `params.override_basecaller_cfg` if present; otherwise use
    // `meta.basecall_models[0]` (there should only be one value in the list because
    // we're running ingress with `allow_multiple_basecall_models: false`; note that
    // `[0]` on an empty list returns `null`)
    String basecall_model = params.override_basecaller_cfg ?: meta.basecall_models[0]
    if (!basecall_model && !params.medaka_model_path) {
        error "Found no basecall model information in the input data for " + \
            "sample '$meta.alias'. Please provide it with the " + \
            "`--override_basecaller_cfg` parameter"
    }
    String model_input =  "${basecall_model}:consensus"
    if (params.medaka_model_path){
        model_input = "${medaka_model}"
    }
    """
    STATUS="Failed to polish assembly with Medaka"
    medaka_consensus -i "${fastq}" -d "${draft}" -m "${model_input}" \
        -o . -t $task.cpus -f -q --bacteria
    echo ">${meta.alias}" >> "${meta.alias}.final.fasta"
    sed "2q;d" consensus.fastq >> "${meta.alias}.final.fasta"
    mv consensus.fastq "${meta.alias}.final.fastq"
    STATUS="Completed successfully"
    """
}


process reorientateFastqAndGetFasta {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path("assembly.fastq"), path(full_reference)
    output:
        tuple val(meta), path(reoriented_fasta), emit: fasta
        tuple val(meta), path(reoriented_fastq), emit: fastq
    script:
        reoriented_fastq = "${meta.alias}.final.fastq"
        reoriented_fasta = "${meta.alias}.final.fasta"
    """
    (
        set -euo pipefail
        if [ -e "OPTIONAL_FILE" ]; then
            # there is no full ref for this sample
            echo "No full reference; skipping reorientation"
            exit 0
        fi

        # duplicate the assembly to avoid potential issues with soft-clipping, repeats etc
        seqkit concat assembly.fastq assembly.fastq -w0 > double-assembly.fastq

        # align ref against the duplicated assembly and get alignment stats
        minimap2 -ax asm5 double-assembly.fastq $full_reference | bamstats - > bamstats.tsv

        # bamstats ignores non-primary alignments; so there is either one alignment (i.e.
        # two lines) or none (one line) ine the file
        if [[ \$(wc -l < bamstats.tsv) == 1 ]]; then
            echo "Reference didn't map against the assembly; skipping reorientation."
            exit 0
        fi

        # get the length of the assembly
        assembly_length=\$(awk 'NR == 2 {print length}' assembly.fastq)

        ref_length=\$(awk 'NR == 2 {print length}' $full_reference)

        # get the strand ('+' or '-') and offset to rotate the assembly
        read -r strand offset <<< "\$(awk -F '\\t' -v assembly_length="\$assembly_length" -v ref_length="\$ref_length" '
            NR == 1 { for (i=1; i<=NF; i++) { indices[\$i] = i } }
            NR == 2 {
                strand = \$indices["direction"]
                rstart = \$indices["rstart"]
            }
            END {
                rstart += 1
                if (strand == "-") {
                    # For negative strand, move to the end of the reference
                    rstart = rstart + ref_length
                    if (rstart > assembly_length) rstart -= assembly_length
                }
                if (rstart > assembly_length) rstart -= assembly_length
                print strand, rstart
            }
        ' bamstats.tsv)"

        # reorient the assembly
        if [[ \$strand == "+" ]]; then
            seqkit restart -w0 -i \$offset assembly.fastq > $reoriented_fastq
        elif [[ \$strand == "-" ]]; then
            # the reference aligned against the negative strand of the assembly; so
            # first rotate the sequence and then reverse complement it
            seqkit restart -w0 -i \$offset assembly.fastq |
                seqkit seq -rpvt dna - > $reoriented_fastq
        else
            echo "Found unexpected strand '\$strand'; skipping reorientation"
            exit 0
        fi
    )

    if [[ ! -f $reoriented_fastq ]]; then
        echo "Reorientation failed"
        mv assembly.fastq $reoriented_fastq
    fi
    # rename sequence to sample alias
    sed -i '1s/^.*\$/@$meta.alias/' $reoriented_fastq
    seqkit fq2fa $reoriented_fastq > $reoriented_fasta
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
    cat !{sequence} | seqkit amplicon -p !{primers} -m !{params.primer_mismatch} -j !{task.cpus} --bed >> !{meta.alias}.interim
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
    cache false
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
    memory params.large_construct ? "7GB" : "2GB"
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
    memory "4GB"
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
    workflow-glue find_inserts --output $output --assemblies "assemblies" \$inserts $ref $large_construct
    """
}


process insert_qc {
    label "wfplasmid"
    cpus 1
    memory "2GB"
    input:
         tuple val(meta), path("insert_assembly.fasta"), path("insert_ref.fasta")
    output:
        tuple val(meta), path("${meta.alias}.insert.calls.bcf"), path("${meta.alias}.insert.calls.bcf.csi"), path("${meta.alias}.insert.stats"), optional: true, emit: insert_stats
        tuple val(meta.alias), env(STATUS), emit: status
         path("${meta.alias}.insert.bam.stats"), emit: insert_align_stats
    script:
    """
    STATUS="Insert found but does not align with provided reference"
    minimap2 -t $task.cpus -y -ax map-ont insert_ref.fasta "insert_assembly.fasta" | samtools sort -o output.bam -
    bamstats -t ${task.cpus} -s "${meta.alias}" "output.bam" > "${meta.alias}.insert.bam.stats"
    mapped=\$(samtools view -F 4 -c  output.bam)
    if [ \$mapped != 0 ]; then
        # -m1 reduces evidence for indels to 1 read instead of the default 2
        bcftools mpileup -m1 -Ou -f insert_ref.fasta output.bam | bcftools call -mv -Ob -o "${meta.alias}.insert.calls.bcf"
        bcftools index "${meta.alias}.insert.calls.bcf"
        bcftools stats "${meta.alias}.insert.calls.bcf" > "${meta.alias}.insert.stats"
        STATUS="Completed successfully"
    fi
    """
}


process align_assembly {
    label "wfplasmid"
    cpus 3
    memory "7GB"
    input:
        tuple val(meta), path("full_assembly.fasta"), path("full_reference.fasta")
    output:
        tuple val(meta), path("${meta.alias}.bam.stats"), emit: bamstats
        tuple val(meta), path("${meta.alias}.bam"), path("${meta.alias}.bam.bai"), path("full_reference.fasta"), emit: ref_aligned, optional: true
    script:
    // Align assembly with the reference. The assembly has previously been attempted to be reorientated to match with the reference.
    """
    minimap2 -t ${task.cpus - 1} -y -ax asm5 full_reference.fasta "full_assembly.fasta" > ${meta.alias}_unsorted.bam
    mapped=\$(samtools view -F 4 -c ${meta.alias}_unsorted.bam)
    bamstats -t ${task.cpus} -s "${meta.alias}" "${meta.alias}_unsorted.bam" > "${meta.alias}.bam.stats"
    if [ \$mapped != 0 ]; then
        samtools sort -@ ${task.cpus} --write-index -o ${meta.alias}.bam##idx##${meta.alias}.bam.bai ${meta.alias}_unsorted.bam
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
        tuple val(meta), path("${meta.alias}.full_construct.calls.bcf"), path("${meta.alias}.full_construct.calls.bcf.csi"), path("${meta.alias}.full_construct.stats"), emit: full_assembly_stats
    script:
    // Also get variants report
    """
    # -m1 reduces evidence for indels to 1 read instead of the default 2
    # -P 1 increases the sensitivity by setting a very high mutation rate prior
    bcftools mpileup -m1 -Ou -f full_reference.fasta "${meta.alias}.bam" | bcftools call -mv -Ob -P 1 -o "${meta.alias}.full_construct.calls.bcf"
    bcftools index "${meta.alias}.full_construct.calls.bcf"
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
        tuple val(meta), path("${meta.alias}.assembly_stats.tsv")
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
        path full_assembly_variants, stageAs: "full_assembly_variants/*"
        path assembly_bamstats, stageAs: "assembly_bamstats/*"
        path cutsite_csv, stageAs: "cutsite_csv/*"
        path insert_bamstats, stageAs: "insert_bamstats/*"
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
        String insert_bamstats = insert_bamstats.fileName.name == "OPTIONAL_FILE" ? "" : "--insert_alignment_bamstats insert_bamstats/*"
        String assembly_bamstats = assembly_bamstats.fileName.name == "OPTIONAL_FILE" ? "" : "--reference_alignment_bamstats assembly_bamstats/*"
        String full_assembly_variants = full_assembly_variants.fileName.name == "OPTIONAL_FILE" ? "" : "--full_assembly_variants full_assembly_variants"
    """
    echo '${metadata}' > metadata.json
  

   
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
    $assembly_bamstats \
    $full_assembly_variants \
    $cutsite \
    $insert_bamstats \
    --expected_coverage ${params.expected_coverage} \
    --expected_identity ${params.expected_identity} 
    """
}

workflow pipeline {
    take:
        samples
        database
        primers
        cutsite_csv // alias, read_counts, cutsite_count
        user_medaka_model

    main:

        // Optionally filter the data, removing reads mapping to
        // the host or background genome

        samples
        | branch {meta, fastq, stats ->
            host: meta.host_reference != "NO_HOST_REF"
                return [meta, fastq, stats, meta.host_reference, meta.regions_bedfile]
            no_host: meta.host_reference == "NO_HOST_REF"
                return [meta, fastq, stats]
        }
        | set{ host_samples }
        
        filtered = filterHostReads(host_samples.host)
        samples_filtered = filtered.unmapped.concat(host_samples.no_host)
        host_bams = filtered.host_bam
        host_status = filtered.status
        filtered_stats = filtered.host_filter_stats | collect | ifEmpty(OPTIONAL_FILE)

        // drop samples with too low coverage
        sample_fastqs = checkIfEnoughReads(samples_filtered)


        // Core assembly and reconciliation
        assemblies = assembleCore(sample_fastqs.sample)
        named_drafts = assemblies.assembly.groupTuple()
        named_samples = assemblies.downsampled.groupTuple()
        named_drafts_samples = named_drafts
        .join(named_samples, failOnMismatch: false, remainder: true)
        .filter{alias, assembly, fastq -> (assembly != null)}
        .combine(user_medaka_model)

        // Polish draft assembly
        polished = medakaPolishAssembly(named_drafts_samples)
        reorientated = reorientateFastqAndGetFasta(
            polished.fastq | map { meta, assembly ->
                [meta, assembly, meta.full_reference ?: OPTIONAL_FILE]
            }
        )

        downsampled_stats = downsampledStats(assemblies.downsampled)

        primer_beds = findPrimers(primers, reorientated.fasta)
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


        // Assembly QC
        assembly_quality = assembly_qc(reorientated.fastq)
        
        // Find inserts
        insert = inserts(ref_groups)
        // Inserts are processed together (to create MSA) so flatten and add key 

        insert_tuple = insert.inserts.flatten().map{ assembly -> tuple(assembly.simpleName, assembly)}
        // Keep assemblies that have an insert_reference and an insert assembly for insert_qc
        insert_meta = reorientated.fasta.map{ 
            meta, polished -> [meta.alias, meta]}
            .join(insert_tuple, failOnMismatch: false, remainder: true)
            .filter{ key, meta, insert -> (insert != null) && meta.insert_reference}
            .map{ key, meta, insert -> [meta, insert, meta.insert_reference]} 
        // Insert QC
        insert_qc_tuple = insert_qc(insert_meta)

        qc_insert = insert_qc_tuple.insert_stats.map{meta, bcf, idx, stats -> stats}
        bcf_insert = insert_qc_tuple.insert_stats.map{meta, bcf, idx, stats -> bcf}
        insert_status = insert_qc_tuple.status
        
        // Full reference QC
        qc_full = align_assembly( reorientated.fasta
        | filter {meta, assembly -> meta.full_reference}
        | map {meta, assembly -> [meta, assembly, meta.full_reference ]})
        assembly_comparison_output = assembly_comparison( qc_full.ref_aligned )
        ref_bamstats = qc_full.bamstats
        full_assembly_stats = assembly_comparison_output.full_assembly_stats.map{meta, bcf, idx, stats -> stats}
        bcf = assembly_comparison_output.full_assembly_stats.map{meta, bcf, idx, stats -> bcf}
        bam = qc_full.ref_aligned.map{meta, bam, bai, ref -> [meta, bam, bai]}


        // Concat statuses and keep the last of each
        final_status = host_status.concat(sample_fastqs.status)
        .concat(assemblies.status).concat(polished.status).concat(insert_status).groupTuple()
        .map { it -> it[0].toString() + ',' + it[1][-1].toString() }
        final_status = final_status.collectFile(name: 'final_status.csv', newLine: true)

        annotation = runPlannotate(
            database, reorientated.fasta.map { it -> it[1] }.collect().ifEmpty(OPTIONAL_FILE),
            final_status)

        mafs = assemblyMafs(reorientated.fasta)

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
            assembly_quality.map{ meta, stats -> stats }.collect().ifEmpty(OPTIONAL_FILE),
            mafs.map{ meta, maf -> maf}.collect().ifEmpty(OPTIONAL_FILE),
            client_fields,
            workflow.manifest.version,
            full_assembly_stats.collect().ifEmpty(OPTIONAL_FILE),
            ref_bamstats.map{ meta, bamstats -> bamstats }.collect().ifEmpty(OPTIONAL_FILE),
            cutsite_csv,
            insert_qc_tuple.insert_align_stats.collect().ifEmpty(OPTIONAL_FILE)
            )

    emit:
        report = report | concat // aggregated channel
        annotations = annotation | concat // aggregated channel
        assemblies = reorientated | join // per-sample channel
        inserts = insert.inserts // per-insert channel
        insert_qc = insert_qc_tuple.insert_stats // per-sample channel
        bams = bam //per-sample channel
        full_assembly_stats = assembly_comparison_output.full_assembly_stats // per-sample channel
        bamstats = ref_bamstats // per-sample channel
        assembly_stats = assembly_quality // per-sample channel
        mafs = mafs // per-sample channel
        host = host_bams // per-sample channel
        telemetry = workflow_params // run specific
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish {
    // publish inputs to output directory
    label "wfplasmid"
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
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.containsKey("reference")) {
        throw new Exception("--reference is deprecated, use --insert_reference instead.")
    }

    // warn the user if overriding the basecall models found in the inputs
    if (params.override_basecaller_cfg) {
        log.warn \
            "Overriding basecall model with '${params.override_basecaller_cfg}'."
    }

    // warn the user if overriding the medaka model
    // which will also override the override_basecaller_cfg param
    if (params.medaka_model_path) {
            log.warn "Overriding the automatically selected medaka model with '${params.medaka_model_path}'."
            user_medaka_model = Channel.fromPath(params.medaka_model_path, type: "file", checkIfExists: true)
    } else {
        user_medaka_model = Channel.fromPath(OPTIONAL_FILE)
    }

    if (params.assembly_tool == 'canu'){
        log.warn "Assembly tool Canu. This may result in suboptimal performance on ARM devices."
    }

    Map ingress_args = [
        "sample": params.sample,
        "sample_sheet": params.sample_sheet,
        "analyse_unclassified": params.analyse_unclassified,
        "stats": true,
        "allow_multiple_basecall_models": false,
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

        // Host reference checks
        if (meta.host_reference && params.host_reference){
            throw new Exception("Parameter --host_reference cannot be used in conjunction with a sample sheet with column 'host_reference'")
        }
        // Host regions bed
        if (meta.regions_bedfile && params.regions_bedfile){
            throw new Exception("Parameter --regions_bedfile cannot be used in conjunction with a sample sheet with column 'regions_bedfile'")
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
	    full_reference = (meta.full_reference ?: params.full_reference) 
        insert_reference = (meta.insert_reference ?: params.insert_reference)
        host_reference = (meta.host_reference ?: params.host_reference)
        regions_bedfile = (meta.regions_bedfile ?: params.regions_bedfile)
        approx_size = (meta.approx_size ?: params.approx_size) as int
        // Check if null before adding file around it ()
        new_keys = ["full_reference": full_reference ? file(full_reference, checkIfExists: true) : null,
            "insert_reference": insert_reference ? file(insert_reference, checkIfExists: true) : null,
            "host_reference": host_reference ? file(host_reference, checkIfExists: true) : "NO_HOST_REF",
            "regions_bedfile": regions_bedfile ? file(regions_bedfile, checkIfExists: true) : file("NO_REG_BED", checkIfExists: false),
            "approx_size": approx_size
        ]
        return [ meta + new_keys, read, stats ]
    }


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
        database,
        primer_file,
        cutsite_csv,
        user_medaka_model
        )
        
    results.report
    | concat(
        results.annotations | flatten ,
        results.inserts | collect ,
        results.assemblies | map { meta, fa, fq -> [fa, fq] }| flatten ,
        results.insert_qc | map { meta, bcf, idx, stats -> [bcf, idx, stats] } | flatten , 
        results.full_assembly_stats | map { meta, bcf, idx, stats -> [bcf, idx, stats ] } | flatten ,
        results.bamstats | map { meta, bamstats -> bamstats }, 
        results.assembly_stats | map { meta, stats -> stats },
        results.mafs | map { meta, mafs -> mafs }, 
        results.host | map { meta, bam, bai -> [ bam, bai ] } | flatten ,
        results.bams | map { meta, bam, bai -> [ bam, bai ] } | flatten )
    | publish
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
