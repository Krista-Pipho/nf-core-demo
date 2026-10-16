//
// SUBWORKFLOW: VARIANT_CALLING
//
// Takes the coordinate-sorted, indexed BAM produced by the ALIGNMENT
// subworkflow, marks duplicate reads with GATK MarkDuplicates, and calls
// germline short variants (SNVs/indels) per-sample with GATK HaplotypeCaller.
//
// The GATK sequence dictionary (.dict) is built once per pipeline run and
// reused for every sample, unless the user supplies a pre-built one via
// `--dict`.
//
// To call variants with a different tool, swap GATK4_HAPLOTYPECALLER below
// for the equivalent module (e.g. deepvariant) - the mark-duplicates step and
// the rest of the pipeline are unaffected.
//

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP } from '../../../modules/nf-core/samtools/index/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow VARIANT_CALLING {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]  - from ALIGNMENT.out.bam_bai
    ch_fasta // channel: [ val(meta_ref), path(fasta) ]
    ch_fasta_fai // channel: [ val(meta_ref), path(fai) ]          - from ALIGNMENT.out.fasta_fai

    main:

    //
    // Build (or reuse) the GATK/Picard sequence dictionary for the reference genome.
    // This runs once for the whole pipeline run, not once per sample.
    //
    if (params.dict) {
        ch_dict = ch_fasta.map { meta, _fasta -> [meta, file(params.dict, checkIfExists: true)] }
    }
    else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    }

    //
    // MODULE: Mark PCR/optical duplicates on the coordinate-sorted BAM before
    // variant calling. GATK4_MARKDUPLICATES takes the reference as plain
    // (non-tupled) paths, hence the `.map` below.
    //
    GATK4_MARKDUPLICATES(
        ch_bam_bai.map { meta, bam, _bai -> [meta, bam] },
        ch_fasta.map { _meta, fasta -> fasta },
        ch_fasta_fai.map { _meta, fai -> fai },
    )

    //
    // MODULE: Index the deduplicated BAM - this is the final, analysis-ready
    // alignment that HaplotypeCaller calls variants from.
    //
    SAMTOOLS_INDEX_DEDUP(GATK4_MARKDUPLICATES.out.bam)

    ch_dedup_bam_bai = GATK4_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX_DEDUP.out.index)

    //
    // MODULE: Call germline short variants (SNVs/indels) per-sample.
    // `--intervals` (the exome capture BED/interval_list) restricts calling to
    // the covered regions - essential for exome data, where the vast majority
    // of the genome has no coverage at all. The trailing [[:], []] arguments
    // are HaplotypeCaller's optional dbSNP inputs, left empty here.
    //
    def ch_intervals = params.intervals ? file(params.intervals, checkIfExists: true) : []
    GATK4_HAPLOTYPECALLER(
        ch_dedup_bam_bai.map { meta, bam, bai -> [meta, bam, bai, ch_intervals, []] },
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        [[:], []],
        [[:], []],
    )

    emit:
    bam             = ch_dedup_bam_bai // channel: [ val(meta), path(bam), path(bai) ] - final, analysis-ready BAM
    markdup_metrics = GATK4_MARKDUPLICATES.out.metrics // channel: [ val(meta), path(metrics) ]
    vcf             = GATK4_HAPLOTYPECALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi             = GATK4_HAPLOTYPECALLER.out.tbi // channel: [ val(meta), path(tbi) ]
}
