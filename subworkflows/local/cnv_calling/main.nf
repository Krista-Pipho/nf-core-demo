//
// SUBWORKFLOW: CNV_CALLING
//
// Calls copy-number variants from the deduplicated, analysis-ready BAM
// produced by the variant_calling subworkflow, using CNVkit in
// germline/tumor-only mode: no matched normal sample and no pre-built
// reference.cnn, so CNVkit builds a "flat" reference internally (assumes
// uniform coverage across the capture targets). This is the standard
// approach for single-sample or small-cohort exome CNV calling when you
// don't have a panel of normals - see CNVkit's docs on "tumor-only" mode.
//
// Requires --intervals (the exome capture regions, as a BED file - CNVkit
// needs BED, not a GATK .interval_list) to build target/antitarget bins
// that make sense for capture data; without it, CNVkit would fall back to
// whole-genome-style bins that don't reflect where exome data actually has
// coverage.
//
// To call CNVs with a panel of normals instead of a flat reference (more
// accurate, but requires a cohort of normal samples), build a reference.cnn
// with `cnvkit/reference` and pass it via CNVKIT_BATCH's `reference` input
// instead of `[[:], []]` below.
//

include { CNVKIT_BATCH } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_EXPORT } from '../../../modules/nf-core/cnvkit/export/main'

workflow CNV_CALLING {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]  - from VARIANT_CALLING.out.bam (deduplicated BAM)
    ch_fasta // channel: [ val(meta_ref), path(fasta) ]
    ch_fasta_fai // channel: [ val(meta_ref), path(fai) ]          - from ALIGNMENT.out.fasta_fai

    main:

    if (!params.intervals) {
        error("CNV_CALLING requires --intervals (the exome capture regions, as a BED file) to build sensible target/antitarget bins.")
    }
    // A plain file path (not a channel): params.intervals is a CLI param available
    // synchronously, so it can just be captured by the closures below directly.
    def targets_file = file(params.intervals, checkIfExists: true)

    //
    // MODULE: CNVkit batch - target/antitarget binning, flat reference, coverage,
    // ratio ("fix") and segmentation ("segment" + "call") all in one wrapped command.
    //
    // Pair every sample's BAM with the single reference fasta+fai via an explicit
    // `.combine()` + re-split `.map()`s, rather than passing a `.join()`-derived
    // reference channel straight in as a separate process input: a `.join()` of two
    // single-item ("value") channels does NOT reliably broadcast across multiple
    // per-sample invocations - in testing it silently ran the process for only the
    // first sample, with no error at all. `.combine()` (a genuine cross-product) does
    // not have this problem. See the identical fix/comment on ALIGNMENT's
    // SAMTOOLS_STATS call.
    //
    def ch_with_ref = ch_bam_bai
        .combine(ch_fasta.join(ch_fasta_fai))
        .map { meta, bam, bai, _ref_meta, fasta, fai -> [meta, bam, bai, fasta, fai] }
    CNVKIT_BATCH(
        ch_with_ref.map { meta, bam, bai, _fasta, _fai -> [meta, bam, bai, [], []] },
        ch_with_ref.map { meta, _bam, _bai, fasta, fai -> [[id: meta.id], fasta, fai] },
        ch_with_ref.map { meta, _bam, _bai, _fasta, _fai -> [[id: meta.id], targets_file] },
        [[:], []],
        false,
    )

    //
    // CNVKIT_BATCH's `*.cns` glob matches both the raw segments ("<id>.cns") and the
    // integer-copy-number calls ("<id>.call.cns") in the same emission - pick out the
    // call file specifically for VCF export.
    //
    def ch_call_cns = CNVKIT_BATCH.out.cns.map { meta, cns_files ->
        def call_cns = (cns_files instanceof List ? cns_files : [cns_files]).find { file -> file.name.endsWith('.call.cns') }
        [meta, call_cns]
    }

    //
    // MODULE: Export the integer copy-number calls to VCF, for consistency with the
    // variant_calling subworkflow's VCF output.
    //
    CNVKIT_EXPORT(ch_call_cns)

    emit:
    cnr = CNVKIT_BATCH.out.cnr // channel: [ val(meta), path(cnr) ]      - bin-level copy-number ratios
    cns = CNVKIT_BATCH.out.cns // channel: [ val(meta), path(cns) ]      - segments (raw + integer-call)
    vcf = CNVKIT_EXPORT.out.output // channel: [ val(meta), path(vcf) ]      - CNV calls as VCF
}
