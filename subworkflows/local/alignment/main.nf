//
// SUBWORKFLOW: ALIGNMENT
//
// Aligns paired-end Illumina reads to a reference genome with BWA-MEM and
// produces alignment QC metrics for every sample.
//
// Reference preparation (the BWA index and the samtools .fai index) is built
// once per pipeline run and reused for every sample, unless the user already
// supplies pre-built files via `--bwa_index` / `--fasta_fai` - handy on an HPC
// system where the same reference is reused across many runs.
//
// To add a different aligner later, swap the BWA_INDEX/BWA_MEM calls below
// for the equivalent modules (e.g. bwamem2/index + bwamem2/mem) and keep the
// rest of the subworkflow (indexing + QC) unchanged.
//

include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGN } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'
include { FASTQC } from '../../../modules/nf-core/fastqc/main'

workflow ALIGNMENT {
    take:
    ch_reads // channel: [ val(meta), [ fastq_1, fastq_2 ] ] - paired-end reads, one entry per sample
    ch_fasta // channel: [ val(meta_ref), path(fasta) ]       - reference genome, a single value reused for every sample

    main:

    //
    // Build (or reuse) the BWA index for the reference genome.
    // This runs once for the whole pipeline run, not once per sample.
    //
    if (params.bwa_index) {
        ch_bwa_index = ch_fasta.map { meta, _fasta -> [meta, file(params.bwa_index, checkIfExists: true)] }
    }
    else {
        BWA_INDEX(ch_fasta)
        ch_bwa_index = BWA_INDEX.out.index
    }

    //
    // Build (or reuse) the samtools .fai index for the reference genome.
    // Also reused downstream by the variant_calling subworkflow, so it is
    // emitted below rather than rebuilt a second time.
    //
    if (params.fasta_fai) {
        ch_fasta_fai = ch_fasta.map { meta, _fasta -> [meta, file(params.fasta_fai, checkIfExists: true)] }
    }
    else {
        SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [meta, fasta, []] }, false)
        ch_fasta_fai = SAMTOOLS_FAIDX.out.fai
    }

    //
    // MODULE: Align paired-end reads with BWA-MEM.
    // sort_bam = true pipes the output straight into `samtools sort`, so the
    // BAM emitted here is already coordinate-sorted - no separate sort step.
    //
    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, true)

    //
    // MODULE: Index the sorted BAM.
    //
    SAMTOOLS_INDEX_ALIGN(BWA_MEM.out.bam)

    ch_bam_bai = BWA_MEM.out.bam.join(SAMTOOLS_INDEX_ALIGN.out.index)

    //
    // Alignment QC: FastQC run directly on the aligned BAM (it auto-detects
    // the BAM format), plus samtools flagstat/stats for mapping-rate,
    // insert-size and coverage-style metrics. All three are collected into
    // the pipeline-wide MultiQC report and the flat QC metrics table.
    //
    FASTQC(BWA_MEM.out.bam)

    SAMTOOLS_FLAGSTAT(ch_bam_bai)

    //
    // Pair every sample's BAM with the single reference fasta+fai via an explicit
    // `.combine()` + re-split `.map()`s, rather than passing a `.join()`-derived
    // reference channel straight in as SAMTOOLS_STATS's second input: a `.join()` of
    // two single-item ("value") channels does NOT reliably broadcast across multiple
    // per-sample invocations - in testing it silently ran the process for only the
    // first sample, with no error at all. `.combine()` (a genuine cross-product) does
    // not have this problem.
    //
    def ch_stats_input = ch_bam_bai
        .combine(ch_fasta.join(ch_fasta_fai))
        .map { meta, bam, bai, ref_meta, fasta, fai -> [meta, bam, bai, ref_meta, fasta, fai] }
    SAMTOOLS_STATS(
        ch_stats_input.map { meta, bam, bai, _ref_meta, _fasta, _fai -> [meta, bam, bai] },
        ch_stats_input.map { _meta, _bam, _bai, ref_meta, fasta, fai -> [ref_meta, fasta, fai] },
    )

    emit:
    bam_bai    = ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    fastqc_zip = FASTQC.out.zip // channel: [ val(meta), path(zip) ]
    flagstat   = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), path(flagstat) ]
    stats      = SAMTOOLS_STATS.out.stats // channel: [ val(meta), path(stats) ]
    fasta_fai  = ch_fasta_fai // channel: [ val(meta_ref), path(fai) ] - reused by variant_calling
}
