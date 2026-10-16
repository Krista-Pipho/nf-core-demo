//
// MODULE: QC_METRICS_TABLE (local)
//
// Publishes MultiQC's "general stats" table - one row per sample, with the
// FastQC / samtools flagstat / samtools stats / GATK MarkDuplicates columns
// MultiQC already collected - as a plain, flat TSV file. This is the
// pipeline's tabular QC metrics deliverable, alongside the final BAM and VCF.
//
// This is a thin wrapper around a file MultiQC already writes, not a new
// metrics calculation - if you add another QC tool upstream and MultiQC picks
// it up, its columns will appear here automatically with no changes needed.
//
process QC_METRICS_TABLE {
    tag "${meta.id}"
    label 'process_single'

    // Fully-qualified so it resolves from Docker Hub regardless of the
    // `singularity.registry` / `docker.registry` default set in nextflow.config.
    container 'docker.io/library/ubuntu:22.04'

    input:
    tuple val(meta), path(multiqc_data)

    output:
    tuple val(meta), path("*.qc_metrics.tsv"), emit: tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp '${multiqc_data}/multiqc_general_stats.txt' '${prefix}.qc_metrics.tsv'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.qc_metrics.tsv
    """
}
