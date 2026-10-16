/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGNMENT } from '../subworkflows/local/alignment'
include { VARIANT_CALLING } from '../subworkflows/local/variant_calling'
include { CNV_CALLING } from '../subworkflows/local/cnv_calling'
include { QC_METRICS_TABLE } from '../modules/local/qc_metrics_table/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_varcall_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARCALL {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Reference genome, read once and reused (as a Nextflow value channel)
    // for every sample by both subworkflows below.
    //
    def ref_meta = [id: file(params.fasta).baseName]
    def ch_fasta = channel.value([ref_meta, file(params.fasta, checkIfExists: true)])

    //
    // SUBWORKFLOW: Align paired-end reads to the reference and generate
    // alignment QC metrics (FastQC + samtools flagstat/stats).
    //
    ALIGNMENT(
        ch_samplesheet,
        ch_fasta,
    )
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.fastqc_zip.map { _meta, qc_file -> qc_file })
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.flagstat.map { _meta, qc_file -> qc_file })
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.stats.map { _meta, qc_file -> qc_file })

    //
    // SUBWORKFLOW: Mark duplicates and call germline short variants with GATK.
    //
    VARIANT_CALLING(
        ALIGNMENT.out.bam_bai,
        ch_fasta,
        ALIGNMENT.out.fasta_fai,
    )
    ch_multiqc_files = ch_multiqc_files.mix(VARIANT_CALLING.out.markdup_metrics.map { _meta, qc_file -> qc_file })

    //
    // SUBWORKFLOW: Call copy-number variants with CNVkit (requires --intervals).
    //
    CNV_CALLING(
        VARIANT_CALLING.out.bam,
        ch_fasta,
        ALIGNMENT.out.fasta_fai,
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'varcall_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'varcall'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    //
    // MODULE: Write the flat, tabular QC metrics file - one of this
    // pipeline's three headline outputs, alongside the BAM and VCF.
    //
    QC_METRICS_TABLE(MULTIQC.out.data)

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    bam            = VARIANT_CALLING.out.bam // channel: [ val(meta), path(bam), path(bai) ] - final, analysis-ready BAM
    vcf            = VARIANT_CALLING.out.vcf // channel: [ val(meta), path(vcf) ]
    cnv_vcf        = CNV_CALLING.out.vcf // channel: [ val(meta), path(vcf) ] - CNV calls
    qc_metrics     = QC_METRICS_TABLE.out.tsv // channel: [ val(meta), path(tsv) ]
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
