# kpipho/varcall: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory (`--outdir`) after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Alignment](#alignment) (`subworkflows/local/alignment`) - BWA-MEM alignment to the reference genome, plus alignment QC (FastQC, samtools flagstat/stats)
- [Variant calling](#variant-calling) (`subworkflows/local/variant_calling`) - GATK MarkDuplicates and GATK HaplotypeCaller germline variant calling
- [CNV calling](#cnv-calling) (`subworkflows/local/cnv_calling`) - CNVkit copy-number calling
- [QC metrics table](#qc-metrics-table) - a single flat, tabular QC summary file
- [MultiQC](#multiqc) - Aggregate HTML report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `reference/`
  - `bwa/`: BWA index for the reference genome (only written if `--bwa_index` was not supplied).
  - `*.fai`: samtools FASTA index (only written if `--fasta_fai` was not supplied).
  - `*.dict`: GATK/Picard sequence dictionary (only written if `--dict` was not supplied).
- `alignment/aligned_bam/`
  - `*.bam`, `*.bam.bai`: coordinate-sorted BAM from BWA-MEM, before duplicate marking.
- `alignment/fastqc/`
  - `*_fastqc.html`, `*_fastqc.zip`: FastQC report run directly on the aligned BAM.
- `alignment/samtools_flagstat/`
  - `*.flagstat`: read mapping/pairing counts from `samtools flagstat`.
- `alignment/samtools_stats/`
  - `*.stats`: detailed alignment statistics (insert size, error rate, coverage, etc.) from `samtools stats`.

</details>

Paired-end reads are aligned to the reference genome with [BWA-MEM](https://github.com/lh3/bwa), piped directly into `samtools sort` so the output BAM is already coordinate-sorted. [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and `samtools flagstat`/`stats` then summarise alignment quality; all three feed into the MultiQC report and the flat QC metrics table below.

### Variant calling

<details markdown="1">
<summary>Output files</summary>

- **`bam/`** (headline output)
  - `*.bam`, `*.bam.bai`: final, deduplicated, analysis-ready BAM used for variant calling.
- `variant_calling/markduplicates_metrics/`
  - `*.metrics`: duplication metrics from GATK MarkDuplicates.
- **`vcf/`** (headline output)
  - `*.vcf.gz`, `*.vcf.gz.tbi`: germline short variants (SNVs/indels) called by GATK HaplotypeCaller.

</details>

[GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard) flags PCR/optical duplicate reads on the aligned BAM. The deduplicated BAM is then used by [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) to call germline SNVs and indels per sample. If `--intervals` is set, calling is restricted to those regions (+100bp padding) - essential for exome data, where most of the genome has no coverage at all.

### CNV calling

<details markdown="1">
<summary>Output files</summary>

- **`cnv/`** (headline output)
  - `*.targetcoverage.cnn`, `*.antitargetcoverage.cnn`, `*.cnr`: per-bin coverage and copy-number ratios.
  - `*.cns`: raw copy-number segments.
  - `*.call.cns`: segments with integer copy-number calls.
  - `*.vcf`: the integer copy-number calls, converted to VCF.

</details>

[CNVkit](https://cnvkit.readthedocs.io/) calls copy-number variants from the deduplicated BAM, using the `--intervals` exome capture regions to build target/antitarget bins. It runs in germline/tumor-only mode (no matched normal, no panel of normals) - CNVkit builds a "flat" reference internally that assumes uniform coverage, which is the standard approach for single-sample or small-cohort exome CNV calling without a cohort of normal samples to build a proper panel-of-normals reference from. This subworkflow requires `--intervals` to be set; it errors out clearly if it isn't.

### QC metrics table

<details markdown="1">
<summary>Output files</summary>

- **`qc_metrics/`** (headline output)
  - `*.qc_metrics.tsv`: flat, tabular QC summary - one row per sample, columns from FastQC, samtools flagstat/stats and GATK MarkDuplicates (MultiQC's "general stats" table, republished as a plain TSV).

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline (this is where `qc_metrics/*.qc_metrics.tsv` is copied from).
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC, samtools and GATK. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.
  - Collated software versions: `varcall_software_mqc_versions.yml`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
