# kpipho/varcall: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

Initial release of kpipho/varcall, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `alignment` subworkflow: BWA-MEM alignment to a reference genome, with FastQC + samtools flagstat/stats alignment QC
- `variant_calling` subworkflow: GATK MarkDuplicates and HaplotypeCaller germline short-variant calling
- Flat, tabular QC metrics file (one row per sample) alongside the final BAM and VCF

### `Fixed`

### `Dependencies`

### `Deprecated`
