# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An nf-core-style [Nextflow](https://www.nextflow.io/) DSL2 pipeline (`kpipho/varcall`) for human **whole-exome sequencing (WES)** data against **hg38**: it aligns paired-end Illumina reads with BWA-MEM, calls germline short variants with GATK, and calls copy-number variants with CNVkit. Scaffolded with `nf-core pipelines create` (nf-core/tools v4.1.0) and built out from there — it is a standalone pipeline, not part of the nf-core organisation (`is_nfcore: false` in `.nf-core.yml`). Resource requests throughout `conf/base.config`/`conf/modules.config` are sized for this specific case (exome data volume, hg38-scale reference) — see "Resource sizing" below before reusing them for WGS or a different organism.

Headline outputs, one set per sample, in `--outdir`: `bam/<sample>.markdup.bam(.bai)`, `vcf/<sample>.vcf.gz(.tbi)`, `cnv/<sample>.vcf`, `qc_metrics/<sample>.qc_metrics.tsv`. See `docs/output.md` for the full layout.

## Environment

Nextflow and Java are not on `PATH` by default in a fresh shell in this sandbox; a working toolchain (OpenJDK 17, Nextflow, an `nf-core` Python venv, and SRA Toolkit) was set up once under `~/.local/opt`, `~/.local/bin`, and `/tmp/nf-core-venv`. Source `~/.nf_env` (if present) to put them on `PATH`, or otherwise ensure `nextflow`, `nf-core`, and a JDK are installed before running any commands below. Singularity (`singularity-ce`) is available system-wide.

## Commands

```bash
# Run the pipeline (HPC/Singularity is the primary target profile). --intervals (exome
# capture BED) is optional for alignment/variant_calling but required for cnv_calling.
nextflow run . -profile singularity --input samplesheet.csv --fasta genome.fa --intervals capture.bed --outdir results

# Lint (nf-core/tools, not plain `nextflow lint`)
nf-core pipelines lint

# Install/remove/update an nf-core module or subworkflow (updates modules.json + conf/modules.config entries)
nf-core modules install <tool/subtool>
nf-core modules remove <tool/subtool>
```

```bash
# Fast, self-contained correctness test suite (synthetic 2-sample data, generated on
# the fly - needs only nextflow/singularity/python3, no nf-test). Run this after any
# change to subworkflow wiring - it's what caught the channel-broadcast bug below.
tests/run_tests.sh
```

There is no nf-test suite (the `nf-test`/`ci` template features were deliberately
skipped — see `skip_features` in `.nf-core.yml`). Instead there's `tests/` (see
`tests/README.md`): `tests/generate_test_data.py` deterministically generates a tiny
two-sample synthetic dataset (with known injected variants) and `tests/run_tests.sh`
runs the full pipeline against it via `-profile test,singularity` (`conf/test.config`)
plus two negative tests (bad samplesheet, missing `--intervals`). This has also been
validated once against real data: hg38 **chr22 only** + a 500k-read-pair subsample of a
public SRA exome run (`SRR30504930`, pulled via `prefetch`/`fasterq-dump`), with a
`conf/modules.config`-shaped resource override scaled down to fit a small dev machine.
That dataset was built under `/tmp` and is not part of the repo - re-derive a similar
one if you need to repeat that specific check; for everyday wiring changes, prefer
`tests/run_tests.sh` since it's fast and self-contained.

## Architecture

```
main.nf                                    entrypoint; calls workflows/varcall.nf
workflows/varcall.nf                       wires ALIGNMENT -> VARIANT_CALLING -> CNV_CALLING -> MultiQC -> QC_METRICS_TABLE
subworkflows/local/alignment/               BWA_INDEX, BWA_MEM, SAMTOOLS_FAIDX/INDEX/FLAGSTAT/STATS, FASTQC
subworkflows/local/variant_calling/         GATK4_CREATESEQUENCEDICTIONARY, GATK4_MARKDUPLICATES, SAMTOOLS_INDEX, GATK4_HAPLOTYPECALLER
subworkflows/local/cnv_calling/              CNVKIT_BATCH (germline/tumor-only, flat reference), CNVKIT_EXPORT (-> VCF)
subworkflows/local/utils_nfcore_varcall_pipeline/   stock nf-core template plumbing (samplesheet parsing, help text, versions)
subworkflows/nf-core/, modules/nf-core/     unmodified nf-core modules/subworkflows (installed via `nf-core modules install`)
modules/local/qc_metrics_table/             the one hand-written process: republishes MultiQC's general-stats table as a flat TSV
conf/base.config                            per-label (process_low/medium/high/...) default resources
conf/modules.config                         per-process publishDir layout + ext.args/ext.prefix/resource overrides
assets/schema_input.json                    samplesheet JSON schema (paired-end only; 5 required columns)
nextflow_schema.json                        `--fasta`/`--bwa_index`/`--fasta_fai`/`--dict`/`--intervals` live under `reference_genome_options`
```

Both local subworkflows that need the reference genome take it as a `[meta, path]` value channel (built once in `workflows/varcall.nf`, reused for every sample) and read `params.bwa_index` / `params.fasta_fai` / `params.dict` / `params.intervals` directly to decide whether to build each reference artifact or reuse a pre-built one — see the `if (params.xxx) {...} else { BUILD_IT(...) }` pattern in each subworkflow's `main.nf`. `ALIGNMENT` emits the `.fai` index so `VARIANT_CALLING` doesn't rebuild it. `CNV_CALLING` hard-requires `params.intervals` (`error(...)` if unset) since CNVkit's target/antitarget binning is meaningless without capture regions; `VARIANT_CALLING` treats it as optional (restricts HaplotypeCaller + adds 100bp padding when set, scans the whole genome when not).

## Resource sizing

`conf/base.config`'s labels and the per-process overrides in `conf/modules.config` (`BWA_INDEX`, `GATK4_MARKDUPLICATES`, `GATK4_HAPLOTYPECALLER`, `CNVKIT_BATCH`, `MULTIQC`) are deliberately tuned for _exome data against the full hg38 reference_, not WGS and not a smaller/larger genome. The two things pulling in opposite directions (don't "fix" one without the other): exome per-sample read volume is small, but BWA_MEM/CNVKIT_BATCH still load or scan the whole reference regardless of how little data is aligned to it, so their memory floor is set by genome size. If you repurpose this pipeline for WGS, revisit `process_medium`/`process_high` and those explicit overrides — they scale with depth across the whole genome, not just the exome. The rationale is written out in comments at the top of both config files; don't delete those comments when editing the numbers.

When testing on a dev machine instead of an HPC cluster, don't edit `conf/base.config`/`conf/modules.config` in place — pass a `-c /path/to/local_resources.config` override on the command line instead (same `withLabel`/`withName` selectors, smaller numbers) so the real HPC-sized defaults stay intact in the repo.

## Non-obvious gotchas

- **`GATK4_MARKDUPLICATES` needs `ext.prefix` set to something other than `"${meta.id}.bam"`.** The upstream aligned BAM is also named `"${meta.id}.bam"`; if MarkDuplicates' output prefix collides with its input filename, GATK silently writes through the staged-input symlink instead of creating a new file in its own work dir, and Nextflow (correctly) refuses to re-collect an input file as a process output — the `bam` emit channel just comes up empty, with no error, and everything downstream (indexing, HaplotypeCaller) silently never runs. Fixed via `ext.prefix = { "${meta.id}.markdup.bam" }` in `conf/modules.config`; don't remove it.
- **`BWA_MEM` sets `ext.args` to add a read group (`-R '@RG\tID:...'`).** GATK MarkDuplicates/HaplotypeCaller require one; BWA doesn't add it by default.
- **Never pass a `.join()` of two value ("single-item") channels straight into a process as a broadcast reference input - it silently runs the process for only the first sample, with no error.** Confirmed with a minimal standalone repro (no modules, no stub-run involved) on Nextflow 26.10.x: `ch_fasta.join(ch_fasta_fai)` passed directly as a second process input alongside a 2-item per-sample channel executed the process for sample 1 only; sample 2's invocation never appeared in the trace at all. This bit both `ALIGNMENT`'s `SAMTOOLS_STATS` call and `CNV_CALLING`'s `CNVKIT_BATCH` call during testing (`tests/run_tests.sh`, which uses 2 samples, is what caught it - a single-sample test would not have). The fix in both places: `ch_bam_bai.combine(ch_fasta.join(ch_fasta_fai)).map { ... }`, called once and then `.map()`'d again per process-input-shape needed - `.combine()` is a genuine cross-product and does not have this problem. Passing an _unmodified_ `channel.value(...)` (not run through `.join()` or another operator) directly as a broadcast input is fine and used correctly elsewhere (e.g. `ch_fasta` itself into `BWA_MEM`) - the hazard is specifically operator-chained/derived value channels, not value channels in general. If you add a new reference-broadcast input anywhere, default to the `.combine()` pattern rather than assuming implicit broadcasting will work - it does not have a proven-reliable heuristic here, so verify empirically if you deviate.
- `conf/containers_*.config` files (per-platform container pinning, auto-generated by `nf-core modules install/update`) are intentionally deleted after every module change — they aren't `includeConfig`'d anywhere (every module already resolves its own Singularity/Docker/conda container inline), so they're just regenerate-and-redelete noise. `container_configs` is disabled in `.nf-core.yml`'s lint config accordingly.
- **`CNVKIT_BATCH.out.cns` bundles two different files under one glob.** `cnvkit batch` (run via `CNVKIT_BATCH`) writes both `<id>.cns` (raw segments) and `<id>.call.cns` (integer copy-number calls) into the same task work dir, and the module's `path("*.cns")` output matches both, so that emit channel item is a `[meta, [file1, file2]]` list, not a single file. `subworkflows/local/cnv_calling/main.nf` filters for the `.call.cns` one by filename before passing it to `CNVKIT_EXPORT` — don't swap in a plain `.map { meta, cns -> ... }` that assumes a single bare file.
