# Tests

A small, self-contained correctness test suite for kpipho/varcall - intended to be the
first thing you run after porting this pipeline onto a new HPC system, to confirm
Nextflow + Singularity + the pipeline code actually work together end-to-end before
pointing it at real patient data. It needs nothing beyond what the pipeline itself
needs: `nextflow`, `singularity`, and `python3` (standard library only - no extra
packages, no nf-test, no other language toolchain).

## Running

```bash
tests/run_tests.sh
```

Any extra arguments are passed straight through to every `nextflow run` call, so you can
add HPC-specific config, e.g.:

```bash
tests/run_tests.sh -c /path/to/your/institution.config
```

This takes a few minutes (it pulls containers on first run, same as any real pipeline
run) and prints a `[PASS]`/`[FAIL]` line per check, ending in a summary. It exits
non-zero if anything failed, so it's safe to use as a CI/deployment gate too
(`tests/run_tests.sh && echo "ready for production data"`).

## What it checks

1. **Full pipeline run** - `alignment` → `variant_calling` → `cnv_calling` → QC, against
   a tiny synthetic two-sample dataset (generated fresh by `generate_test_data.py`, see
   below). Confirms all four headline outputs exist for both samples, and - this is the
   part that actually proves variant calling is working, not just "the pipeline didn't
   crash" - that each sample's germline VCF contains the _exact_ SNP that was
   deliberately injected into that sample's synthetic reads.
2. **Samplesheet validation** - runs the pipeline against a samplesheet with a blank
   `fastq_2` value (violates the paired-end-only schema) and confirms the pipeline
   _rejects_ it. This is as important to test as the happy path: a schema check that
   silently stops validating is worse than no schema check, because it fails open.
3. **`cnv_calling` requires `--intervals`** - confirms that running without
   `--intervals` fails fast with the expected error, rather than silently producing
   meaningless CNV calls.

## Test data

`generate_test_data.py` deterministically generates (fixed random seed - re-running it
always produces byte-identical output) everything under `tests/data/` (gitignored,
regenerated on demand, never committed):

- `reference/genome.fa` - a ~30kb random single-contig "reference genome"
- `reference/capture_targets.bed` - ten ~200bp scattered "exome capture" regions on it
- `fastq/<SAMPLE>_R{1,2}.fastq.gz` - paired-end reads for two samples (`SAMPLE_1`,
  `SAMPLE_2`), mostly sampled from the capture regions plus a bit of off-target noise
  (like real exome data), each with one deliberately injected heterozygous SNP inside a
  capture region
- `expected_variants.tsv` - the exact position/ref/alt of each sample's injected SNP,
  used by `run_tests.sh` to check the called VCF against a known ground truth
- `samplesheet.csv` / `samplesheet_invalid.csv` - the valid and intentionally-broken
  samplesheets used by tests 1 and 2 above

Run it on its own (e.g. while developing) with:

```bash
python3 tests/generate_test_data.py
```

## Running the test profile manually

`conf/test.config` (the `-profile test` profile) just points `--input`/`--fasta`/
`--intervals` at `tests/data/` and caps every process's resources low enough for a
laptop or a small test allocation - the real HPC-sized defaults in `conf/base.config`/
`conf/modules.config` would be wasteful for a 30kb test genome. If you want to run just
the happy-path case yourself, without the full `run_tests.sh` script:

```bash
python3 tests/generate_test_data.py
nextflow run . -profile test,singularity --outdir test_results
```

## Extending

Adding a new tool to one of the subworkflows? Add an assertion for its output file to
the relevant loop/section in `run_tests.sh` - the existing `check_file` helper covers
the common "file exists and is non-empty" case. If a new check needs different test
data (e.g. a third sample, a different capture layout), extend
`generate_test_data.py` rather than hand-crafting new fixtures, so everything stays
deterministic and regenerable.
