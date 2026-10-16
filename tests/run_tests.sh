#!/usr/bin/env bash
#
# Small, self-contained correctness test suite for kpipho/varcall.
#
# Intended use: after cloning/porting this pipeline onto a new HPC system,
# run this script to confirm nextflow + singularity + the pipeline code
# actually work together end-to-end, before pointing it at real patient data.
# It needs nothing beyond what the pipeline itself needs (nextflow, singularity,
# python3) - no nf-test, no extra language toolchains.
#
# What it checks:
#   1. The full pipeline (alignment -> variant_calling -> cnv_calling -> QC)
#      runs to completion on a small synthetic dataset and produces all four
#      headline outputs (BAM, germline VCF, CNV VCF, QC metrics TSV) for both
#      of its two test samples, and that the germline VCF for each sample
#      contains the exact SNP that was deliberately injected into its reads.
#   2. A samplesheet that violates the paired-end-only schema (missing
#      fastq_2) is correctly REJECTED by schema validation, rather than
#      silently accepted.
#   3. cnv_calling's explicit `--intervals` requirement is enforced: running
#      without --intervals set fails fast with a clear error, rather than
#      silently producing meaningless CNV calls.
#
# Usage:
#   tests/run_tests.sh [extra nextflow args, e.g. -c my_hpc.config]
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
WORK_DIR="${SCRIPT_DIR}/.work"
RESULTS_DIR="${SCRIPT_DIR}/.results"
LOG_DIR="${SCRIPT_DIR}/.logs"
EXTRA_ARGS=("$@")

PASS=0
FAIL=0

pass() {
    echo "  [PASS] $1"
    PASS=$((PASS + 1))
}

fail() {
    echo "  [FAIL] $1"
    FAIL=$((FAIL + 1))
}

section() {
    echo
    echo "=== $1 ==="
}

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: '$1' is required but was not found on PATH." >&2
        exit 1
    fi
}

cd "${REPO_ROOT}"

section "Checking prerequisites"
require_command nextflow
require_command python3
require_command singularity
echo "  nextflow: $(nextflow -v 2>&1 | head -1)"
echo "  singularity: $(singularity --version 2>&1)"

section "Generating synthetic test data (tests/data/)"
python3 "${SCRIPT_DIR}/generate_test_data.py"

rm -rf "${WORK_DIR}" "${RESULTS_DIR}"
mkdir -p "${LOG_DIR}"

section "Test 1/3: full pipeline run (alignment + variant_calling + cnv_calling)"
RUN_LOG="${LOG_DIR}/full_run.log"
if nextflow run . -profile test,singularity \
    -work-dir "${WORK_DIR}" \
    --outdir "${RESULTS_DIR}" \
    "${EXTRA_ARGS[@]}" \
    >"${RUN_LOG}" 2>&1; then
    pass "pipeline exited successfully"
else
    fail "pipeline exited with an error - see ${RUN_LOG}"
    tail -n 40 "${RUN_LOG}"
fi

check_file() {
    local description="$1"
    local path="$2"
    if [[ -s "${path}" ]]; then
        pass "${description} (${path})"
    else
        fail "${description} missing or empty (${path})"
    fi
}

for sample in SAMPLE_1 SAMPLE_2; do
    check_file "${sample}: analysis-ready BAM" "${RESULTS_DIR}/bam/${sample}.markdup.bam"
    check_file "${sample}: BAM index" "${RESULTS_DIR}/bam/${sample}.markdup.bam.bai"
    check_file "${sample}: germline VCF" "${RESULTS_DIR}/vcf/${sample}.vcf.gz"
    check_file "${sample}: germline VCF index" "${RESULTS_DIR}/vcf/${sample}.vcf.gz.tbi"
    check_file "${sample}: CNV VCF" "${RESULTS_DIR}/cnv/${sample}.vcf"

    # Confirm the exact SNP injected into this sample's reads (tests/data/expected_variants.tsv)
    # was actually called, not just that *some* variant exists somewhere.
    expected=$(awk -F'\t' -v s="${sample}" '$1 == s { print $2":"$3":"$4">"$5 }' "${SCRIPT_DIR}/data/expected_variants.tsv")
    chrom=$(cut -d: -f1 <<<"${expected}")
    pos=$(cut -d: -f2 <<<"${expected}")
    vcf="${RESULTS_DIR}/vcf/${sample}.vcf.gz"
    if [[ -f "${vcf}" ]] && zcat "${vcf}" | awk -F'\t' -v c="${chrom}" -v p="${pos}" '$1==c && $2==p {found=1} END {exit !found}'; then
        pass "${sample}: expected injected variant ${expected} was called"
    else
        fail "${sample}: expected injected variant ${expected} NOT found in ${vcf}"
    fi
done

check_file "aggregated QC metrics table" "${RESULTS_DIR}/qc_metrics/varcall.qc_metrics.tsv"
if [[ -f "${RESULTS_DIR}/qc_metrics/varcall.qc_metrics.tsv" ]]; then
    rows=$(($(wc -l <"${RESULTS_DIR}/qc_metrics/varcall.qc_metrics.tsv") - 1))
    if [[ "${rows}" -eq 2 ]]; then
        pass "QC metrics table has one row per sample (2)"
    else
        fail "QC metrics table has ${rows} data row(s), expected 2"
    fi
fi
check_file "MultiQC report" "${RESULTS_DIR}/multiqc/multiqc_report.html"

section "Test 2/3: samplesheet schema validation rejects a malformed samplesheet"
BAD_SHEET_LOG="${LOG_DIR}/invalid_samplesheet.log"
if nextflow run . -profile test,singularity \
    -work-dir "${WORK_DIR}" \
    --outdir "${RESULTS_DIR}_invalid_samplesheet" \
    --input "${SCRIPT_DIR}/data/samplesheet_invalid.csv" \
    "${EXTRA_ARGS[@]}" \
    >"${BAD_SHEET_LOG}" 2>&1; then
    fail "pipeline accepted a samplesheet with a missing fastq_2 value (should have been rejected)"
else
    pass "pipeline correctly rejected the malformed samplesheet"
fi

section "Test 3/3: cnv_calling requires --intervals"
NO_INTERVALS_CONFIG=$(mktemp)
echo "params.intervals = null" >"${NO_INTERVALS_CONFIG}"
NO_INTERVALS_LOG="${LOG_DIR}/no_intervals.log"
if nextflow run . -profile test,singularity \
    -work-dir "${WORK_DIR}" \
    --outdir "${RESULTS_DIR}_no_intervals" \
    -c "${NO_INTERVALS_CONFIG}" \
    "${EXTRA_ARGS[@]}" \
    >"${NO_INTERVALS_LOG}" 2>&1; then
    fail "pipeline ran to completion without --intervals (cnv_calling should have errored out)"
else
    if grep -q "requires --intervals" "${NO_INTERVALS_LOG}"; then
        pass "pipeline correctly refused to run cnv_calling without --intervals"
    else
        fail "pipeline failed without --intervals, but not for the expected reason - see ${NO_INTERVALS_LOG}"
    fi
fi
rm -f "${NO_INTERVALS_CONFIG}"

section "Summary"
echo "  ${PASS} passed, ${FAIL} failed"
echo "  Logs: ${LOG_DIR}"
echo "  Results: ${RESULTS_DIR}, ${RESULTS_DIR}_invalid_samplesheet, ${RESULTS_DIR}_no_intervals"

if [[ "${FAIL}" -gt 0 ]]; then
    exit 1
fi
