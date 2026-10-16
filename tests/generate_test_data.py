#!/usr/bin/env python3
"""
Generates a small, deterministic, synthetic paired-end exome dataset for
testing kpipho/varcall without needing network access or real patient data.

Everything here is made up: a random ~30kb "reference", ten scattered
~200bp "capture target" regions on it, and two samples' worth of paired-end
reads sampled mostly from those regions (plus a bit of off-target noise, like
real exome data). Each sample gets one deliberately injected heterozygous SNP
inside a target region, so the test suite has a concrete, predictable variant
to look for in the output VCF.

Re-run this any time to regenerate tests/data/ (e.g. after editing it) - the
fixed random seed below makes every regeneration produce byte-identical
output, so there is nothing here to commit to git; tests/data/ is gitignored
and (re-)created on demand by tests/run_tests.sh.
"""

import gzip
import random
from pathlib import Path

SEED = 1234
REF_LENGTH = 30_000
CHROM_NAME = "chrTEST"
N_TARGETS = 10
TARGET_LEN = 200
READ_LEN = 100
# Deliberately < 2 * READ_LEN: R1 covers the first READ_LEN bases of the fragment and R2
# (reverse-complemented) covers the last READ_LEN, so together they cover the *whole*
# fragment with some overlap in the middle and no gap - every base in a sequenced fragment
# is covered by at least one read. A longer fragment (more realistic, but with a real
# sequencing gap between R1 and R2) risks landing the injected SNP in that gap, where
# neither read would ever see it.
FRAGMENT_LEN = 180
ON_TARGET_FRAGMENTS = 300
OFF_TARGET_FRAGMENTS = 60

OUT_DIR = Path(__file__).parent / "data"

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def make_reference(rng):
    return "".join(rng.choice(BASES) for _ in range(REF_LENGTH))


def make_targets(rng):
    """Ten non-overlapping ~200bp regions, spread across the reference."""
    span = REF_LENGTH // N_TARGETS
    targets = []
    for i in range(N_TARGETS):
        start = i * span + rng.randint(100, span - TARGET_LEN - 100)
        targets.append((start, start + TARGET_LEN))
    return targets


def inject_snp(ref_seq, targets, rng, avoid_positions):
    """Pick a position inside a target region and substitute one base there."""
    while True:
        target = rng.choice(targets)
        pos = rng.randint(target[0] + 20, target[1] - 20)
        if pos not in avoid_positions:
            break
    ref_base = ref_seq[pos]
    alt_base = rng.choice([b for b in BASES if b != ref_base])
    alt_seq = ref_seq[:pos] + alt_base + ref_seq[pos + 1 :]
    return pos, ref_base, alt_base, alt_seq


def quality_string(length):
    return "I" * length  # flat Phred 40 - good enough for a functional test


def write_fastq_pair(path_r1, path_r2, sample, ref_seq, alt_seq, variant_pos, targets, rng):
    """Writes paired-end reads: mostly on-target (within +/-30bp of a target
    region), a handful off-target, and - for fragments spanning the injected
    SNP - a 50/50 mix of the ref/alt sequence to produce a heterozygous call."""
    r1_records = []
    r2_records = []

    def add_fragment(source_seq, read_index):
        max_start = len(source_seq) - FRAGMENT_LEN - 1
        start = rng.randint(0, max_start)
        fragment = source_seq[start : start + FRAGMENT_LEN]
        read1 = fragment[:READ_LEN]
        read2 = revcomp(fragment[-READ_LEN:])
        name = f"{sample}:sim:{read_index}"
        r1_records.append((name, read1))
        r2_records.append((name, read2))

    read_index = 0
    for _ in range(ON_TARGET_FRAGMENTS):
        target = rng.choice(targets)
        window_start = max(0, target[0] - 30)
        window_end = min(len(ref_seq), target[1] + 30)
        frag_start = rng.randint(window_start, max(window_start, window_end - FRAGMENT_LEN))
        use_alt = (frag_start <= variant_pos < frag_start + FRAGMENT_LEN) and rng.random() < 0.5
        source = alt_seq if use_alt else ref_seq
        fragment = source[frag_start : frag_start + FRAGMENT_LEN]
        read1 = fragment[:READ_LEN]
        read2 = revcomp(fragment[-READ_LEN:])
        name = f"{sample}:sim:{read_index}"
        r1_records.append((name, read1))
        r2_records.append((name, read2))
        read_index += 1

    for _ in range(OFF_TARGET_FRAGMENTS):
        add_fragment(ref_seq, read_index)
        read_index += 1

    with gzip.open(path_r1, "wt") as f1, gzip.open(path_r2, "wt") as f2:
        qual = quality_string(READ_LEN)
        for name, seq in r1_records:
            f1.write(f"@{name}/1\n{seq}\n+\n{qual}\n")
        for name, seq in r2_records:
            f2.write(f"@{name}/2\n{seq}\n+\n{qual}\n")


def main():
    rng = random.Random(SEED)

    ref_dir = OUT_DIR / "reference"
    fastq_dir = OUT_DIR / "fastq"
    ref_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir.mkdir(parents=True, exist_ok=True)

    ref_seq = make_reference(rng)
    targets = make_targets(rng)

    with open(ref_dir / "genome.fa", "w") as f:
        f.write(f">{CHROM_NAME}\n")
        for i in range(0, len(ref_seq), 70):
            f.write(ref_seq[i : i + 70] + "\n")

    with open(ref_dir / "capture_targets.bed", "w") as f:
        for start, end in targets:
            f.write(f"{CHROM_NAME}\t{start}\t{end}\n")

    samples = {}
    used_positions = set()
    for sample in ("SAMPLE_1", "SAMPLE_2"):
        pos, ref_base, alt_base, alt_seq = inject_snp(ref_seq, targets, rng, used_positions)
        used_positions.add(pos)
        samples[sample] = {"pos": pos, "ref": ref_base, "alt": alt_base}
        write_fastq_pair(
            fastq_dir / f"{sample}_R1.fastq.gz",
            fastq_dir / f"{sample}_R2.fastq.gz",
            sample,
            ref_seq,
            alt_seq,
            pos,
            targets,
            rng,
        )

    # 1-based coordinate for human-readable reporting (VCF/BED above stay 0-based/BED convention).
    with open(OUT_DIR / "expected_variants.tsv", "w") as f:
        f.write("sample\tchrom\tpos_1based\tref\talt\n")
        for sample, v in samples.items():
            f.write(f"{sample}\t{CHROM_NAME}\t{v['pos'] + 1}\t{v['ref']}\t{v['alt']}\n")

    def samplesheet_row(sample, fastq_2_path):
        run_date = "2026-01-15"
        sample_type = "tumor" if sample == "SAMPLE_1" else "normal"
        return f"{sample},{fastq_dir / f'{sample}_R1.fastq.gz'},{fastq_2_path},{run_date},{sample_type}\n"

    header = "sample,fastq_1,fastq_2,run_date,sample_type\n"

    with open(OUT_DIR / "samplesheet.csv", "w") as f:
        f.write(header)
        for sample in samples:
            f.write(samplesheet_row(sample, fastq_dir / f"{sample}_R2.fastq.gz"))

    # Negative-test samplesheet: SAMPLE_2's fastq_2 is blank, violating the
    # "paired-end only" schema requirement - used by run_tests.sh to confirm
    # schema validation actually rejects a bad samplesheet instead of silently
    # accepting it.
    with open(OUT_DIR / "samplesheet_invalid.csv", "w") as f:
        f.write(header)
        f.write(samplesheet_row("SAMPLE_1", fastq_dir / "SAMPLE_1_R2.fastq.gz"))
        f.write(samplesheet_row("SAMPLE_2", ""))

    print(f"Generated test data under {OUT_DIR}")
    for sample, v in samples.items():
        print(f"  {sample}: injected SNP {CHROM_NAME}:{v['pos'] + 1} {v['ref']}>{v['alt']}")


if __name__ == "__main__":
    main()
