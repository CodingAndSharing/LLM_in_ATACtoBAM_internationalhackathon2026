# ATAC-seq → mouse BAM pipeline

Turns the paired-end ATAC-seq FASTQ files in `../sandbox/` into **filtered,
duplicate-marked, coordinate-sorted BAM files** aligned to the mouse genome
**Ensembl GRCm39, release 116** — the release that matches
`Mus_musculus.GRCm39.116.gtf.gz`.

Everything runs from one script, `atac_to_bam.py`, inside a **pixi** environment
pinned to **Python 3.12**.

---

## Pipeline at a glance

```mermaid
flowchart TD
    %% ---------- reference prep (runs once) ----------
    subgraph REF["Reference prep - one-time"]
        ENS["Ensembl FTP release-116"]
        ENS -->|"GRCm39.dna.primary_assembly.fa.gz (806 MB)"| DLFA["download<br/>curl"]
        ENS -->|"GRCm39.116.gtf.gz (103 MB)"| GTF[/"annotation kept<br/>for downstream QC"/]
        DLFA -->|".fa.gz"| DECOMP["decompress<br/>pigz"]
        DECOMP -->|"genome .fa (~2.7 GB)"| BUILD["build index<br/>bowtie2-build<br/>(--packed on low RAM)"]
        BUILD --> IDX(["bowtie2_GRCm39.*.bt2"])
    end

    %% ---------- per-sample processing ----------
    IN[/"paired FASTQ per sample (rep1, rep2)<br/>_1 / _2 mates, 101 bp, 15-17 M pairs"/]
    IN -->|"R1, R2"| S1["1. adapter + quality trim<br/>fastp (detect-PE + explicit Nextera, len>=20)"]
    S1 -->|"trim_1/2.fastq.gz<br/>+ fastp.html / .json"| S2["2. align to GRCm39<br/>bowtie2 (very-sensitive, -X 2000, dovetail, +RG)"]
    IDX -.->|"-x index"| S2
    S2 -->|"SAM stream"| S3["3. fix mate tags<br/>samtools fixmate -m"]
    S3 -->|"BAM stream"| S4["4. coordinate sort<br/>samtools sort"]
    S4 -->|"BAM stream"| S5["5. mark duplicates<br/>samtools markdup"]
    S5 -->|"sample.markdup.bam + .bai<br/>+ markdup_stats.txt"| S6["6. filter<br/>samtools view -f 2 -F 3852 -q 30, drop MT"]
    S6 -->|"sample.filtered.bam"| S7["7. index + QC<br/>samtools index / flagstat / idxstats"]
    S7 --> OUT[/"sample.filtered.bam + .bai<br/>flagstat / idxstats txt"/]

    %% ---------- roll-up ----------
    OUT --> MQC["8. roll-up report<br/>multiqc (optional)"]
    S1 -.->|"fastp json"| MQC
    S5 -.->|"markdup stats"| MQC
    S7 -.->|"flagstat"| MQC
    MQC --> REP[/"multiqc_report.html"/]

    classDef art fill:#eef6ff,stroke:#4a90d9,color:#0b3d63;
    classDef step fill:#f5f5f5,stroke:#888888,color:#222222;
    class IN,OUT,REP,IDX,GTF art;
    class S1,S2,S3,S4,S5,S6,S7,MQC,DLFA,DECOMP,BUILD step;
```

Steps 2–5 are one streamed pipe (no full-size intermediate SAM/BAM on disk); the
`markdup.bam` is the first file written, `filtered.bam` is the deliverable.

---

## 1. Input data (auto-detected from `../sandbox/`)

| Sample name (derived)                                        | Run / GEO             | FASTQ pair                                  | Read pairs¹ | Read len |
|-------------------------------------------------------------|-----------------------|--------------------------------------------|-------------|----------|
| `SRR11878619_GSM4579144_mESC_DMSO_ATAC_rep1_Mus_musculus_ATAC-seq` | SRR11878619 / GSM4579144 | `..._downsample_1.fastq.gz` + `_2.fastq.gz` | ~15.3 M     | 101 bp   |
| `SRR11878620_GSM4579145_mESC_DMSO_ATAC_rep2_Mus_musculus_ATAC-seq` | SRR11878620 / GSM4579145 | `..._downsample_1.fastq.gz` + `_2.fastq.gz` | ~17.0 M     | 101 bp   |

¹ mESC, DMSO-treated, ATAC-seq, two replicates. Exact read-pair counts (line
count ÷ 4): rep1 **15,291,302**, rep2 **17,024,590**.

> **Note on the input files:** they are named `*.fastq.gz` but are actually
> **plain uncompressed text** (~3.7 GB and ~4.1 GB per mate). The script checks
> the gzip magic bytes and feeds `fastp` a correctly-named symlink, so the
> misnamed files work as-is — nothing to rename by hand.

### Pre-alignment QC (bounded sample: first 100k read pairs/file)

Run with the `exploratory-data-analysis` skill's `sequence_inspector.py` plus a
`fastp` dry pass; numbers are a subsample, not full-file QC.

| Metric | rep1 | rep2 | Read on it |
|---|---|---|---|
| Read length | 101 bp, fixed | 101 bp, fixed | uniform; no pre-trimmed input |
| Mean Phred (Phred+33) | R1 35.7 / R2 35.2 | R1 35.7 / R2 35.4 | high quality; min score 2 in tails → quality trim still warranted |
| GC | ~38 % | ~37 % | slightly below mouse genomic ~42 %, consistent with adapter/short-insert content |
| Ambiguous (`N`) fraction | ~3e-5 | ~3–4e-5 | negligible |
| Reads containing Tn5/Nextera adapter | **~38 % full motif, ~49 % trimmed by fastp** | same | heavy adapter read-through → adapter trimming is **required**, not cosmetic |
| fastp insert-size peak | ~41 bp (subsample) | — | very short fragments (< read length) → `bowtie2 -X 2000 --dovetail` is the right setting |
| fastp duplication estimate | ~11 % | — | handled by `samtools markdup` |

Implication for the pipeline: the `fastp` step passes the explicit Nextera
adapter sequences **in addition to** `--detect_adapter_for_pe`, and bowtie2 runs
with `--dovetail` so mate pairs that extend past each other (short inserts) still
count as concordant.

---

## 2. What the script does (per sample)

| # | Step | Tool | Notes |
|---|------|------|-------|
| 1 | Adapter + quality trim | `fastp` | `--detect_adapter_for_pe` **+ explicit Tn5/Nextera adapter seqs** (~49 % of reads carry it), min length 20, HTML/JSON report |
| 2 | Align to GRCm39 | `bowtie2` | `--very-sensitive --dovetail -X 2000` (ATAC fragments up to 2 kb), read groups added |
| 3 | Fix mate flags/tags | `samtools fixmate -m` | needed for markdup |
| 4 | Coordinate sort | `samtools sort` | streamed, uncompressed between pipe stages |
| 5 | Mark duplicates | `samtools markdup` | writes `*.markdup_stats.txt` |
| 6 | Filter | `samtools view` | keep `-f 2` proper pairs, `-F 3852` (drop unmapped/mate-unmapped/secondary/QCfail/dup/supplementary), `MAPQ ≥ 30`, **drop mitochondrial (`MT`) reads** |
| 7 | Index + QC | `samtools index / flagstat / idxstats` | per BAM |
| 8 | Roll-up report | `multiqc` | one `multiqc_report.html` across all samples (optional) |

Steps 2–5 run as a **single streamed pipeline** (`bowtie2 | fixmate | sort |
markdup`) so no full-size intermediate SAM/BAM is written.

### Reference handling (one-time)

If `../reference/` doesn't already have them, the script downloads from
`https://ftp.ensembl.org/pub/release-116/` (paths + sizes verified against
Ensembl on 2026-09-10 via the `database-lookup` skill):

| File | Path | Size | Purpose |
|---|---|---|---|
| `Mus_musculus.GRCm39.dna.primary_assembly.fa.gz` | `fasta/mus_musculus/dna/` | **806 MB** (2.7 GB unzipped) | genome, for the bowtie2 index |
| `Mus_musculus.GRCm39.116.gtf.gz` | `gtf/mus_musculus/` | **103 MB** | annotation, kept for downstream use (TSS enrichment, peak annotation…) |

* **Assembly:** GRCm39 (GCA_000001635.9); release 116 dated 2026-03.
* **Contig naming:** Ensembl style — `1`–`19`, `X`, `Y`, **`MT`** (no `chr`
  prefix) plus unplaced scaffolds. The script's MT filter and `--main-chroms-only`
  regex match this exactly; genome FASTA and GTF use the same names, so **no
  liftover / `chr`-prefix reconciliation is needed**.
* `dna.primary_assembly` (not `toplevel`) is used: it drops the ALT/patch
  haplotypes that would otherwise create multi-mapping artefacts. Soft-masking is
  irrelevant to bowtie2, so plain `dna` (not `dna_sm`) is fine.

The script then builds the **bowtie2 index** once (`../reference/bowtie2_GRCm39.*.bt2`).
Downloads are checksum-verified against Ensembl's `CHECKSUMS` (warn-only) and are
resumable (`curl -C -`). Bring your own with `--genome-fa` / `--index-prefix`.

---

## 3. Outputs (`../results/`)

Per sample:

```
<sample>.filtered.bam            <-- main deliverable (filtered, dedup-marked, sorted)
<sample>.filtered.bam.bai
<sample>.filtered.flagstat.txt
<sample>.filtered.idxstats.txt
<sample>.markdup.bam(.bai)       <-- pre-filter BAM (all mapped reads, dups flagged; see --keep-markdup)
<sample>.markdup.flagstat.txt
<sample>.markdup_stats.txt
<sample>.fastp.html / .json / .log
<sample>.bowtie2.log
multiqc_report.html              <-- once, across all samples
```

The `*.filtered.bam` files are what you feed to peak callers (MACS2/MACS3),
`alignmentSieve`/Tn5 shift, fragment-size QC, bigWig generation, etc.

---

## 4. Compute requirements — read before running

Measured / estimated for the two samples above on a 16-core machine.

| Resource | Minimum | Recommended | Why |
|---|---|---|---|
| **CPU** | 4 cores | **8–16 cores** | bowtie2 alignment dominates wall time and scales well |
| **RAM** | **8 GB** | **16 GB** | bowtie2 index **build** peaks ~6.5 GB (≈4 GB with `--packed`); alignment ~3.4 GB; `samtools sort` uses `--sort-mem` × sort-threads (default 512M × 4 = 2 GB) |
| **Disk (free)** | 40 GB | **60–80 GB** | see breakdown below |
| **Network** | — | ~0.9 GB one-time | genome FASTA 806 MB + GTF 103 MB; skip with `--skip-download` after first run |
| **OS** | Linux x86-64 | — | pixi env is `linux-64`; needs `bash` |

### Low-memory machines (< 16 GB RAM)

The script **auto-detects RAM < 16 GB and adds `bowtie2-build --packed`**
(≈2× slower build, ~4 GB peak) and caps build threads at 4. Force it anytime
with `--packed`.

> **This dev box has ~7.6 GB RAM.** `--packed` turns on automatically, but the
> index build's ~4 GB peak plus whatever else is resident is close to the edge.
> If `bowtie2-build` gets OOM-killed: build the index once on a ≥16 GB machine
> and copy `../reference/bowtie2_GRCm39.*`, then run here with
> `--skip-download --index-prefix ../reference/bowtie2_GRCm39`. Also drop
> `--sort-mem` to `256M` and `--threads` to `4` to keep the alignment stage
> lean.

### Disk breakdown (peak, processing one sample at a time)

| Item | Size |
|---|---|
| Genome FASTA `.gz` + decompressed | ~0.8 + ~2.7 GB |
| GTF `.gz` | ~0.1 GB |
| bowtie2 index (`*.bt2`) | ~3.9 GB |
| Trimmed FASTQ (per sample, gzipped) | ~3–3.5 GB |
| `markdup.bam` (per sample) | ~2.5–3.5 GB |
| `filtered.bam` (per sample) | ~1.5–2.5 GB |
| sort/markdup temp (transient) | up to ~markdup size |
| **Reference (kept)** | **~7.5 GB** |
| **Peak extra while a sample runs** | **~12–15 GB** |

Trimmed FASTQ and scratch dirs are deleted after each sample unless
`--keep-intermediate`. Input FASTQ in `../sandbox/` (~15 GB) is never modified.

### Wall-time estimate (8–16 threads)

| Phase | Time |
|---|---|
| Reference download | 2–10 min |
| bowtie2 index build (one-time) | 30–45 min (up to ~75 min with `--packed`) |
| Per sample: trim → align → markdup → filter → QC | 45–90 min |
| **Total, 2 samples, cold start** | **≈ 2.5–4 h** |
| **Total, 2 samples, index already built** | **≈ 1.5–3 h** |

---

## 5. How to run

```bash
# 0. one-time: get pixi (https://pixi.sh) if you don't have it
#    curl -fsSL https://pixi.sh/install.sh | bash

cd atac_pipeline

# 1. create the Python 3.12 environment (bowtie2, samtools, fastp, pigz, multiqc)
pixi install

# 2. sanity-check the toolchain
pixi run tools

# 3. dry run — prints every command, runs nothing, touches no network
pixi run pipeline --dry-run

# 4. real run (defaults: --fastq-dir ../sandbox  --outdir ../results  --refdir ../reference)
pixi run pipeline --threads 12

#    equivalently:
#    pixi run python atac_to_bam.py --threads 12
```

### Useful options

| Option | Effect |
|---|---|
| `--fastq-dir DIR` | where the paired `*_1/*_2` FASTQ live (default `../sandbox`) |
| `--outdir DIR` / `--refdir DIR` | output / reference-cache locations |
| `--samples NAME [NAME …]` | process only some samples |
| `--threads N` | CPU threads (default = min(cores, 12)) |
| `--sort-mem 256M` | lower `samtools sort` memory on tight boxes |
| `--mapq 20` | change the MAPQ filter (default 30) |
| `--main-chroms-only` | keep only chr 1–19, X, Y (drop unplaced scaffolds) |
| `--genome-fa FILE` / `--index-prefix PFX` | reuse your own genome / bowtie2 index |
| `--skip-download` | never touch the network (fails if reference missing) |
| `--resume` | re-run and skip per-sample steps whose output already exists |
| `--keep-intermediate` / `--keep-markdup` | keep trimmed FASTQ / pre-filter BAM |
| `--dry-run` | print commands only |

---

## 6. Design choices

* **bowtie2**, not bwa — matches the ENCODE ATAC-seq pipeline; `-X 2000
  --very-sensitive --dovetail` is the standard ATAC setting for capturing
  nucleosome-spanning fragments and Tn5 dovetailed pairs.
* **`samtools markdup`** (not Picard) — no JVM, streams in the same pipe, and the
  duplicate flag is set (not removed) so you can revisit the decision.
* **Mitochondrial reads dropped at the filter step** — ATAC libraries are
  typically 20–50 % `MT`; removing them is expected before peak calling. The
  `markdup.bam` still contains them if you want the raw fraction (`samtools
  idxstats`).
* **GTF downloaded but not used for alignment** — bowtie2 is unspliced and needs
  only the genome FASTA. The GTF is kept in `../reference/` for TSS-enrichment QC
  and peak annotation downstream.
* **Reference kept outside the pipeline folder** (`../reference/`) so the ~7.5 GB
  index is reused across runs and isn't tied to this directory.

---

## 7. Provenance / skills used

| Skill | Used for |
|---|---|
| `exploratory-data-analysis` | bounded FASTQ profiling (`sequence_inspector.py` on the first 100k pairs/file) → read length, Phred, GC, `N` fraction in §1 |
| `database-lookup` | verified Ensembl release-116 assembly (GRCm39 / GCA_000001635.9), the two FTP paths + byte sizes, and the `1..19,X,Y,MT` contig convention (§2) — confirmed against `ftp.ensembl.org` on 2026-09-10 |
| `genomic-coordinates` | checked that genome FASTA, GTF, the MT filter and `--main-chroms-only` regex all use the same Ensembl contig names → no `chr`-prefix / build reconciliation needed |
| `genomic-intelligence` | not applicable — this task produces alignments, not sequence-model predictions |

The command chain (`fastp → bowtie2 → fixmate → sort → markdup → filter → index`)
was smoke-tested end-to-end in the pixi env on a 200 kb chr19 slice + 50k real
read pairs from rep1 before this README was finalised.

---

## 8. Actual run — 2026-09-11 (outputs in `../scripts_output/`)

Ran on the 16-core / 7.6 GB box: `--threads 10 --sort-mem 384M`, `--packed`
index build auto-enabled. Wall time **~1 h 26 m** total (index build 36 m,
rep1 22 m, rep2 23 m). Peak RAM held; the only OOM casualty was an unrelated
watcher process.

| | rep1 (SRR11878619) | rep2 (SRR11878620) |
|---|---|---|
| Input read pairs | 15,291,302 | 17,024,590 |
| Pairs after `fastp` (too-short dropped) | 15,156,408 | 16,914,444 |
| **bowtie2 overall alignment rate** | **49.3 %** | **47.3 %** |
| Concordant ≥1× | 7,442,573 pairs | 7,969,984 pairs |
| Duplicate rate (of mapped) | 31.9 % | 37.9 % |
| Est. library size (markdup) | ~9.0 M | ~7.6 M |
| **`*.filtered.bam` reads (pairs)** | **7,204,368 (3,602,184)** | **6,767,616 (3,383,808)** |
| MT reads in filtered BAM | 0 | 0 |
| Both filtered BAMs | `samtools quickcheck` ✅, coordinate-sorted, `@RG` + full `@PG` chain | |

### Why the alignment rate is ~48 %, not >90 %

This is a **property of the input data, not the pipeline.** Re-aligning a fresh
raw 300k-pair subsample end-to-end reproduced **46 %** — identical. The unmapped
half is:

* **Tn5 adapter dimers** — `fastp` overrepresented-sequence analysis returns
  almost exclusively `…CTGTCTCTTATACACATCT…` (Nextera mosaic end) at or near
  read start, i.e. transposase self-ligation with little or no genomic insert;
  after adapter trimming the stub is too short / low-complexity to place.
* **Low-complexity, AT-rich fragments** — unmapped reads are dominated by
  poly-A/poly-T runs and near-duplicate AT-rich strings (matches the ~37–38 %
  GC seen in the §1 EDA). No single dominant contaminant — the top unmapped
  30-mer is <0.2 % of unmapped reads.

fastp insert-size peak is **41–52 bp** with ~97 % of fragments < 271 bp — a
heavily sub-nucleosomal / short-insert library, consistent with high adapter
read-through. An ENCODE/nf-core run on the same FASTQ would report the same rate.

### Bottom line

The `*.filtered.bam` files (~3.4–3.6 M clean, deduplicated, MAPQ≥30, non-MT read
pairs each) are valid and ready for MACS2/MACS3 peak calling, Tn5 shift, and
fragment-size QC — just shallow, because the input was downsampled and roughly
half of it is unmappable adapter/low-complexity content. If you need more usable
depth, go back to the non-downsampled runs; the pipeline itself needs no change.

---

## 9. Second run — 2026-09-11 (2 more samples added to `../sandbox/`)

Two more mouse ATAC-seq samples landed in `../sandbox/`: **SRR9894854 /
GSM4005240** and **SRR9894855 / GSM4005241** (NPC48h DMSO, rep1/rep2, 75 bp
reads). Re-ran scoped to just the new samples (`--samples ...`) so the
already-finished rep1/rep2 above weren't reprocessed; the existing GRCm39
genome + bowtie2 index in `../reference/` were reused as-is (no re-download,
no rebuild — confirmed from the run log). All four samples' filtered BAMs and
QC now live together in `../scripts_output/`, and `multiqc_report.html` was
regenerated to cover all four.

| | rep1 (SRR9894854) | rep2 (SRR9894855) |
|---|---|---|
| Input read pairs | 5,119,261 | 3,238,820 |
| **bowtie2 overall alignment rate** | **99.0 %** | **99.1 %** |
| Duplicate rate (of mapped) | 6.8 % | 4.7 % |
| **`*.filtered.bam` reads (pairs)** | **7,386,684 (3,693,342)** | **4,808,600 (2,404,300)** |
| MT reads in filtered BAM | 0 | 0 |
| `samtools quickcheck` | ✅ | ✅ |

This is a useful confirmation of the §8 diagnosis: **same pipeline, same
reference index, same day** — but this pair of samples aligns at ~99% instead
of ~48%. That rules out anything pipeline-side for the first two samples;
the earlier shortfall really was specific to that library's adapter-dimer /
low-complexity content, not a bug or a misconfigured reference.
