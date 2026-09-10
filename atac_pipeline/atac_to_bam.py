#!/usr/bin/env python3
"""
atac_to_bam.py
==============

Turn paired-end ATAC-seq FASTQ files into analysis-ready BAM files aligned to the
mouse genome (Ensembl GRCm39, release 116 -- the same release as
``Mus_musculus.GRCm39.116.gtf.gz``).

Pipeline, per sample
--------------------
  1. adapter + quality trimming .................. fastp (auto Nextera detection)
  2. alignment .................................. bowtie2  --very-sensitive -X 2000 --dovetail
  3. mate fixup ................................. samtools fixmate -m
  4. coordinate sort ........................... samtools sort
  5. duplicate marking ......................... samtools markdup
  6. filtering ................................. properly-paired, primary, MAPQ>=30,
                                                 no dup/QCfail/secondary/supplementary,
                                                 mitochondrial (MT) reads removed
  7. index + QC stats .......................... samtools index / flagstat / idxstats
  (8. MultiQC report across all samples, if multiqc is installed)

Reference handling
------------------
If not already present it downloads, from https://ftp.ensembl.org/pub/release-116 :
  * Mus_musculus.GRCm39.dna.primary_assembly.fa.gz   (genome, for the bowtie2 index)
  * Mus_musculus.GRCm39.116.gtf.gz                    (annotation, kept for downstream use)
and builds the bowtie2 index once. Pass --genome-fa / --index-prefix to reuse your own.

The FASTQ files in ``sandbox/`` are plain text even though they are named
``*.fastq.gz``; this script detects the real compression from the file's magic
bytes and feeds fastp a correctly-named symlink, so misnamed inputs just work.

Usage
-----
    pixi install
    pixi run pipeline --fastq-dir ../sandbox --outdir ../results --threads 12

Run ``python atac_to_bam.py --help`` for all options. Nothing here needs the
network except the one-time reference download; use --skip-download once you have
a local genome + index.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

ENSEMBL_FTP = "https://ftp.ensembl.org/pub"
SPECIES = "mus_musculus"
ASSEMBLY = "GRCm39"

# read1 / read2 suffixes we know how to pair, longest first
_MATE_RE = re.compile(
    r"(.+?)(?:_S\d+)?(?:_L\d{3})?[._](?:R?1)(?:_001)?\.(?:fastq|fq)(?:\.gz)?$",
    re.IGNORECASE,
)


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def ts() -> str:
    return _dt.datetime.now().strftime("%H:%M:%S")


def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", flush=True)


def die(msg: str) -> None:
    print(f"[{ts()}] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


DRY_RUN = False


def run(cmd: list[str], *, log_path: Path | None = None) -> None:
    """Run a plain (non-piped) command, aborting on failure."""
    printable = " ".join(str(c) for c in cmd)
    log(f"$ {printable}" + (f"   (2> {log_path})" if log_path else ""))
    if DRY_RUN:
        return
    stderr = open(log_path, "wb") if log_path else None
    try:
        subprocess.run([str(c) for c in cmd], check=True, stderr=stderr)
    except subprocess.CalledProcessError as exc:
        if log_path:
            die(f"command failed (exit {exc.returncode}); see {log_path}")
        die(f"command failed (exit {exc.returncode}): {printable}")
    finally:
        if stderr:
            stderr.close()


def sh(pipeline: str) -> None:
    """Run a bash pipeline with `set -Eeuo pipefail` so mid-pipe failures abort."""
    log("$ " + pipeline)
    if DRY_RUN:
        return
    try:
        subprocess.run(
            ["bash", "-c", "set -Eeuo pipefail\n" + pipeline], check=True
        )
    except subprocess.CalledProcessError as exc:
        die(f"pipeline failed (exit {exc.returncode})")


def need(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        die(
            f"'{tool}' not found on PATH. Run this inside the pixi env: "
            f"`pixi run pipeline ...` (or `pixi shell`)."
        )
    return path


def is_gzip(path: Path) -> bool:
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def human_bytes(n: int) -> str:
    f = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if f < 1024 or unit == "TB":
            return f"{f:.1f} {unit}"
        f /= 1024
    return f"{f:.1f} TB"


def total_ram_gb() -> float:
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1024**3
    except (ValueError, OSError):
        return 0.0


# --------------------------------------------------------------------------- #
# sample discovery
# --------------------------------------------------------------------------- #
def discover_samples(fastq_dir: Path, only: set[str] | None) -> list[dict]:
    if not fastq_dir.is_dir():
        die(f"--fastq-dir does not exist: {fastq_dir}")

    pairs: list[dict] = []
    for r1 in sorted(fastq_dir.iterdir()):
        if not r1.is_file():
            continue
        m = _MATE_RE.match(r1.name)
        if not m:
            continue
        # derive the read2 name by swapping the mate token
        r2_name = re.sub(
            r"([._])(R?)1((?:_001)?\.(?:fastq|fq)(?:\.gz)?)$",
            r"\g<1>\g<2>2\g<3>",
            r1.name,
            flags=re.IGNORECASE,
        )
        r2 = r1.with_name(r2_name)
        if not r2.exists():
            log(f"skip {r1.name}: no matching read-2 file ({r2_name})")
            continue

        sample = m.group(1)
        sample = re.sub(r"_downsample$", "", sample)
        sample = re.sub(r"[^A-Za-z0-9._-]+", "_", sample).strip("_")
        if only and sample not in only:
            continue
        pairs.append({"sample": sample, "r1": r1, "r2": r2})

    if not pairs:
        die(
            f"no paired FASTQ files found in {fastq_dir}. Expected names like "
            f"'<sample>_1.fastq.gz' / '<sample>_2.fastq.gz'."
        )
    return pairs


def staged_fastq(src: Path, work: Path) -> Path:
    """
    Return a path whose extension matches the real compression of *src*.
    The sandbox files are plain text named '*.fastq.gz'; fastp trusts the
    extension, so we hand it a correctly-named symlink instead of the original.
    """
    gz = is_gzip(src)
    stem = src.name
    for ext in (".gz", ".fastq", ".fq"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    target = work / (f"{stem}.fastq.gz" if gz else f"{stem}.fastq")
    if DRY_RUN:
        return target
    if target.exists() or target.is_symlink():
        target.unlink()
    target.symlink_to(src.resolve())
    if gz != src.name.lower().endswith(".gz"):
        log(
            f"note: {src.name} is {'gzip' if gz else 'plain text'} despite its "
            f"name -> staged as {target.name}"
        )
    return target


# --------------------------------------------------------------------------- #
# reference
# --------------------------------------------------------------------------- #
def ensembl_urls(release: int) -> dict[str, str]:
    base = f"{ENSEMBL_FTP}/release-{release}"
    return {
        "genome": f"{base}/fasta/{SPECIES}/dna/"
        f"Mus_musculus.{ASSEMBLY}.dna.primary_assembly.fa.gz",
        "gtf": f"{base}/gtf/{SPECIES}/Mus_musculus.{ASSEMBLY}.{release}.gtf.gz",
        "checksums": f"{base}/fasta/{SPECIES}/dna/CHECKSUMS",
    }


def download(url: str, dest: Path, force: bool) -> None:
    if dest.exists() and not force:
        log(f"have {dest.name} ({human_bytes(dest.stat().st_size)}) -- skip download")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "curl", "-fL", "--retry", "5", "--retry-delay", "5",
            "-C", "-", "-o", str(dest), url,
        ]
    )


def verify_checksum(dest: Path, checksums_url: str) -> None:
    """Best-effort check against Ensembl's BSD `sum` CHECKSUMS file (warn only)."""
    try:
        text = subprocess.run(
            ["curl", "-fsL", checksums_url],
            check=True, capture_output=True, text=True,
        ).stdout
        got = subprocess.run(
            ["sum", str(dest)], check=True, capture_output=True, text=True
        ).stdout.split()
    except (subprocess.CalledProcessError, FileNotFoundError):
        log("checksum: could not verify (non-fatal)")
        return
    want = {
        parts[2]: (parts[0], parts[1])
        for line in text.splitlines()
        if len(parts := line.split()) == 3
    }
    exp = want.get(dest.name)
    if exp and tuple(got[:2]) == exp:
        log(f"checksum OK: {dest.name}")
    elif exp:
        log(f"WARNING: checksum mismatch for {dest.name} (expected {exp}, got {tuple(got[:2])})")
    else:
        log(f"checksum: {dest.name} not listed in CHECKSUMS (skipped)")


def prepare_reference(args) -> tuple[Path, Path | None]:
    refdir: Path = args.refdir
    if not DRY_RUN:
        refdir.mkdir(parents=True, exist_ok=True)

    # ---- genome FASTA (decompressed, for bowtie2-build) --------------------
    if args.genome_fa:
        fa = Path(args.genome_fa)
        if not fa.exists():
            die(f"--genome-fa not found: {fa}")
    else:
        urls = ensembl_urls(args.ensembl_release)
        fa_gz = refdir / Path(urls["genome"]).name
        fa = fa_gz.with_suffix("")  # strip .gz
        if args.skip_download and not (fa.exists() or fa_gz.exists()):
            die("--skip-download set but no genome FASTA present in --refdir")
        if not fa.exists():
            download(urls["genome"], fa_gz, args.force_download)
            if not DRY_RUN:
                verify_checksum(fa_gz, urls["checksums"])
            log(f"decompressing {fa_gz.name} -> {fa.name}")
            sh(f"pigz -d -k -f -p {args.threads} {shq(fa_gz)}")
        # keep the annotation next to the genome for downstream tools
        gtf_gz = refdir / Path(urls["gtf"]).name
        if not args.skip_download:
            download(urls["gtf"], gtf_gz, args.force_download)

    # ---- bowtie2 index ---------------------------------------------------- #
    if args.index_prefix:
        idx = Path(args.index_prefix)
        if not Path(f"{idx}.1.bt2").exists() and not Path(f"{idx}.1.bt2l").exists():
            die(f"--index-prefix has no .bt2 files: {idx}")
        return idx, fa

    idx = refdir / f"bowtie2_{ASSEMBLY}"
    if Path(f"{idx}.1.bt2").exists() or Path(f"{idx}.1.bt2l").exists():
        log(f"bowtie2 index present: {idx}.*.bt2")
        return idx, fa

    if args.skip_download:
        die("--skip-download set but no bowtie2 index present; drop the flag once")

    build_threads = max(1, min(args.threads, 4))  # memory scales with threads
    packed = "--packed " if (args.packed or (total_ram_gb() and total_ram_gb() < 16)) else ""
    if packed:
        log(f"detected {total_ram_gb():.1f} GB RAM -> building index with --packed (slower, lower memory)")
    log("building bowtie2 index -- mouse genome, expect ~30-75 min")
    sh(
        f"bowtie2-build {packed}--threads {build_threads} "
        f"{shq(fa)} {shq(idx)} > {shq(args.refdir / 'bowtie2_build.log')} 2>&1"
    )
    return idx, fa


def shq(p) -> str:
    """minimal shell quoting for paths used in sh() pipelines"""
    s = str(p)
    return "'" + s.replace("'", "'\\''") + "'" if re.search(r"[^\w./-]", s) else s


# --------------------------------------------------------------------------- #
# per-sample processing
# --------------------------------------------------------------------------- #
def process_sample(s: dict, idx: Path, args) -> dict:
    name = s["sample"]
    out: Path = args.outdir
    work: Path = args.workdir / name
    tmp: Path = args.tmpdir / name
    if not DRY_RUN:
        for d in (out, work, tmp):
            d.mkdir(parents=True, exist_ok=True)

    log(f"=== sample {name} ===")

    r1 = staged_fastq(s["r1"], work)
    r2 = staged_fastq(s["r2"], work)

    # 1. trimming ---------------------------------------------------------- #
    t1 = work / f"{name}.trim_1.fastq.gz"
    t2 = work / f"{name}.trim_2.fastq.gz"
    fastp_threads = min(args.threads, 16)
    if args.resume and t1.exists() and t2.exists():
        log("trimmed FASTQ present -- skip fastp (--resume)")
    else:
        run(
            [
                "fastp",
                "-i", r1, "-I", r2,
                "-o", t1, "-O", t2,
                "--detect_adapter_for_pe",
                # Tn5 mosaic-end (Nextera) adapter -- QC of these files showed it
                # in ~38 % of reads, so pass it explicitly as well as auto-detect
                "--adapter_sequence", "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC",
                "--adapter_sequence_r2", "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
                "--length_required", "20",
                "--thread", str(fastp_threads),
                "--json", out / f"{name}.fastp.json",
                "--html", out / f"{name}.fastp.html",
            ],
            log_path=out / f"{name}.fastp.log",
        )

    # 2-5. align | fixmate | sort | markdup ------------------------------- #
    markdup_bam = out / f"{name}.markdup.bam"
    sort_threads = max(1, min(args.threads, 4))
    if args.resume and markdup_bam.exists():
        log("markdup BAM present -- skip alignment (--resume)")
        if not DRY_RUN and not Path(f"{markdup_bam}.bai").exists():
            run(["samtools", "index", "-@", str(sort_threads), markdup_bam])
    else:
        sh(
            f"bowtie2 --very-sensitive --dovetail -X 2000 --mm "
            f"-p {args.threads} -x {shq(idx)} -1 {shq(t1)} -2 {shq(t2)} "
            f"--rg-id {shq(name)} --rg SM:{shq(name)} --rg PL:ILLUMINA --rg LB:{shq(name)} "
            f"2> {shq(out / f'{name}.bowtie2.log')} "
            f"| samtools fixmate -m -u -@ {sort_threads} - - "
            f"| samtools sort -u -@ {sort_threads} -m {shq(args.sort_mem)} "
            f"      -T {shq(tmp / 'sort')} - "
            f"| samtools markdup -@ {sort_threads} -T {shq(tmp / 'mdup')} "
            f"      -f {shq(out / f'{name}.markdup_stats.txt')} - {shq(markdup_bam)}"
        )
        run(["samtools", "index", "-@", str(sort_threads), markdup_bam])
        run(
            ["bash", "-c",
             f"samtools flagstat -@ {sort_threads} {shq(markdup_bam)} "
             f"> {shq(out / f'{name}.markdup.flagstat.txt')}"]
        )

    # 6. filtering ------------------------------------------------------- #
    # keep every reference sequence except the mitochondrion (Ensembl: 'MT')
    keep_bed = work / f"{name}.keep.bed"
    if not DRY_RUN:
        idxstats = subprocess.run(
            ["samtools", "idxstats", str(markdup_bam)],
            check=True, capture_output=True, text=True,
        ).stdout
        with open(keep_bed, "w") as fh:
            for line in idxstats.splitlines():
                chrom, length, *_ = line.split("\t")
                if chrom in ("*", "MT", "chrM"):
                    continue
                if args.main_chroms_only and not re.fullmatch(r"(chr)?([0-9]{1,2}|X|Y)", chrom):
                    continue
                fh.write(f"{chrom}\t0\t{length}\n")

    filtered_bam = out / f"{name}.filtered.bam"
    # -f 2   : properly paired
    # -F 3852: drop unmapped, mate-unmapped, secondary, QCfail, duplicate, supplementary
    sh(
        f"samtools view -@ {sort_threads} -b -f 2 -F 3852 -q {args.mapq} "
        f"-L {shq(keep_bed)} -o {shq(filtered_bam)} {shq(markdup_bam)}"
    )
    run(["samtools", "index", "-@", str(sort_threads), filtered_bam])
    run(["bash", "-c",
         f"samtools flagstat -@ {sort_threads} {shq(filtered_bam)} "
         f"> {shq(out / f'{name}.filtered.flagstat.txt')}"])
    run(["bash", "-c",
         f"samtools idxstats {shq(filtered_bam)} "
         f"> {shq(out / f'{name}.filtered.idxstats.txt')}"])

    if not args.keep_intermediate:
        shutil.rmtree(work, ignore_errors=True)
        shutil.rmtree(tmp, ignore_errors=True)
        if not args.keep_markdup:
            for p in (markdup_bam, Path(f"{markdup_bam}.bai")):
                p.unlink(missing_ok=True)

    log(f"--- {name} done -> {filtered_bam}")
    return {"sample": name, "filtered_bam": filtered_bam, "markdup_bam": markdup_bam}


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def parse_args(argv=None) -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__),
    )
    io = p.add_argument_group("inputs / outputs")
    io.add_argument("--fastq-dir", type=Path, default=here.parent / "sandbox",
                    help="directory with paired *_1/*_2 FASTQ files (default: ../sandbox)")
    io.add_argument("--outdir", type=Path, default=here.parent / "results",
                    help="where BAMs + QC land (default: ../results)")
    io.add_argument("--refdir", type=Path, default=here.parent / "reference",
                    help="genome / index cache (default: ../reference)")
    io.add_argument("--workdir", type=Path, default=None,
                    help="scratch for trimmed FASTQ etc (default: <outdir>/_work)")
    io.add_argument("--tmpdir", type=Path, default=None,
                    help="samtools sort/markdup temp (default: <workdir>/_tmp)")
    io.add_argument("--samples", nargs="*", default=None,
                    help="only process these sample names")

    ref = p.add_argument_group("reference")
    ref.add_argument("--ensembl-release", type=int, default=116)
    ref.add_argument("--genome-fa", default=None,
                     help="use this (decompressed) genome FASTA instead of downloading")
    ref.add_argument("--index-prefix", default=None,
                     help="use this existing bowtie2 index prefix")
    ref.add_argument("--skip-download", action="store_true",
                     help="fail rather than touch the network (needs local genome+index)")
    ref.add_argument("--force-download", action="store_true")
    ref.add_argument("--packed", action="store_true",
                     help="force bowtie2-build --packed (auto-on when RAM < 16 GB)")

    tune = p.add_argument_group("tuning")
    tune.add_argument("--threads", type=int, default=min(multiprocessing.cpu_count(), 12))
    tune.add_argument("--sort-mem", default="512M", help="samtools sort -m (per thread)")
    tune.add_argument("--mapq", type=int, default=30, help="minimum MAPQ in the filtered BAM")
    tune.add_argument("--main-chroms-only", action="store_true",
                      help="keep only chr 1-19, X, Y (drop unplaced scaffolds)")
    tune.add_argument("--keep-intermediate", action="store_true",
                      help="keep trimmed FASTQ and scratch dirs")
    tune.add_argument("--keep-markdup", action="store_true",
                      help="keep the pre-filter markdup BAM (kept by default unless --keep-intermediate off)")
    tune.add_argument("--resume", action="store_true",
                      help="skip a per-sample step if its output already exists")
    tune.add_argument("--dry-run", action="store_true",
                      help="print every command without running it")

    args = p.parse_args(argv)
    args.fastq_dir = args.fastq_dir.resolve()
    args.outdir = args.outdir.resolve()
    args.refdir = args.refdir.resolve()
    args.workdir = (args.workdir or args.outdir / "_work").resolve()
    args.tmpdir = (args.tmpdir or args.workdir / "_tmp").resolve()
    if args.keep_intermediate:
        args.keep_markdup = True
    return args


def main(argv=None) -> None:
    global DRY_RUN
    args = parse_args(argv)
    DRY_RUN = args.dry_run

    for t in ("bowtie2", "bowtie2-build", "samtools", "fastp", "pigz", "curl", "sum"):
        need(t)

    ram = total_ram_gb()
    log(f"threads={args.threads}  RAM={ram:.1f} GB  fastq-dir={args.fastq_dir}")
    if ram and ram < 8:
        log("WARNING: <8 GB RAM -- bowtie2 index build may fail; see README compute notes")

    samples = discover_samples(args.fastq_dir, set(args.samples) if args.samples else None)
    log(f"{len(samples)} sample(s): " + ", ".join(s['sample'] for s in samples))

    if not DRY_RUN:
        args.outdir.mkdir(parents=True, exist_ok=True)
    idx, _fa = prepare_reference(args)

    results = [process_sample(s, idx, args) for s in samples]

    # optional MultiQC roll-up
    if shutil.which("multiqc") and not DRY_RUN:
        log("multiqc: summarising fastp + samtools stats")
        subprocess.run(
            ["multiqc", "-f", "-o", str(args.outdir), "-n", "multiqc_report",
             str(args.outdir)],
            check=False,
        )

    print("\n" + "=" * 64)
    log("DONE")
    for r in results:
        print(f"  {r['sample']:<28}  {r['filtered_bam']}")
    print("=" * 64)


if __name__ == "__main__":
    main()
