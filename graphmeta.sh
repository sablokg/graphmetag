#!/usr/bin/env bash
#
# graph_metagenomics_pipeline.sh Gaurav Sablok gsablok@proton.me
#
# Graph-based metagenome assembly + binning pipeline.
#
# Steps:
#   1. QC / adapter trimming        (fastp)
#   2. Metagenomic assembly with an explicit assembly graph (MEGAHIT -> .fastg, or metaSPAdes -> .gfa)
#   3. Read mapping back to contigs for coverage/abundance    (bwa + samtools)
#   4. Initial contig binning                                  (MetaBAT2)
#   5. Graph-aware bin refinement using the assembly graph      (GraphBin)
#
# Requirements (install via conda/mamba is easiest):
#   fastp, megahit (or spades.py), bwa, samtools, metabat2, graphbin
#
# Usage:
#   ./graph_metagenomics_pipeline.sh -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results -t 16
#
set -euo pipefail

# ---------- defaults ----------
THREADS=8
OUTDIR="results"
ASSEMBLER="megahit"   # megahit | metaspades

usage() {
  echo "Usage: $0 -1 <R1.fastq.gz> -2 <R2.fastq.gz> -o <outdir> [-t threads] [-a megahit|metaspades]"
  exit 1
}

while getopts "1:2:o:t:a:h" opt; do
  case "$opt" in
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    a) ASSEMBLER="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${R1:-}" || -z "${R2:-}" ]] && usage

mkdir -p "$OUTDIR"/{qc,assembly,mapping,binning/metabat2,binning/graphbin}

# ---------- 1. QC / trimming ----------
echo "[1/5] Running fastp QC..."
fastp \
  -i "$R1" -I "$R2" \
  -o "$OUTDIR/qc/R1.trim.fastq.gz" -O "$OUTDIR/qc/R2.trim.fastq.gz" \
  --thread "$THREADS" \
  --json "$OUTDIR/qc/fastp.json" --html "$OUTDIR/qc/fastp.html"

CLEAN_R1="$OUTDIR/qc/R1.trim.fastq.gz"
CLEAN_R2="$OUTDIR/qc/R2.trim.fastq.gz"

# ---------- 2. Assembly with graph output ----------
echo "[2/5] Assembling with $ASSEMBLER..."
if [[ "$ASSEMBLER" == "megahit" ]]; then
  rm -rf "$OUTDIR/assembly/megahit_out"
  megahit \
    -1 "$CLEAN_R1" -2 "$CLEAN_R2" \
    -o "$OUTDIR/assembly/megahit_out" \
    -t "$THREADS"

  CONTIGS="$OUTDIR/assembly/megahit_out/final.contigs.fa"

  # Convert MEGAHIT's intermediate contig graph to FASTG for downstream graph tools
  LAST_K=$(ls "$OUTDIR"/assembly/megahit_out/intermediate_contigs/k*.contigs.fa \
            | sed -E 's/.*k([0-9]+)\.contigs\.fa/\1/' | sort -n | tail -1)
  megahit_toolkit contig2fastg "$LAST_K" \
    "$OUTDIR/assembly/megahit_out/intermediate_contigs/k${LAST_K}.contigs.fa" \
    > "$OUTDIR/assembly/assembly_graph.fastg"
  GRAPH="$OUTDIR/assembly/assembly_graph.fastg"
  GRAPH_FORMAT="fastg"

elif [[ "$ASSEMBLER" == "metaspades" ]]; then
  spades.py --meta \
    -1 "$CLEAN_R1" -2 "$CLEAN_R2" \
    -o "$OUTDIR/assembly/metaspades_out" \
    -t "$THREADS"

  CONTIGS="$OUTDIR/assembly/metaspades_out/contigs.fasta"
  GRAPH="$OUTDIR/assembly/metaspades_out/assembly_graph_with_scaffolds.gfa"
  GRAPH_FORMAT="gfa"
else
  echo "Unknown assembler: $ASSEMBLER (use megahit or metaspades)"; exit 1
fi

# ---------- 3. Map reads back for coverage ----------
echo "[3/5] Mapping reads back to contigs for coverage..."
bwa index "$CONTIGS"
bwa mem -t "$THREADS" "$CONTIGS" "$CLEAN_R1" "$CLEAN_R2" \
  | samtools sort -@ "$THREADS" -o "$OUTDIR/mapping/reads.sorted.bam" -
samtools index "$OUTDIR/mapping/reads.sorted.bam"

# ---------- 4. Initial binning ----------
echo "[4/5] Binning contigs with MetaBAT2..."
jgi_summarize_bam_contig_depths \
  --outputDepth "$OUTDIR/binning/metabat2/depth.txt" \
  "$OUTDIR/mapping/reads.sorted.bam"

metabat2 \
  -i "$CONTIGS" \
  -a "$OUTDIR/binning/metabat2/depth.txt" \
  -o "$OUTDIR/binning/metabat2/bin" \
  -t "$THREADS"

# Build the contig->bin mapping GraphBin expects
python3 - "$OUTDIR/binning/metabat2" "$OUTDIR/binning/metabat2/initial_binning.csv" <<'PYEOF'
import sys, glob, os
bindir, outpath = sys.argv[1], sys.argv[2]
with open(outpath, "w") as out:
    for bin_id, fasta in enumerate(sorted(glob.glob(os.path.join(bindir, "bin.*.fa")))):
        with open(fasta) as f:
            for line in f:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    out.write(f"{contig},{bin_id}\n")
PYEOF

# ---------- 5. Graph-aware bin refinement ----------
echo "[5/5] Refining bins using the assembly graph (GraphBin)..."
graphbin \
  --assembler "$ASSEMBLER" \
  --graph "$GRAPH" \
  --contigs "$CONTIGS" \
  --binned "$OUTDIR/binning/metabat2/initial_binning.csv" \
  --output "$OUTDIR/binning/graphbin"

echo "Done. Refined bins are in: $OUTDIR/binning/graphbin"
