#!/bin/bash
#SBATCH --job-name=cutrun
#SBATCH --out="slurm-cutrun-%A_%a.out"
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --array=0-5

# Stop on the first error so a failed step doesn't silently cascade.
set -euo pipefail

# ===========================================================================
# CUT&RUN pipeline (one SLURM array task per sample).
#
# You normally do NOT run this directly. Submit it with the wrapper, which
# discovers the samples and sizes the array for you:
#
#     bash run_cutrun.sh "/home/jz954/<your fastq folder>"
#
# To run it directly, export DATA_ROOT (and optionally OUT_ROOT) first and make
# sure --array matches the sample count, e.g.:
#
#     export DATA_ROOT="/home/jz954/C&R GLI2 EX3"
#     sbatch --array=0-5 cutrun.sh
# ===========================================================================

# ---------------------------------------------------------------------------
# Fixed environment paths -- normally never change these.
# IMPORTANT: use the /home/jz954 path, NOT /vast/palmer/home.mccleary/jz954.
# They resolve to the same files on disk, but Apptainer only bind-mounts
# /home/jz954 into the container (the raw /vast/palmer/home.mccleary path is
# NOT bound), so any path used *inside* an `apptainer exec` call -- the image
# and everything in src_dir -- must live under /home/jz954.
# ---------------------------------------------------------------------------
image=/home/jz954/jz_projects/CUTRUN_JZ/yq        # container built by build.sh
src_dir=/home/jz954/jz_projects/CUTRUN_JZ/src      # Trimmomatic jar, picard.jar, index/mm10

# ---------------------------------------------------------------------------
# Per-run inputs (set by run_cutrun.sh, or exported before `sbatch cutrun.sh`).
# data_root holds one sub-folder per sample, each with Unaligned/*_R1_001.fastq.gz
# (the layout produced by `ycgaFastq -o <data_root> <URL>`).
# ---------------------------------------------------------------------------
data_root=${DATA_ROOT:?Set DATA_ROOT to the fastq folder, or submit via: bash run_cutrun.sh "<folder>"}
# Outputs default to palmer_scratch/cutrun_out_<folder name> (spaces -> underscores).
out_root=${OUT_ROOT:-$HOME/palmer_scratch/cutrun_out_$(basename "$data_root" | tr ' ' '_')}

# The input FASTQs are symlinks into YCGA's /gpfs/ycga sequencer storage (via the
# /ycga-gpfs alias). Apptainer does not bind that filesystem by default, so the
# symlinks resolve to "No such file" inside the container. Expose both paths to
# every `apptainer exec` below so the links resolve.
export APPTAINER_BIND="/gpfs/ycga,/ycga-gpfs"

# Discover samples: every immediate sub-folder of data_root that has an R1 FASTQ
# under Unaligned/. Sorted + de-duped so the array index -> sample mapping is
# stable. This MUST match the discovery in run_cutrun.sh so the indices line up.
shopt -s nullglob
samples=()
for r1 in "$data_root"/*/Unaligned/*_R1_001.fastq.gz; do
	s=${r1%/Unaligned/*}
	samples+=("${s##*/}")
done
mapfile -t samples < <(printf '%s\n' "${samples[@]}" | sort -u)

(( ${#samples[@]} > 0 )) || { echo "No samples found under $data_root/*/Unaligned/*_R1_001.fastq.gz" >&2; exit 1; }
if (( SLURM_ARRAY_TASK_ID >= ${#samples[@]} )); then
	echo "Array task $SLURM_ARRAY_TASK_ID >= sample count ${#samples[@]} -- nothing to do for this task."
	exit 0
fi
sample=${samples[$SLURM_ARRAY_TASK_ID]}

threads=${SLURM_CPUS_PER_TASK:-8}

# Input FASTQs (R1/R2) live under <data_root>/<sample>/Unaligned/
in_dir=$data_root/$sample/Unaligned
r1=$(ls "$in_dir"/*_R1_001.fastq.gz)
r2=$(ls "$in_dir"/*_R2_001.fastq.gz)

# Each sample gets its own clean output folder.
dir=$out_root/$sample
mkdir -p "$dir"

echo "=== Sample $sample (array task $SLURM_ARRAY_TASK_ID of ${#samples[@]}) ==="
echo "Data root: $data_root"
echo "R1: $r1"
echo "R2: $r2"
echo "Output dir: $dir"

# 1. Quality report on the raw reads
apptainer exec "$image" fastqc -o "$dir" "$r1" "$r2"

# 2. Trim leftover adapters and drop reads shorter than 20 bp
apptainer exec "$image" java -jar "$src_dir/Trimmomatic-0.36/trimmomatic-0.36.jar" PE -phred33 \
	"$r1" "$r2" \
	"$dir/paired_r1.fq.gz" "$dir/unpaired_r1.fq.gz" \
	"$dir/paired_r2.fq.gz" "$dir/unpaired_r2.fq.gz" \
	ILLUMINACLIP:"$src_dir/Trimmomatic-0.36/adapters/TruSeq3-PE.fa":2:15:4:1:true MINLEN:20

# 3. Align to mm10, convert SAM->BAM (pipe runs inside one container call)
apptainer exec "$image" bash -c \
	"bowtie2 --dovetail --threads $threads -x '$src_dir/index/mm10' \
	 -1 '$dir/paired_r1.fq.gz' -2 '$dir/paired_r2.fq.gz' \
	 | samtools view -bS - > '$dir/paired.bam'"

# 4. Sort by genome coordinate
apptainer exec "$image" java -jar "$src_dir/picard.jar" SortSam \
	I="$dir/paired.bam" O="$dir/sorted_paired.bam" \
	SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

# 5. Flag PCR duplicates
apptainer exec "$image" java -jar "$src_dir/picard.jar" MarkDuplicates \
	I="$dir/sorted_paired.bam" O="$dir/marked_dup_paired.bam" \
	M="$dir/marked_dup_paired.txt" VALIDATION_STRINGENCY=LENIENT

# 6. Remove duplicates (-F 1024) and keep only properly paired reads (-f 2)
apptainer exec "$image" samtools view -F 1024 -f 2 -b \
	"$dir/marked_dup_paired.bam" -o "$dir/removed_dup_paired.bam"

# 7. Index the final BAM
apptainer exec "$image" samtools index "$dir/removed_dup_paired.bam" "$dir/removed_dup_paired.bai"

# 8. Build the coverage track
apptainer exec "$image" bamCoverage -b "$dir/removed_dup_paired.bam" -o "$dir/deeptools_coverage.bw"

echo "=== Done: $dir/deeptools_coverage.bw ==="
