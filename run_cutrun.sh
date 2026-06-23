#!/bin/bash
# Submit the CUT&RUN pipeline for a folder of samples, with no code edits.
#
# Usage:
#     bash run_cutrun.sh "<data_dir>" ["<out_dir>"]
#
# <data_dir> is a folder holding one sub-folder per sample, each containing
# Unaligned/*_R1_001.fastq.gz and *_R2_001.fastq.gz -- exactly the layout that
#     ycgaFastq -o "<data_dir>" <URL>
# produces. The wrapper discovers the samples, sizes the SLURM array to match,
# and submits cutrun.sh (one array task per sample).
#
# Example:
#     ycgaFastq -f -o "/home/jz954/C&R GLI2 EX3" fcb.ycga.yale.edu:3010/.../sample_dir_XXXX
#     bash run_cutrun.sh "/home/jz954/C&R GLI2 EX3"
#
# Outputs default to: ~/palmer_scratch/cutrun_out_<data_dir name> (spaces -> _).

set -euo pipefail

data_root=${1:?Usage: bash run_cutrun.sh "<data_dir>" ["<out_dir>"]}
out_root=${2:-$HOME/palmer_scratch/cutrun_out_$(basename "$data_root" | tr ' ' '_')}
here=$(cd "$(dirname "$0")" && pwd)

[ -d "$data_root" ] || { echo "Data dir not found: $data_root" >&2; exit 1; }

# Discover samples the SAME way cutrun.sh does, so the array indices line up.
shopt -s nullglob
samples=()
for r1 in "$data_root"/*/Unaligned/*_R1_001.fastq.gz; do
	s=${r1%/Unaligned/*}
	samples+=("${s##*/}")
done
mapfile -t samples < <(printf '%s\n' "${samples[@]}" | sort -u)

n=${#samples[@]}
(( n > 0 )) || { echo "No samples found under $data_root/*/Unaligned/*_R1_001.fastq.gz" >&2; exit 1; }

echo "Found $n sample(s) under: $data_root"
printf '  %s\n' "${samples[@]}"
echo "Outputs -> $out_root"
mkdir -p "$out_root"

# Export the per-run paths into the job's environment (SLURM's default
# --export=ALL propagates them, and this safely carries values with spaces).
export DATA_ROOT="$data_root"
export OUT_ROOT="$out_root"

# Command-line --array overrides the #SBATCH --array directive in cutrun.sh.
sbatch --array=0-$((n - 1)) "$here/cutrun.sh"
