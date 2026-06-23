# CUTRUN

CUT&RUN pipeline for mouse (mm10) paired-end FASTQ data. Each sample goes from
raw reads to a BigWig (`.bw`) coverage track. All tools run inside an Apptainer
(Singularity) container; the pipeline runs as a SLURM array job (one task per
sample) on the HPC cluster.

## One-time setup

```bash
bash download.sh        # picard.jar, Trimmomatic-0.36, samtools-1.17 -> ./src
bash download_index.sh  # mm10 bowtie2 index -> ./src/index
bash build.sh           # build the 'yq' Apptainer container from yq.def
```

The container's `%files` line in `yq.def` and the `image`/`src_dir` paths in
`cutrun.sh` point at this repo under `/home/jz954/jz_projects/CUTRUN_JZ`.

> **Note (McCleary):** use the `/home/jz954/...` path, not the equivalent
> `/vast/palmer/home.mccleary/jz954/...` path. They are the same files on disk,
> but Apptainer only bind-mounts `/home/jz954` into the container, so any path
> used *inside* a container call (the image, the `src` jars and index) must be
> under `/home/jz954`. This is already baked into the scripts.

## Running the pipeline (per dataset, no code edits)

1. Pull the FASTQ symlinks into a per-dataset folder with YCGA's `ycgaFastq`.
   This creates one `<sample>/Unaligned/*_R{1,2}_001.fastq.gz` sub-folder per
   sample:

   ```bash
   ycgaFastq -f -o "/home/jz954/C&R GLI2 EX3" \
     fcb.ycga.yale.edu:3010/<token>/sample_dir_XXXXXXXX
   ```

2. Submit the pipeline for that folder. The wrapper discovers the samples,
   sizes the SLURM array automatically, and submits one task per sample:

   ```bash
   bash run_cutrun.sh "/home/jz954/C&R GLI2 EX3"
   ```

   Outputs go to `~/palmer_scratch/cutrun_out_<folder name>/<sample>/`, with the
   final track at `deeptools_coverage.bw`. To choose the output folder
   explicitly, pass it as a second argument:

   ```bash
   bash run_cutrun.sh "/home/jz954/C&R GLI2 EX3" "/home/jz954/palmer_scratch/my_out"
   ```

`cutrun.sh` is the array job itself; `run_cutrun.sh` is the wrapper you call.
You can still run `cutrun.sh` directly by exporting `DATA_ROOT` (and optionally
`OUT_ROOT`) and passing a matching `--array`, but `run_cutrun.sh` is the easy path.

## Pipeline stages

FastQC → Trimmomatic trim → bowtie2 align to mm10 → samtools SAM→BAM →
Picard SortSam → Picard MarkDuplicates → samtools filter (`-F 1024 -f 2`: drop
duplicates, keep proper pairs) → samtools index → deeptools `bamCoverage` → `.bw`.
