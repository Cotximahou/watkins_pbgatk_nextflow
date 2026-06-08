# Watkins pbgatk Nextflow Workflow

This folder is intended to host a portable Nextflow implementation of the Watkins pipeline currently run with Slurm scripts.

The original workflow stages are:
1. Per-sample variant calling with Parabricks (`pbrun germline`) -> `cram/*.cram` and `vcf/*.vcf`
2. VCF compression and indexing -> `vcfgz/*.vcf.gz` and `.csi`
3. Per-contig extraction and staged merge with `bcftools merge` -> `chr/watkins_<contig>.vcf.gz`
4. Optional CRAM QC with `samtools flagstat`
5. Preflight GPU profile report from samplesheet (`gpu_profile_counts.tsv`)

## Goal

Make the workflow runnable by anyone with:
- Nextflow
- A container engine (Apptainer/Singularity or Docker)
- Slurm access (for cluster execution)

## Required Software

- Nextflow: 24.x or newer
- Java: 17+
- One container engine:
  - Apptainer/Singularity (recommended on HPC)
  - Docker/Podman (local testing)

## Required Containers

Use these as the default tool containers.

| Tool | Purpose | Required local file |
|---|---|---|
| Parabricks 4.3.2-1 | Germline calling | `parabricks-4.3.2-1.sif` |
| bcftools 1.21 | VCF view/merge/index | `bcftools-1.21.sif` |
| samtools 1.18 | flagstat/QC | `samtools-1.18.sif` |

Notes:
- Place all `.sif` files in one local directory and pass it with `--local_container_dir`, or keep the default `./containers`.
- Offline mode uses `pullMode = 'never'`, so missing `.sif` files will fail fast.

## Cluster Container Definitions

Container behavior is controlled by profile:
- `singularity_offline`: use local pre-staged `.sif` files only.

Optional advanced override:
- If your site has custom container locations, copy [nextflow/cluster-containers.example.yml](nextflow/cluster-containers.example.yml) to `cluster-containers.yml` and provide your own `container_*` paths.
- Use that file with `-params-file cluster-containers.yml`.

Slurm tuning override:
- Copy [nextflow/cluster-slurm.example.yml](nextflow/cluster-slurm.example.yml) to `cluster-slurm.yml`.
- Tune queues, GPU type/count, scratch, CPUs, memory, and walltime.
- Use it with `-params-file cluster-slurm.yml`.
- GPU count is automatic per sample from `gpu_profile` in the samplesheet (`1gpu`, `2gpu`, `4gpu`).
- `slurm_gpu_request` is only a fallback when `gpu_profile` is missing.
- Invalid `gpu_profile` values fail early during input validation.

## Offline Runtime Flags

Set these variables when running in strict offline mode:

```bash
export NXF_SINGULARITY_CACHEDIR=/hpc-home/luk24miz/nextflow/watkins_pbgatk_nextflow/.singularity
export NXF_OFFLINE=true
```

## Expected Inputs

Minimum required inputs for the workflow:
- Reference FASTA: for example `../ref/19082024_paragon_v3.fa`
- Sample sheet CSV describing samples and read files

Recommended sample sheet columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample name |
| `read1` | Absolute path to read 1 FASTQ |
| `read2` | Absolute path to read 2 FASTQ |
| `source` | Optional, e.g. `ecs`, `mybookduo`, `missing` |
| `gpu_profile` | Optional, e.g. `1gpu`, `2gpu`, `4gpu` |

Validation rule:
- `gpu_profile` must be one of `1gpu`, `2gpu`, `4gpu`.
- If `gpu_profile` is omitted, the workflow defaults to `1gpu`.

Example `samplesheet.csv`:

```csv
sample_id,read1,read2,source,gpu_profile
SAMPLE_A,/path/to/SAMPLE_A.R1.fastq.gz,/path/to/SAMPLE_A.R2.fastq.gz,mybookduo,1gpu
SAMPLE_B,/path/to/SAMPLE_B_1.fq.gz,/path/to/SAMPLE_B_2.fq.gz,mybookduo,1gpu
SAMPLE_C,/path/to/SAMPLE_C.f1.fastq.gz,/path/to/SAMPLE_C.r2.fastq.gz,missing,4gpu
```

## Recommended Nextflow Profiles

Define at least these profiles in `nextflow.config`:
- `slurm`: Slurm executor settings
- `singularity`: container engine settings
- `singularity_offline`: explicit offline mode (local `.sif` files)
- `docker`: optional local execution profile

Typical settings to include:
- Queue/partition selection for GPU and CPU jobs
- GPU resource requests for Parabricks processes
- High memory settings for merge steps
- Local scratch settings for heavy merge tasks

## Run Commands

Use this profile for cluster execution:

| Environment | Profile | Container source |
|---|---|---|
| Cluster/compute nodes without internet | `slurm,singularity_offline` | Local `.sif` files only |

### 1. Offline Cluster Mode (no internet required)

Copy/paste:

```bash
NXF_SINGULARITY_CACHEDIR=/hpc-home/luk24miz/nextflow/watkins_pbgatk_nextflow/.singularity \
NXF_OFFLINE=true \
nextflow run main.nf \
  -profile slurm,singularity_offline \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --ref /absolute/path/to/19082024_paragon_v3.fa \
  --local_container_dir /absolute/path/to/containers \
  --outdir /absolute/path/to/results
```

Preflight behavior:
- The workflow writes a GPU profile count report to `results/preflight/gpu_profile_counts.tsv`.
- If `gpu_profile` contains a value other than `1gpu`, `2gpu`, or `4gpu`, the run fails before heavy jobs start.

Resume:

```bash
NXF_SINGULARITY_CACHEDIR=/hpc-home/luk24miz/nextflow/watkins_pbgatk_nextflow/.singularity \
NXF_OFFLINE=true \
nextflow run main.nf -profile slurm,singularity_offline --local_container_dir /absolute/path/to/containers -resume
```

### 2. Optional: custom container mapping via params file

Use this only if you need non-default names/locations for containers.

```bash
cp cluster-containers.example.yml cluster-containers.yml
# edit cluster-containers.yml to match your site
nextflow run main.nf \
  -profile slurm,singularity_offline \
  -params-file cluster-containers.yml \
  -params-file cluster-slurm.yml \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --ref /absolute/path/to/19082024_paragon_v3.fa \
  --outdir /absolute/path/to/results
```

### 3. Optional: Slurm tuning via params file

```bash
cp cluster-slurm.example.yml cluster-slurm.yml
# edit cluster-slurm.yml for your cluster
nextflow run main.nf \
  -profile slurm,singularity_offline \
  -params-file cluster-slurm.yml \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --ref /absolute/path/to/19082024_paragon_v3.fa \
  --outdir /absolute/path/to/results
```

### 4. Optional: run subset contigs

```bash
nextflow run main.nf \
  -profile slurm,singularity_offline \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --ref /absolute/path/to/19082024_paragon_v3.fa \
  --contig_subset 1A,1B,1D \
  --outdir /absolute/path/to/results
```

## Outputs

Expected output structure:

- `results/cram/<sample>.cram`
- `results/vcf/<sample>.vcf`
- `results/vcfgz/<sample>.vcf.gz`
- `results/vcfgz/<sample>.vcf.gz.csi`
- `results/chr/watkins_<contig>.vcf.gz`
- `results/chr/watkins_<contig>.vcf.gz.csi`
- `results/flagstat/<sample>.txt` (if QC enabled)
- `results/preflight/gpu_profile_counts.tsv`

## Reproducibility Checklist

- Pin container images by explicit tag (or digest)
- Commit `main.nf`, `nextflow.config`, and module scripts
- Keep a frozen `samplesheet.csv` per run
- Keep `nextflow.log`, `timeline.html`, `trace.txt`, and `report.html`
- Always use `-resume` for retries

## Troubleshooting

- Container pull/auth failures:
  - Pre-pull SIF files and reference local `.sif` paths in config.
- Merge jobs OOM:
  - Increase merge memory and reduce per-merge chunk size.
- GPU queue pressure:
  - Use sample-level routing (`gpu_profile`) to send heavy samples to larger GPU profiles.
- Missing output files:
  - Re-run with `-resume`; Nextflow will execute only missing or failed tasks.

## Status

This README documents the required runtime and containers, and the workflow code is now present in this folder:
- `main.nf`
- `nextflow.config`
- `modules/local/*.nf`

Quick start:

```bash
NXF_SINGULARITY_CACHEDIR=/hpc-home/luk24miz/nextflow/watkins_pbgatk_nextflow/.singularity \
NXF_OFFLINE=true \
nextflow run main.nf \
  -profile slurm,singularity_offline \
  --samplesheet /absolute/path/to/samplesheet.csv \
  --ref /absolute/path/to/19082024_paragon_v3.fa \
  --local_container_dir /absolute/path/to/containers \
  --outdir /absolute/path/to/results
```
