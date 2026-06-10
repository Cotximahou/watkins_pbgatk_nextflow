process PBGATK_GERMLINE {

    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram'
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'

    label 'gpu_call'

    queue params.slurm_gpu_queue

    clusterOptions """
        --gres=ssd,gpu:${params.slurm_gpu_type}:${params.slurm_gpu_request}
        --localscratch=ssd:${params.slurm_gpu_localscratch_gb}
    """

    input:
    tuple val(sample_id), path(read1), path(read2), val(gpu_profile)
    tuple path(ref), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_fai)

    output:
    tuple val(sample_id), path("${sample_id}.cram"), emit: cram
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:

    def gpus = (gpu_profile ==~ /[124]gpu/) ? gpu_profile.replace('gpu','') as int : params.call_default_gpus
    def runPartition = gpus > 1 ? '--run-partition' : ''

    def r1 = read1 instanceof List ? read1 : [read1]
    def r2 = read2 instanceof List ? read2 : [read2]

    def inFqArgs = r1.indices.collect { i ->
        "--in-fq ${r1[i]} ${r2[i]}"
    }.join(' ')

    """
    set -euo pipefail

    echo "HOST=$(hostname)"
    echo "JOB=${SLURM_JOB_ID:-NA}"

    SCRATCH_BASE="\${SLURM_LOCAL_SCRATCH:-\${SLURM_TMPDIR:-/tmp}}"
    SCRATCH_DIR="\$SCRATCH_BASE/pbgatk_${sample_id}"
    mkdir -p "\$SCRATCH_DIR"

    echo "SCRATCH=\$SCRATCH_DIR"

    # IMPORTANT FIX: use symlinks, NOT copies (avoids IO + memory spikes)
    ln -s ${ref}     \$SCRATCH_DIR/${ref.name}
    ln -s ${ref_amb} \$SCRATCH_DIR/${ref.name}.amb
    ln -s ${ref_ann} \$SCRATCH_DIR/${ref.name}.ann
    ln -s ${ref_bwt} \$SCRATCH_DIR/${ref.name}.bwt
    ln -s ${ref_pac} \$SCRATCH_DIR/${ref.name}.pac
    ln -s ${ref_sa}  \$SCRATCH_DIR/${ref.name}.sa
    ln -s ${ref_fai} \$SCRATCH_DIR/${ref.name}.fai

    pbrun germline \
        --ref \$SCRATCH_DIR/${ref.name} \
        ${inFqArgs} \
        --out-bam \$SCRATCH_DIR/${sample_id}.cram \
        --tmp-dir \$SCRATCH_DIR \
        --out-variants \$SCRATCH_DIR/${sample_id}.vcf \
        ${runPartition} \
        --num-gpus ${gpus} \
        --num-htvc-threads ${params.call_htvc_threads} \
        --read-group-sm ${sample_id} \
        --read-group-pl ILLUMINA \
        --read-group-id-prefix ${sample_id} \
        --keep-tmp \
        --x3

    cp \$SCRATCH_DIR/${sample_id}.cram .
    cp \$SCRATCH_DIR/${sample_id}.vcf .
    """
}