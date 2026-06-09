process PBGATK_GERMLINE {

    tag "$sample_id"

    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram'
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'

    queue { params.slurm_gpu_queue }

    clusterOptions {
        def req = (gpu_profile ==~ /[124]gpu/) ?
            (gpu_profile.replace('gpu','') as int) :
            (params.slurm_gpu_request as int)

        "--gres=ssd,gpu:${params.slurm_gpu_type}:${req} --localscratch=ssd:${params.slurm_gpu_localscratch_gb}"
    }

    input:
    tuple val(sample_id), val(read1), val(read2), val(gpu_profile)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.cram"), emit: cram
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    def gpus = (gpu_profile ==~ /[124]gpu/) ?
        (gpu_profile.replace('gpu','') as int) :
        (params.call_default_gpus as int)

    def runPartition = gpus > 1 ? '--run-partition' : ''

    """
    SCRATCH_DIR=\${TMPDIR:-\$PWD/tmp}
    mkdir -p \$SCRATCH_DIR

    pbrun germline \
        --ref ${ref} \
        --in-fq ${read1.join(' ')} ${read2.join(' ')} \
        --out-bam \$PWD/${sample_id}.cram \
        --out-variants \$PWD/${sample_id}.vcf \
        --tmp-dir \$SCRATCH_DIR \
        ${runPartition} \
        --num-gpus ${gpus} \
        --num-htvc-threads ${params.call_htvc_threads} \
        --read-group-sm ${sample_id} \
        --read-group-pl ILLUMINA \
        --read-group-id-prefix ${sample_id} \
        --keep-tmp \
        --x3
    """
}