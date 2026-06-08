process PBGATK_GERMLINE {

    tag "$sample_id"

    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram'
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'

    // =========================
    // SLURM GPU QUEUE
    // =========================
    queue { params.slurm_gpu_queue }

    // =========================
    // FIXED GPU + SLURM OPTIONS
    // =========================
    clusterOptions {
        def req = (gpu_profile ==~ /[124]gpu/) ?
            (gpu_profile.replace('gpu','') as int) :
            (params.slurm_gpu_request as int)

        // IMPORTANT:
        // Removed localscratch completely (this was breaking sbatch)
        "--gres=gpu:${params.slurm_gpu_type}:${req}"
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
    mkdir -p tmp

    pbrun germline \
        --ref ${ref} \
        --in-fq ${read1.join(' ')} ${read2.join(' ')} \
        --out-bam ${sample_id}.cram \
        --out-variants ${sample_id}.vcf \
        --tmp-dir \${PWD}/tmp \
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