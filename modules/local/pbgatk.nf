process PBGATK_GERMLINE {

    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram'
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'

    queue { params.slurm_gpu_queue }

    clusterOptions {
        def req = (gpu_profile ==~ /[124]gpu/) ? (gpu_profile.replace('gpu','') as int) : (params.slurm_gpu_request as int)
        "--gres=ssd,gpu:${params.slurm_gpu_type}:${req} --localscratch=ssd:${params.slurm_gpu_localscratch_gb}"
    }

    input:
    tuple val(sample_id), path(read1), path(read2), val(gpu_profile)
    tuple path(ref), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_fai)

    output:
    tuple val(sample_id), path("${sample_id}.cram"), emit: cram
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:

    def gpus = (gpu_profile ==~ /[124]gpu/) ? gpu_profile.replace('gpu','') as int : (params.call_default_gpus as int)
    def runPartition = gpus > 1 ? '--run-partition' : ''

    def r1List = (read1 instanceof List) ? read1 : [read1]
    def r2List = (read2 instanceof List) ? read2 : [read2]

    if( r1List.size() != r2List.size() )
        error "Sample ${sample_id}: read1 count (${r1List.size()}) != read2 count (${r2List.size()})"

    if( r1List.isEmpty() )
        error "Sample ${sample_id}: no FASTQ pairs resolved"

    def sortedR1 = r1List.sort { it.name }
    def sortedR2 = r2List.sort { it.name }

    def inFqArgs = sortedR1.indices.collect { i ->
        "--in-fq ${sortedR1[i]} ${sortedR2[i]}"
    }.join(' ')

    """
    echo "========== ENV =========="
    echo "HOSTNAME=\$(hostname)"
    echo "SLURM_JOB_ID=\$SLURM_JOB_ID"
    echo "SLURM_LOCAL_SCRATCH=\$SLURM_LOCAL_SCRATCH"
    echo "SLURM_TMPDIR=\$SLURM_TMPDIR"
    echo "TMPDIR=\$TMPDIR"
    nvidia-smi || true
    echo "========================="

    if [[ -z "\${SLURM_LOCAL_SCRATCH:-}" ]]; then
        echo "ERROR: SLURM_LOCAL_SCRATCH is not set"
        exit 1
    fi

    SCRATCH_DIR="\${SLURM_LOCAL_SCRATCH}/pbgatk_${sample_id}"
    mkdir -p "\$SCRATCH_DIR"

    echo "Using scratch directory: \$SCRATCH_DIR"

    pbrun germline \
        --ref ${ref} \
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

    ls -lh \$SCRATCH_DIR

    cp \$SCRATCH_DIR/${sample_id}.cram .
    cp \$SCRATCH_DIR/${sample_id}.vcf .
    """
}