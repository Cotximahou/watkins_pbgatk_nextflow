process PBGATK_GERMLINE {
    publishDir "${params.outdir}/cram", mode: 'copy', pattern: '*.cram'
    publishDir "${params.outdir}/vcf", mode: 'copy', pattern: '*.vcf'
  queue { params.slurm_gpu_queue }
  clusterOptions {
    def req = (gpu_profile ==~ /[124]gpu/) ? (gpu_profile.replace('gpu', '') as int) : (params.slurm_gpu_request as int)
    "--gres=ssd,gpu:${params.slurm_gpu_type}:${req} --localscratch=ssd:${params.slurm_gpu_localscratch_gb}"
  }

    input:
    tuple val(sample_id), path(read1), path(read2), val(gpu_profile)
    tuple path(ref), path(ref_amb), path(ref_ann), path(ref_bwt), path(ref_pac), path(ref_sa), path(ref_fai)

    output:
    tuple val(sample_id), path("${sample_id}.cram"), emit: cram
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    def gpus = (gpu_profile ==~ /[124]gpu/) ? gpu_profile.replace('gpu', '') as int : (params.call_default_gpus as int)
    def runPartition = gpus > 1 ? '--run-partition' : ''
    def r1List = (read1 instanceof List) ? read1 : [read1]
    def r2List = (read2 instanceof List) ? read2 : [read2]

    if( r1List.size() != r2List.size() ) {
        error "Sample ${sample_id}: read1 count (${r1List.size()}) != read2 count (${r2List.size()})"
    }

    if( r1List.isEmpty() ) {
        error "Sample ${sample_id}: no FASTQ pairs resolved from read1/read2"
    }

    def sortedR1 = r1List.sort { it.name }
    def sortedR2 = r2List.sort { it.name }
    def inFqArgs = sortedR1.indices.collect { i -> "--in-fq ${sortedR1[i]} ${sortedR2[i]}" }.join(' ')

    """
    SCRATCH_DIR="\${SLURM_LOCAL_SCRATCH:-}"
    if [[ -z "\$SCRATCH_DIR" ]]; then
      SCRATCH_DIR="\${SLURM_TMPDIR:-}"
    fi
    if [[ -z "\$SCRATCH_DIR" ]]; then
      SCRATCH_DIR="\${TMPDIR:-\$PWD/tmp}"
    fi
    SCRATCH_DIR="\$SCRATCH_DIR/pbgatk_${sample_id}"
    mkdir -p "\$SCRATCH_DIR"

    # Copy staged ref + index files into scratch so Parabricks finds them together
    cp -f ${ref}     \$SCRATCH_DIR/${ref.name}
    cp -f ${ref_amb} \$SCRATCH_DIR/${ref.name}.amb
    cp -f ${ref_ann} \$SCRATCH_DIR/${ref.name}.ann
    cp -f ${ref_bwt} \$SCRATCH_DIR/${ref.name}.bwt
    cp -f ${ref_pac} \$SCRATCH_DIR/${ref.name}.pac
    cp -f ${ref_sa}  \$SCRATCH_DIR/${ref.name}.sa
    cp -f ${ref_fai} \$SCRATCH_DIR/${ref.name}.fai

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
