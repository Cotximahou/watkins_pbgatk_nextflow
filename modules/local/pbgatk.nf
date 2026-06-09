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
    path ref

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

    REF_SRC="$(readlink -f ${ref})"
    REF_NAME="$(basename ${ref})"
    REF_LOCAL="\$SCRATCH_DIR/\$REF_NAME"

    cp -f "\$REF_SRC" "\$REF_LOCAL"

    REQUIRED_BWA_EXTS=(amb ann bwt pac sa)
    for ext in "\${REQUIRED_BWA_EXTS[@]}"; do
      if [[ ! -r "\$REF_SRC.\$ext" ]]; then
        echo "ERROR: Missing or unreadable reference index file: \$REF_SRC.\$ext" >&2
        exit 101
      fi
      cp -f "\$REF_SRC.\$ext" "\$REF_LOCAL.\$ext"
    done

    if [[ -r "\$REF_SRC.fai" ]]; then
      cp -f "\$REF_SRC.fai" "\$REF_LOCAL.fai"
    fi

    pbrun germline \
      --ref \$REF_LOCAL \
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
