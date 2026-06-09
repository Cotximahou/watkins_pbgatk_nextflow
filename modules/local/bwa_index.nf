// Builds BWA index + samtools .fai for a reference FASTA if not already present.
// Uses storeDir so the index is built only once and reused across runs.

process BUILD_BWA_INDEX {
    label 'cpu_medium'
    container params.container_parabricks

    storeDir "${params.outdir}/ref_index"

    input:
    path ref

    output:
    tuple path("${ref.name}"),
          path("${ref.name}.amb"),
          path("${ref.name}.ann"),
          path("${ref.name}.bwt"),
          path("${ref.name}.pac"),
          path("${ref.name}.sa"),
          path("${ref.name}.fai"), emit: ref_with_index

    script:
    """
    # Copy ref into work dir for indexing
    cp -f ${ref} ${ref.name}

    # Build BWA index
    bwa index ${ref.name}

    # Build FASTA index
    samtools faidx ${ref.name}
    """
}
