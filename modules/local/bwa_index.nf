// Builds BWA index + samtools .fai for a reference FASTA if not already present.
// Uses storeDir so the index is built only once and reused across runs.

process BWA_INDEX {
    label 'cpu_medium'
    container params.container_bwa

    storeDir "${params.outdir}/ref_index"

    input:
    path ref

    output:
    tuple path("${ref.name}"),
          path("${ref.name}.amb"),
          path("${ref.name}.ann"),
          path("${ref.name}.bwt"),
          path("${ref.name}.pac"),
          path("${ref.name}.sa"), emit: bwa_indexed

    script:
    """
    # ref is already staged by Nextflow; index it directly
    bwa index ${ref}
    """
}

process SAMTOOLS_FAIDX {
    label 'cpu_medium'
    container params.container_samtools

    storeDir "${params.outdir}/ref_index"

    input:
    path ref

    output:
    path("${ref.name}.fai"), emit: fai

    script:
    """
    # ref is already staged by Nextflow; index it directly
    samtools faidx ${ref}
    """
}

workflow BUILD_BWA_INDEX {
    take:
    ref

    main:
    BWA_INDEX(ref)
    SAMTOOLS_FAIDX(ref)

    ref_with_index = BWA_INDEX.out.bwa_indexed
        .combine(SAMTOOLS_FAIDX.out.fai)
        .map { fasta, amb, ann, bwt, pac, sa, fai ->
            tuple(fasta, amb, ann, bwt, pac, sa, fai)
        }

    emit:
    ref_with_index = ref_with_index
}
