process FLAGSTAT_CRAM {
    publishDir "${params.outdir}/flagstat", mode: 'copy', pattern: '*.flagstat.txt'

    input:
    tuple val(sample_id), path(cram)

    output:
    tuple val(sample_id), path("${sample_id}.flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat -@ ${task.cpus} ${cram} > ${sample_id}.flagstat.txt
    """
}
