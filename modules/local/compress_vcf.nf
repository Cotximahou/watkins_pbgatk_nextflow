process COMPRESS_AND_INDEX_VCF {
    publishDir "${params.outdir}/vcfgz", mode: 'copy', pattern: '*.vcf.gz'
    publishDir "${params.outdir}/vcfgz", mode: 'copy', pattern: '*.vcf.gz.csi'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi"), emit: vcfgz

    script:
    """
    bgzip --threads ${task.cpus} -c ${vcf} > ${sample_id}.vcf.gz
    bcftools index --threads ${task.cpus} ${sample_id}.vcf.gz
    """
}
