process GET_CONTIGS {
    input:
    path ref

    output:
    path 'contigs.txt', emit: contigs

    script:
    """
    awk '{print \$1}' *.fai > contigs.txt
    """
}
