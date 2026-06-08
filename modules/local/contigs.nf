process GET_CONTIGS {
    input:
    path ref

    output:
    path 'contigs.txt', emit: contigs

    script:
    """
    grep '^>' ${ref} | sed 's/^>//' > contigs.txt
    """
}
