process PRECHECK_GPU_PROFILE_COUNTS {

    label 'cpu_small'

    publishDir "${params.outdir}/preflight", mode: 'copy'

    input:
    path samplesheet

    output:
    path 'gpu_profile_counts.tsv', emit: counts_tsv
    path 'gpu_profile_summary.txt', emit: summary

    script:
    """
    awk -F',' '
    NR==1 {
        idx=0
        for(i=1;i<=NF;i++) {
            col=\$i
            gsub(/^ +| +\$/, "", col)

            if(tolower(col)=="gpu_profile")
                idx=i
        }

        if(idx==0) {
            print "ERROR: samplesheet is missing gpu_profile column" > "/dev/stderr"
            exit 1
        }
        next
    }

    {
        g=\$idx
        gsub(/^ +| +\$/, "", g)

        if(g=="")
            g="1gpu"

        if(g!="1gpu" && g!="2gpu" && g!="4gpu") {
            printf("ERROR: invalid gpu_profile %s on CSV line %d\\n", g, NR) > "/dev/stderr"
            bad=1
        }

        c[g]++
    }

    END {
        if(bad==1)
            exit 1

        print "gpu_profile,count" > "gpu_profile_counts.tsv"

        for(k in c)
            print k "," c[k] > "gpu_profile_counts.unsorted.csv"
    }
    ' ${samplesheet}

    sort gpu_profile_counts.unsorted.csv > gpu_profile_counts.sorted.csv

    {
        echo "gpu_profile,count"
        cat gpu_profile_counts.sorted.csv
    } > gpu_profile_counts.tsv

    cp gpu_profile_counts.tsv gpu_profile_summary.txt

    echo "Preflight GPU profile summary:"
    cat gpu_profile_summary.txt

    rm -f gpu_profile_counts.unsorted.csv gpu_profile_counts.sorted.csv
    """
}