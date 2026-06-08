nextflow.enable.dsl = 2

include { PBGATK_GERMLINE } from './modules/local/pbgatk'
include { COMPRESS_AND_INDEX_VCF } from './modules/local/compress_vcf'
include { GET_CONTIGS } from './modules/local/contigs'
include { EXTRACT_CONTIG_SAMPLE; MERGE_CONTIG_CHUNK; MERGE_CONTIG_FINAL } from './modules/local/merge_contigs'
include { FLAGSTAT_CRAM } from './modules/local/flagstat'
include { PRECHECK_GPU_PROFILE_COUNTS } from './modules/local/preflight'

workflow {

    // =========================
    // INPUT VALIDATION (MUST be here)
    // =========================
    if( !params.samplesheet )
        error "Missing required parameter: --samplesheet"

    if( !params.ref )
        error "Missing required parameter: --ref"

    samplesheet_path = file(params.samplesheet, checkIfExists: true)

    PRECHECK_GPU_PROFILE_COUNTS(samplesheet_path)

    ch_samples = Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true)
        .map { row ->
        
            def sampleId = row.sample_id?.toString()?.trim()
        
            def read1 = row.read1
                ?.toString()
                ?.trim()
                ?.split(';')
                ?.collect { file(it.trim(), checkIfExists: true) }
        
            def read2 = row.read2
                ?.toString()
                ?.trim()
                ?.split(';')
                ?.collect { file(it.trim(), checkIfExists: true) }
        
            def gpuProfile = row.gpu_profile?.toString()?.trim() ?: '1gpu'
        
            def validGpuProfiles = ['1gpu', '2gpu', '4gpu'] as Set
        
            if( !sampleId || !read1 || !read2 )
                error "Missing required fields in samplesheet"
        
            if( read1.size() != read2.size() )
                error "R1/R2 mismatch for ${sampleId}"
        
            if( !validGpuProfiles.contains(gpuProfile) )
                error "Invalid gpu_profile '${gpuProfile}'"
        
            tuple(sampleId, read1, read2, gpuProfile)
        }
     

    // =========================
    // PIPELINE
    // =========================

    pbgatk_out = PBGATK_GERMLINE(ch_samples, file(params.ref, checkIfExists: true))
    compressed_out = COMPRESS_AND_INDEX_VCF(pbgatk_out.vcf)

    contig_file = GET_CONTIGS(file(params.ref, checkIfExists: true))

    ch_contigs = contig_file.contigs
        .splitText()
        .map { it.trim() }
        .filter { it }

    if( params.contig_subset ) {
        def selected = params.contig_subset.toString().split(',') as Set
        ch_contigs = ch_contigs.filter { selected.contains(it) }
    }

    ch_contig_sample_vcfgz = ch_contigs
        .combine(compressed_out.vcfgz)
        .map { contig, sampleTuple -> tuple(contig, sampleTuple[0], sampleTuple[1], sampleTuple[2]) }

    extracted_out = EXTRACT_CONTIG_SAMPLE(ch_contig_sample_vcfgz)

    ch_chunk_inputs = extracted_out.contig_vcfgz
        .groupTuple(by: 0)
        .flatMap { contig, vcfgzList, csiList ->
            def sorted = vcfgzList.sort { a, b -> a.name <=> b.name }
            sorted.collate(params.merge_chunk_size).withIndex().collect { chunk, idx ->
                tuple(contig, idx + 1, chunk)
            }
        }

    chunk_out = MERGE_CONTIG_CHUNK(ch_chunk_inputs)

    ch_final_merge_inputs = chunk_out.chunk_vcfgz
        .groupTuple(by: 0)
        .map { contig, chunkVcfs, chunkCsis ->
            def sorted = chunkVcfs.sort { a, b -> a.name <=> b.name }
            tuple(contig, sorted)
        }

    final_out = MERGE_CONTIG_FINAL(ch_final_merge_inputs)

    if( params.run_flagstat ) {
        FLAGSTAT_CRAM(pbgatk_out.cram)
    }
}