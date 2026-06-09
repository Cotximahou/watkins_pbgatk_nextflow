nextflow.enable.dsl = 2

include { PBGATK_GERMLINE } from './modules/local/pbgatk'
include { COMPRESS_AND_INDEX_VCF } from './modules/local/compress_vcf'
include { GET_CONTIGS } from './modules/local/contigs'
include { EXTRACT_CONTIG_SAMPLE; MERGE_CONTIG_CHUNK; MERGE_CONTIG_FINAL } from './modules/local/merge_contigs'
include { FLAGSTAT_CRAM } from './modules/local/flagstat'
include { PRECHECK_GPU_PROFILE_COUNTS } from './modules/local/preflight'
include { BUILD_BWA_INDEX } from './modules/local/bwa_index'

params.samplesheet = params.samplesheet ?: null
params.ref = params.ref ?: null
params.outdir = params.outdir ?: 'results'
params.contig_subset = params.contig_subset ?: ''
params.merge_chunk_size = (params.merge_chunk_size ?: 250) as int
params.run_flagstat = (params.run_flagstat ?: false) as boolean

workflow {
    if( !params.samplesheet ) {
        error "Missing required parameter: --samplesheet"
    }

    if( !params.ref ) {
        error "Missing required parameter: --ref"
    }

    samplesheet_path = file(params.samplesheet, checkIfExists: true)

    PRECHECK_GPU_PROFILE_COUNTS(samplesheet_path)

    ref_file = file(params.ref, checkIfExists: true)
    ref_indexed = BUILD_BWA_INDEX(ref_file)

    ch_samples = Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true)
        .map { row ->
            def sampleId   = row.sample_id?.toString()?.trim()
            def read1Raw   = row.read1?.toString()?.trim()
            def read2Raw   = row.read2?.toString()?.trim()
            def gpuProfile = row.gpu_profile ? row.gpu_profile.toString().trim() : '1gpu'
            def validGpuProfiles = ['1gpu', '2gpu', '4gpu'] as Set

            if( !sampleId || !read1Raw || !read2Raw ) {
                error "Each samplesheet row must contain sample_id, read1, and read2"
            }

            if( !validGpuProfiles.contains(gpuProfile) ) {
                error "Invalid gpu_profile '${gpuProfile}' for sample '${sampleId}'. Allowed values: 1gpu, 2gpu, 4gpu"
            }

            // Support semicolon-delimited multi-lane paths in read1/read2 columns
            def r1Files = read1Raw.split(';').collect { it.trim() }.findAll { it }.collect { file(it, checkIfExists: true) }
            def r2Files = read2Raw.split(';').collect { it.trim() }.findAll { it }.collect { file(it, checkIfExists: true) }

            if( r1Files.size() != r2Files.size() ) {
                error "Sample '${sampleId}': read1 has ${r1Files.size()} file(s) but read2 has ${r2Files.size()} file(s)"
            }

            tuple(sampleId, r1Files, r2Files, gpuProfile)
        }

    pbgatk_out = PBGATK_GERMLINE(ch_samples, ref_indexed.ref_with_index)
    compressed_out = COMPRESS_AND_INDEX_VCF(pbgatk_out.vcf)

    contig_file = GET_CONTIGS(ref_indexed.ref_with_index.map { it[0] })

    ch_contigs = contig_file.contigs
        .splitText()
        .map { it.trim() }
        .filter { it }

    if( params.contig_subset ) {
        def selected = params.contig_subset
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .toSet()
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
