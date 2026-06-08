nextflow.enable.dsl = 2

include { PBGATK_GERMLINE } from './modules/local/pbgatk'
include { COMPRESS_AND_INDEX_VCF } from './modules/local/compress_vcf'
include { GET_CONTIGS } from './modules/local/contigs'
include { EXTRACT_CONTIG_SAMPLE; MERGE_CONTIG_CHUNK; MERGE_CONTIG_FINAL } from './modules/local/merge_contigs'
include { FLAGSTAT_CRAM } from './modules/local/flagstat'
include { PRECHECK_GPU_PROFILE_COUNTS } from './modules/local/preflight'

params.samplesheet = params.samplesheet ?: null
params.ref = params.ref ?: null
params.outdir = params.outdir ?: 'results'
params.contig_subset = params.contig_subset ?: ''
params.merge_chunk_size = (params.merge_chunk_size ?: 250) as int
params.run_flagstat = (params.run_flagstat ?: false) as boolean

if( !params.samplesheet ) error "Missing required parameter: --samplesheet"
if( !params.ref ) error "Missing required parameter: --ref"

workflow {

    samplesheet_path = file(params.samplesheet, checkIfExists: true)

    // Optional preflight (can disable if debugging)
    PRECHECK_GPU_PROFILE_COUNTS(samplesheet_path)

    ch_samples = Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true)
        .map { row ->

            def sampleId = row.sample_id?.toString()?.trim()
            def gpuProfile = row.gpu_profile ? row.gpu_profile.toString().trim() : '1gpu'

            def validGpuProfiles = ['1gpu', '2gpu', '4gpu'] as Set

            if( !sampleId ) {
                error "Missing sample_id"
            }

            if( !validGpuProfiles.contains(gpuProfile) ) {
                error "Invalid gpu_profile '${gpuProfile}' for sample '${sampleId}'"
            }

            def read1 = row.read1?.toString()?.trim()
                ?.split(';')
                ?.findAll { it }
                ?.collect { file(it.trim(), checkIfExists: true) }

            def read2 = row.read2?.toString()?.trim()
                ?.split(';')
                ?.findAll { it }
                ?.collect { file(it.trim(), checkIfExists: true) }

            if( !read1 || !read2 ) {
                error "Missing FASTQs for sample ${sampleId}"
            }

            tuple(sampleId, read1, read2, gpuProfile)
        }

    // Parabricks germline calling
    pbgatk_out = PBGATK_GERMLINE(ch_samples, file(params.ref, checkIfExists: true))

    compressed_out = COMPRESS_AND_INDEX_VCF(pbgatk_out.vcf)

    // Contig handling
    contig_file = GET_CONTIGS(file(params.ref, checkIfExists: true))

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
        .map { contig, sampleTuple ->
            tuple(contig, sampleTuple[0], sampleTuple[1], sampleTuple[2])
        }

    extracted_out = EXTRACT_CONTIG_SAMPLE(ch_contig_sample_vcfgz)

    ch_chunk_inputs = extracted_out.contig_vcfgz
        .groupTuple(by: 0)
        .flatMap { contig, vcfgzList, csiList ->

            def sorted = vcfgzList.sort { a, b -> a.name <=> b.name }

            sorted.collate(params.merge_chunk_size)
                .withIndex()
                .collect { chunk, idx ->
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