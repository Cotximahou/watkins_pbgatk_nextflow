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

    if (!params.samplesheet) error "Missing --samplesheet"
    if (!params.ref) error "Missing --ref"

    samplesheet_path = file(params.samplesheet, checkIfExists: true)

    PRECHECK_GPU_PROFILE_COUNTS(samplesheet_path)

    ref_file = file(params.ref, checkIfExists: true)
    ref_indexed = BUILD_BWA_INDEX(ref_file)

    ref_fa = ref_indexed.ref_with_index.first()

    // ---------------------------
    // SAMPLESHEET CHANNEL
    // ---------------------------

    ch_samples = Channel
        .fromPath(samplesheet_path)
        .splitCsv(header: true)
        .map { row ->

            def sampleId = row.sample_id?.trim()
            def read1Raw = row.read1?.trim()
            def read2Raw = row.read2?.trim()
            def gpuProfile = row.gpu_profile?.trim() ?: '1gpu'

            if (!sampleId || !read1Raw || !read2Raw)
                error "Missing required fields in samplesheet"

            def r1 = read1Raw.split(';').collect { file(it.trim(), checkIfExists: true) }
            def r2 = read2Raw.split(';').collect { file(it.trim(), checkIfExists: true) }

            if (r1.size() != r2.size())
                error "Lane mismatch for ${sampleId}"

            tuple(sampleId, r1, r2, gpuProfile)
        }

    // DEBUG
    ch_samples.view { "DEBUG CH_SAMPLES = ${it}" }

    // ---------------------------
    // GPU STEP
    // ---------------------------

    pbgatk_out = PBGATK_GERMLINE(ch_samples, ref_fa)

    // DEBUG
    pbgatk_out.vcf.view { "DEBUG PBGATK_VCF = ${it}" }

    compressed_out = COMPRESS_AND_INDEX_VCF(pbgatk_out.vcf)

    // DEBUG
    compressed_out.vcfgz.view { "DEBUG COMPRESSED = ${it}" }

    // ---------------------------
    // CONTIG EXTRACTION
    // ---------------------------

    contig_file = GET_CONTIGS(ref_fa)

    ch_contigs = contig_file.contigs
        .splitText()
        .map { it.trim() }
        .filter { it }

    if (params.contig_subset) {
        def allowed = params.contig_subset.split(/\s*,\s*/).toSet()
        ch_contigs = ch_contigs.filter { allowed.contains(it) }
    }

    // DEBUG
    ch_contigs.view { "DEBUG CONTIG = ${it}" }

    // ---------------------------
    // SAFE CONTIG-SAMPLE CROSS
    // ---------------------------

    ch_contig_sample =
        ch_contigs
            .combine(compressed_out.vcfgz)
            .map { contig, sample_id, vcfgz, csi ->
                tuple(contig, sample_id, vcfgz, csi)
            }

    // DEBUG
    ch_contig_sample.view { "DEBUG CONTIG_SAMPLE = ${it}" }

    extracted = EXTRACT_CONTIG_SAMPLE(ch_contig_sample)

    // DEBUG
    extracted.contig_vcfgz.view { "DEBUG EXTRACTED = ${it}" }

    // ---------------------------
    // CHUNK MERGING
    // ---------------------------

    ch_chunks = extracted.contig_vcfgz
        .groupTuple()
        .flatMap { contig, vcfList, csiList ->

            // DEBUG
            println "DEBUG GROUPED CONTIG=${contig}"
            println "DEBUG VCF COUNT=${vcfList.size()}"
            println "DEBUG CSI COUNT=${csiList.size()}"

            def pairs = [vcfList, csiList].transpose()

            pairs.collate(params.merge_chunk_size)
                .withIndex()
                .collect { chunk, idx ->

                    def vcfs = chunk.collect { it[0] }
                    def csis = chunk.collect { it[1] }

                    tuple(
                        contig,
                        idx + 1,
                        vcfs.size(),
                        vcfs,
                        csis
                    )
                }
        }

    // DEBUG
    ch_chunks.view { "DEBUG CHUNK = ${it}" }

    chunk_out = MERGE_CONTIG_CHUNK(ch_chunks)

    // DEBUG
    chunk_out.chunk_vcfgz.view { "DEBUG CHUNK_OUT = ${it}" }

    final_in = chunk_out.chunk_vcfgz
        .groupTuple()
        .map { contig, vcfs, csis ->

            // DEBUG
            println "DEBUG FINAL_INPUT CONTIG=${contig}"
            println "DEBUG FINAL_INPUT N_VCFS=${vcfs.size()}"

            tuple(contig, vcfs.sort { it.name })
        }

    // DEBUG
    final_in.view { "DEBUG FINAL_IN = ${it}" }

    final_out = MERGE_CONTIG_FINAL(final_in)

    // DEBUG
    final_out.contig_vcf.view { "DEBUG FINAL_OUT = ${it}" }

    if (params.run_flagstat)
        FLAGSTAT_CRAM(pbgatk_out.cram)
}