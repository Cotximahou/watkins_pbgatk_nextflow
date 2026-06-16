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

    // ✅ FIX 1: extract real FASTA channel
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

    // ---------------------------
    // GPU STEP
    // ---------------------------
    pbgatk_out = PBGATK_GERMLINE(ch_samples, ref_fa)

    compressed_out = COMPRESS_AND_INDEX_VCF(pbgatk_out.vcf)

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

    // ---------------------------
    // SAFE CONTIG-SAMPLE CROSS
    // ---------------------------
    
    ch_contig_sample =
        ch_contigs
            .combine(compressed_out.vcfgz)
            .map { contig, sample_id, vcfgz, csi ->
                tuple(contig, sample_id, vcfgz, csi)
            }

    extracted = EXTRACT_CONTIG_SAMPLE(ch_contig_sample)

    // ---------------------------
    // CHUNK MERGING
    // ---------------------------
    ch_chunks = extracted.contig_vcfgz
        .map { contig, vcf, csi ->
            tuple(contig, vcf)
        }
        .groupTuple()
        .flatMap { contig, vcfList ->
    
            vcfList.collate(params.merge_chunk_size)
                .withIndex()
                .collect { chunk, idx ->
                    tuple(contig, idx + 1, chunk.size(), chunk)
                }
        }

    chunk_out = MERGE_CONTIG_CHUNK(ch_chunks)

    final_in = chunk_out.chunk_vcfgz
        .groupTuple()
        .map { contig, vcfs, csis ->
            tuple(contig, vcfs.sort { it.name })
        }

    final_out = MERGE_CONTIG_FINAL(final_in)

    if (params.run_flagstat)
        FLAGSTAT_CRAM(pbgatk_out.cram)
}