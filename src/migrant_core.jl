using FASTX
using BioSequences

struct MigrationResult
    containment_score::Float64
    virulence_profile::Dict{Tuple{Int, Int}, Float64}
end

function analyze_migration(oral_fastq::String, tumor_fastq::String, host_fasta::String, bacterial_fasta::String, bacterial_gff3::String; k::Int=31, host_threshold::Float64=0.1)
    # Parse GFF3
    cds_coords = parse_cds_coordinates(bacterial_gff3)
    
    # Extract ncDNA coordinates
    reader = FASTA.Reader(open(bacterial_fasta, "r"))
    record = first(reader)
    close(reader)
    genome_len = length(FASTA.sequence(LongDNA{4}, record))
    ncdna_coords = extract_ncdna_regions(cds_coords, genome_len)
    
    # Host Depletion
    host_idx = build_host_index(host_fasta, k)
    
    clean_oral = "temp_oral.fastq"
    clean_tumor = "temp_tumor.fastq"
    
    filter_host_reads(oral_fastq, clean_oral, host_idx, k, host_threshold)
    filter_host_reads(tumor_fastq, clean_tumor, host_idx, k, host_threshold)
    
    # Migration Tracking
    oral_ncdna_kmers = build_ncdna_footprint(bacterial_fasta, ncdna_coords, k)
    score = calculate_kmer_containment(clean_tumor, oral_ncdna_kmers, k)
    
    # Virulence Profiling
    vir_idx = build_virulence_index(bacterial_fasta, cds_coords, k)
    v_profile = profile_virulence(clean_tumor, vir_idx, k)
    
    # Cleanup
    rm(clean_oral, force=true)
    rm(clean_tumor, force=true)
    
    return MigrationResult(score, v_profile)
end
