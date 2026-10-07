using FASTX
using BioSequences

struct MigrationResult
    containment_score::Float64
    virulence_profile::Dict{Tuple{Int, Int}, Float64}
end

function analyze_migration(oral_fastq::String, tumor_fastq::String, bacterial_fasta::String, bacterial_gff3::String; k::Int=31, target_threshold::Float64=0.05)
    # Parse GFF3
    cds_coords = parse_cds_coordinates(bacterial_gff3)
    
    # Extract ncDNA coordinates
    reader = FASTA.Reader(open(bacterial_fasta, "r"))
    record = first(reader)
    close(reader)
    genome_len = length(FASTA.sequence(String, record))
    ncdna_coords = extract_ncdna_regions(cds_coords, genome_len)
    
    # Target Extraction
    target_idx = build_target_index(bacterial_fasta, k)
    
    clean_oral = tempname() * ".fastq"
    clean_tumor = tempname() * ".fastq"
    
    extract_target_reads(oral_fastq, clean_oral, target_idx, k, target_threshold)
    extract_target_reads(tumor_fastq, clean_tumor, target_idx, k, target_threshold)
    
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
