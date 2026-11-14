using FASTX
using BioSequences

struct MigrationResult
    containment_score::Float64
    cds_profile::Dict{Tuple{Int, Int}, Float64}
end

function analyze_migration(source_fastq::String, sink_fastq::String, bacterial_fasta::String, bacterial_gff3::String; k::Int=31, target_threshold::Float64=0.05)
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
    
    clean_source = tempname() * ".fastq"
    clean_sink = tempname() * ".fastq"
    
    extract_target_reads(source_fastq, clean_source, target_idx, k, target_threshold)
    extract_target_reads(sink_fastq, clean_sink, target_idx, k, target_threshold)
    
    # Migration Tracking
    source_ncdna_kmers = build_ncdna_footprint(bacterial_fasta, ncdna_coords, k)
    score = calculate_kmer_containment(clean_sink, source_ncdna_kmers, k)
    
    # CDS Functional Profiling
    cds_idx = build_cds_index(bacterial_fasta, cds_coords, k)
    cds_prof = profile_cds(clean_sink, cds_idx, k)
    
    # Cleanup
    rm(clean_source, force=true)
    rm(clean_sink, force=true)
    
    return MigrationResult(score, cds_prof)
end
