module migrant

include("target_extract.jl")
include("gff_parser.jl")
include("strain_matcher.jl")
include("cds_tracker.jl")
include("migrant_core.jl")

export extract_kmers, build_target_index, extract_target_reads
export parse_cds_coordinates, extract_ncdna_regions
export build_ncdna_footprint, calculate_kmer_containment
export build_cds_index, profile_cds
export MigrationResult, analyze_migration

end
