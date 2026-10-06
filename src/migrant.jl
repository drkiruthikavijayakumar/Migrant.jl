module migrant

include("host_deplete.jl")
include("gff_parser.jl")
include("strain_matcher.jl")
include("virulence_tracker.jl")
include("migrant_core.jl")

export extract_kmers, build_host_index, filter_host_reads
export parse_cds_coordinates, extract_ncdna_regions
export build_ncdna_footprint, calculate_kmer_containment
export build_virulence_index, profile_virulence
export MigrationResult, analyze_migration

end
