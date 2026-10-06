module migrant

include("host_deplete.jl")
include("gff_parser.jl")
include("strain_matcher.jl")
include("virulence_tracker.jl")
include("migrant_core.jl")
include("ml_predictor.jl")

include("report_generator.jl")
include("cdst_pipeline.jl")

export extract_kmers, build_host_index, filter_host_reads
export parse_cds_coordinates, extract_ncdna_regions
export build_ncdna_footprint, calculate_kmer_containment
export build_virulence_index, profile_virulence
export MigrationResult, analyze_migration
export train_tropism_model, predict_tropism
export generate_clinical_report, run_cdst_pipeline, extract_features

end
