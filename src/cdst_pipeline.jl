function extract_features(result::MigrationResult, sample_id::String, class_label::Int)
    df = DataFrame(sample_id = sample_id, class_label = class_label, containment_score = result.containment_score)
    for (coords, coverage) in result.virulence_profile
        col_name = "Gene_$(coords[1])_$(coords[2])"
        df[!, col_name] = [coverage]
    end
    return df
end

function run_cdst_pipeline(sample_id::String, oral_fastq::String, tumor_fastq::String, host_fasta::String, bacterial_fasta::String, bacterial_gff3::String, trained_mach::Machine, report_out::String; k::Int=31, host_threshold::Float64=0.1)
    # Step A: Analyze migration
    result = analyze_migration(oral_fastq, tumor_fastq, host_fasta, bacterial_fasta, bacterial_gff3; k=k, host_threshold=host_threshold)
    
    # Step B: Extract features into DataFrame
    df = extract_features(result, sample_id, 0)
    
    # Step C: Predict risk
    risk_score = predict_tropism(trained_mach, df)
    
    # Step D: Generate report
    generate_clinical_report(sample_id, result, risk_score, report_out)
    
    return risk_score
end
