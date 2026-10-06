using Pkg
Pkg.activate(".")
using migrant
using DataFrames
using CSV

# Mock function for feature extraction
function extract_features(result::MigrationResult, patient_id::String, class_label::Int)
    # Flattens the migration result into a single DataFrame row
    df = DataFrame(PatientID = patient_id, Class = class_label, ContainmentScore = result.containment_score)
    
    # Add virulence profile as separate columns
    for (coords, coverage) in result.virulence_profile
        col_name = "Gene_$(coords[1])_$(coords[2])"
        df[!, col_name] = [coverage]
    end
    
    return df
end

# Mock function to export the matrix
function export_feature_matrix(features::Vector{DataFrame}, output_csv::String)
    if isempty(features)
        println("No features to export.")
        return
    end
    # Combine all patient records. We use cols=:union in case certain patients 
    # were mapped against slightly varying profiles.
    combined_df = vcat(features..., cols=:union)
    
    # Write the combined feature matrix to CSV
    CSV.write(output_csv, combined_df)
    println("Successfully exported full cohort feature matrix to: ", output_csv)
end

function process_cohort(data_dir::String, host_ref::String, bacterial_ref::String, gff3_ref::String, output_csv::String)
    features_vec = Vector{DataFrame}()
    
    if !isdir(data_dir)
        println("Data directory not found: ", data_dir)
        return
    end
    
    patient_folders = filter(p -> isdir(joinpath(data_dir, p)), readdir(data_dir))
    
    for patient_id in patient_folders
        patient_path = joinpath(data_dir, patient_id)
        oral_fastq = joinpath(patient_path, "oral.fastq")
        tumor_fastq = joinpath(patient_path, "tumor.fastq")
        
        if isfile(oral_fastq) && isfile(tumor_fastq)
            println("Processing Patient: ", patient_id)
            
            try
                result = analyze_migration(
                    oral_fastq, 
                    tumor_fastq, 
                    host_ref, 
                    bacterial_ref, 
                    gff3_ref; 
                    k=31, 
                    host_threshold=0.1
                )
                
                # Assume class_label = 1 for matched tumor samples
                class_label = 1
                
                df_row = extract_features(result, patient_id, class_label)
                push!(features_vec, df_row)
                println("  -> Successfully analyzed and extracted features.")
            catch e
                println("  -> Error processing patient ", patient_id, ": ", e)
            end
        else
            println("Skipping ", patient_id, ": Missing oral.fastq or tumor.fastq")
        end
    end
    
    export_feature_matrix(features_vec, output_csv)
end

# Example trigger (commented out until data exists):
# process_cohort("data/clinical_cohort", "data/ref/host.fasta", "data/ref/bacteria.fasta", "data/ref/bacteria.gff3", "cohort_features.csv")
