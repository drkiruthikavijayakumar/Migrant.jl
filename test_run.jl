using Pkg
Pkg.activate(".")
#Pkg.activate("migrant.jl") # Activate your core package environment
using migrant

# Define the paths to the files we just prepared
oral_fastq  = "combined_oral.fastq"
tumor_fastq = "combined_stool.fastq" # Stool acting as the tumor proxy
host_fasta  = "human_reference.fasta"
bact_fasta  = "f_nucleatum.fna"
bact_gff3   = "f_nucleatum.gff"

println("Starting MIGRANT.jl on CRC Paired WGS Data...")
println("This may take some time depending on file size. Initiating Host Depletion and K-mer Mapping...")

# Use @time to see exactly how fast your Julia package is!
@time result = analyze_migration(
    oral_fastq, 
    tumor_fastq, 
    host_fasta, 
    bact_fasta, 
    bact_gff3;
    k = 31,
    host_threshold = 0.1
)

println("\n--- TEST RESULTS ---")
println("Strain Migration Confidence (ncDNA Containment): ", round(result.containment_score * 100, digits=2), "%")

println("\nVirulence Coverage (CDS Profiling):")
for (coords, coverage) in result.virulence_profile
    println("Gene at $coords: ", round(coverage * 100, digits=2), "%")
end
