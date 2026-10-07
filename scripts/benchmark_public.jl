using Pkg
Pkg.activate(".")
using migrant
using Downloads

# -------------------------------------------------------------------------
# POINT THESE VARIABLES TO DOWNLOADED PUBLIC SRA/NCBI FILES
# Note: You can use Downloads.download(url, dest) to fetch files programmatically.
# -------------------------------------------------------------------------
oral_fastq      = "path/to/public_oral.fastq"
tumor_fastq     = "path/to/public_tumor.fastq"
bacterial_fasta = "path/to/fusobacterium_nucleatum.fasta"
bacterial_gff3  = "path/to/fusobacterium_nucleatum.gff3"

println("Starting benchmark of analyze_migration...")

# Wrap the function execution in @time to benchmark CPU time and RAM allocation
try
    @time result = analyze_migration(
        oral_fastq, 
        tumor_fastq, 
        bacterial_fasta, 
        bacterial_gff3; 
        k=31, 
        target_threshold=0.05
    )

    println("\n========================================")
    println("      MIGRATION RESULT SUMMARY")
    println("========================================")
    println("Containment Score (ncDNA): ", round(result.containment_score * 100, digits=2), "%")
    println("\nVirulence Profile Coverage:")
    for (coords, coverage) in result.virulence_profile
        println("  Gene at ", coords, " => ", round(coverage * 100, digits=2), "%")
    end
    println("========================================")
catch e
    println("\n[Note] Benchmarking aborted: Please ensure the placeholder paths point to valid valid biological files.")
    println("Error details: ", e)
end
