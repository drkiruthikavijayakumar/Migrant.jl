using Pkg
Pkg.activate(".")
using Migrant

oral_fastq      = "path/to/public_oral.fastq"
tumor_fastq     = "path/to/public_tumor.fastq"
bacterial_fasta = "examples/f_nucleatum.fna"
bacterial_gff3  = "examples/f_nucleatum.gff"

println("Starting benchmark of analyze_migration...")

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
    println("\nCDS Profile Coverage:")
    for (coords, coverage) in result.cds_profile
        println("  CDS at ", coords, " => ", round(coverage * 100, digits=2), "%")
    end
    println("========================================")
catch e
    println("\n[Note] Benchmarking aborted: Please ensure the FASTQ paths point to valid sequencing files.")
    println("Error details: ", e)
end
