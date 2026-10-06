using Test
using migrant
using FASTX
using BioSequences

@testset "migrant.jl tests" begin
    # 1. Generate dummy host .fasta
    host_fasta = tempname() * ".fasta"
    open(host_fasta, "w") do io
        write(io, ">host_seq\n")
        write(io, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
    end
    
    # 2. Generate dummy .fastq with 2 human-like and 2 bacterial-like reads
    reads_fastq = tempname() * ".fastq"
    open(reads_fastq, "w") do io
        # human 1
        write(io, "@human1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        # human 2
        write(io, "@human2\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        # bac 1
        write(io, "@bac1\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        # bac 2
        write(io, "@bac2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    end
    
    out_fastq = tempname() * ".fastq"
    
    k = 15
    threshold = 0.5
    
    # Run host depletion
    host_index = build_host_index(host_fasta, k)
    filter_host_reads(reads_fastq, out_fastq, host_index, k, threshold)
    
    # Assert output contains exactly 2 reads
    # And specifically they should be bac1 and bac2
    out_reads = []
    reader = FASTQ.Reader(open(out_fastq, "r"))
    for record in reader
        push!(out_reads, FASTQ.identifier(record))
    end
    close(reader)
    
    @test length(out_reads) == 2
    @test "bac1" in out_reads
    @test "bac2" in out_reads
    
    # Cleanup
    rm(host_fasta, force=true)
    rm(reads_fastq, force=true)
    rm(out_fastq, force=true)
end

@testset "GFF Parser" begin
    # Generate dummy .gff3 file
    gff_path = tempname() * ".gff3"
    open(gff_path, "w") do io
        write(io, "##gff-version 3\n")
        write(io, "seq1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene1\n")
        write(io, "seq1\ttest\tCDS\t100\t200\t.\t+\t0\tID=cds1;Parent=gene1\n")
        write(io, "seq1\ttest\tmRNA\t100\t600\t.\t+\t.\tID=mrna1;Parent=gene1\n")
        write(io, "seq1\ttest\tCDS\t500\t600\t.\t+\t0\tID=cds2;Parent=gene1\n")
        write(io, "\n") # empty line
        write(io, "# another comment\n")
    end
    
    # Test parse_cds_coordinates
    coords = parse_cds_coordinates(gff_path)
    @test coords == [(100, 200), (500, 600)]
    
    # Test extract_ncdna_regions
    genome_length = 1000
    ncdna = extract_ncdna_regions(coords, genome_length)
    @test ncdna == [(1, 99), (201, 499), (601, 1000)]
    
    # Cleanup
    rm(gff_path, force=true)
end

@testset "Strain Matcher" begin
    k = 15
    
    # 1. Create dummy bacterial .fasta file
    bac_fasta = tempname() * ".fasta"
    open(bac_fasta, "w") do io
        write(io, ">bac1\n")
        # 100 bp long genome
        write(io, "A"^50 * "C"^50 * "\n")
    end
    
    # Non-coding coordinates: (1, 50) - this corresponds to "A"^50
    coords = [(1, 50)]
    
    # 2. Build ncDNA footprint
    oral_kmers = build_ncdna_footprint(bac_fasta, coords, k)
    
    # 3. Create dummy tumor .fastq containing reads perfectly matching "A"^50
    tumor_fastq = tempname() * ".fastq"
    open(tumor_fastq, "w") do io
        write(io, "@read1\n")
        write(io, "A"^50 * "\n")
        write(io, "+\n")
        write(io, "!"^50 * "\n")
    end
    
    # 4. Calculate k-mer containment
    containment = calculate_kmer_containment(tumor_fastq, oral_kmers, k)
    
    @test containment == 1.0
    
    # Cleanup
    rm(bac_fasta, force=true)
    rm(tumor_fastq, force=true)
end

@testset "Virulence Tracker" begin
    k = 15
    
    # 1. Create dummy bacterial .fasta file
    bac_fasta = tempname() * ".fasta"
    open(bac_fasta, "w") do io
        write(io, ">bac1\n")
        # 100 bp long genome
        write(io, "G"^50 * "T"^50 * "\n")
    end
    
    # CDS coordinates: (1, 50) - corresponding to "G"^50
    cds_coords = [(1, 50)]
    
    # 2. Build virulence index
    virulence_index = build_virulence_index(bac_fasta, cds_coords, k)
    
    # 3. Create dummy tumor .fastq containing reads perfectly matching "G"^50
    tumor_fastq = tempname() * ".fastq"
    open(tumor_fastq, "w") do io
        write(io, "@read1\n")
        write(io, "G"^50 * "\n")
        write(io, "+\n")
        write(io, "!"^50 * "\n")
    end
    
    # 4. Profile virulence
    results = profile_virulence(tumor_fastq, virulence_index, k)
    
    @test haskey(results, (1, 50))
    @test results[(1, 50)] == 1.0
    
    # Cleanup
    rm(bac_fasta, force=true)
    rm(tumor_fastq, force=true)
end

@testset "End-to-End Core" begin
    k = 15
    
    host_fasta = tempname() * ".fasta"
    open(host_fasta, "w") do io
        write(io, ">host\n", "A"^100, "\n")
    end
    
    bac_fasta = tempname() * ".fasta"
    open(bac_fasta, "w") do io
        write(io, ">bac\n", "C"^50 * "G"^50, "\n")
    end
    
    bac_gff = tempname() * ".gff3"
    open(bac_gff, "w") do io
        write(io, "seq\t.\tCDS\t1\t50\t.\t+\t.\tID=cds1\n")
    end
    
    oral_fastq = tempname() * ".fastq"
    open(oral_fastq, "w") do io
        write(io, "@o1\n", "C"^50, "\n+\n", "!"^50, "\n")
    end
    
    tumor_fastq = tempname() * ".fastq"
    open(tumor_fastq, "w") do io
        write(io, "@t1\n", "G"^50, "\n+\n", "!"^50, "\n")
    end
    
    res = analyze_migration(oral_fastq, tumor_fastq, host_fasta, bac_fasta, bac_gff; k=15, host_threshold=0.1)
    
    @test res.containment_score isa Float64
    @test res.virulence_profile isa Dict
    
    rm(host_fasta, force=true)
    rm(bac_fasta, force=true)
    rm(bac_gff, force=true)
    rm(oral_fastq, force=true)
    rm(tumor_fastq, force=true)
end

@testset "ML Predictor" begin
    using DataFrames
    
    # 1. Create dummy DataFrame
    dummy_df = DataFrame(
        sample_id = ["S1", "S2", "S3", "S4"],
        class_label = [0, 1, 0, 1],
        ncdna_containment = [0.1, 0.8, 0.2, 0.9],
        virulence_1_50 = [0.0, 1.0, 0.1, 0.85]
    )
    
    # 2. Train Tropism Model
    mach = train_tropism_model(dummy_df)
    
    # 3. Simulate new patient
    new_patient = DataFrame(
        sample_id = ["S1"],
        class_label = [0],
        ncdna_containment = [0.1],
        virulence_1_50 = [0.0]
    )
    
    # 4. Predict
    score = predict_tropism(mach, new_patient)
    
    # 5. Assertions
    @test score isa Float64
    @test 0.0 <= score <= 1.0
end

@testset "Clinical Reporting" begin
    dummy_result = MigrationResult(0.95, Dict((1, 50) => 0.88, (100, 200) => 0.45))
    output_html = "test_report.html"
    
    success = generate_clinical_report("Test_001", dummy_result, 0.85, output_html)
    
    @test success == true
    @test isfile(output_html)
    
    # Cleanup
    rm(output_html, force=true)
end
