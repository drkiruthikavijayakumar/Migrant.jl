using Test
using migrant
using FASTX
using BioSequences

@testset "Target Extraction" begin
    # 1. Generate dummy target .fasta
    target_fasta = tempname() * ".fasta"
    open(target_fasta, "w") do io
        write(io, ">target_seq\n")
        write(io, "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")
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
    
    # Run target extraction
    target_index = build_target_index(target_fasta, k)
    extract_target_reads(reads_fastq, out_fastq, target_index, k, threshold)
    
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
    rm(target_fasta, force=true)
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
    source_kmers = build_ncdna_footprint(bac_fasta, coords, k)
    
    # 3. Create dummy sink .fastq containing reads perfectly matching "A"^50
    sink_fastq = tempname() * ".fastq"
    open(sink_fastq, "w") do io
        write(io, "@read1\n")
        write(io, "A"^50 * "\n")
        write(io, "+\n")
        write(io, "!"^50 * "\n")
    end
    
    # 4. Calculate k-mer containment
    containment = calculate_kmer_containment(sink_fastq, source_kmers, k)
    
    @test containment == 1.0
    
    # Cleanup
    rm(bac_fasta, force=true)
    rm(sink_fastq, force=true)
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
    
    # 3. Create dummy sink .fastq containing reads perfectly matching "G"^50
    sink_fastq = tempname() * ".fastq"
    open(sink_fastq, "w") do io
        write(io, "@read1\n")
        write(io, "G"^50 * "\n")
        write(io, "+\n")
        write(io, "!"^50 * "\n")
    end
    
    # 4. Profile virulence
    results = profile_virulence(sink_fastq, virulence_index, k)
    
    @test haskey(results, (1, 50))
    @test results[(1, 50)] == 1.0
    
    # Cleanup
    rm(bac_fasta, force=true)
    rm(sink_fastq, force=true)
end

@testset "End-to-End Core" begin
    k = 15
    
    bac_fasta = tempname() * ".fasta"
    open(bac_fasta, "w") do io
        write(io, ">bac\n", "C"^50 * "G"^50, "\n")
    end
    
    bac_gff = tempname() * ".gff3"
    open(bac_gff, "w") do io
        write(io, "seq\t.\tCDS\t1\t50\t.\t+\t.\tID=cds1\n")
    end
    
    source_fastq = tempname() * ".fastq"
    open(source_fastq, "w") do io
        write(io, "@o1\n", "C"^50, "\n+\n", "!"^50, "\n")
    end
    
    sink_fastq = tempname() * ".fastq"
    open(sink_fastq, "w") do io
        write(io, "@t1\n", "G"^50, "\n+\n", "!"^50, "\n")
    end
    
    res = analyze_migration(source_fastq, sink_fastq, bac_fasta, bac_gff; k=15, target_threshold=0.05)
    
    @test res.containment_score isa Float64
    @test res.virulence_profile isa Dict
    
    rm(bac_fasta, force=true)
    rm(bac_gff, force=true)
    rm(source_fastq, force=true)
    rm(sink_fastq, force=true)
end
