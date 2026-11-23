using FASTX
using BioSequences
using Kmers

"""
    build_cds_index(fasta_path::String, cds_coords::Array{Tuple{Int, Int}}, k::Int) -> Dict{Tuple{Int, Int}, Set{DNAKmer{k}}}

Build a dictionary mapping each `(start, end)` CDS coordinate pair in `cds_coords` to its set of canonical `k`-mers from `fasta_path`.
"""
function build_cds_index(fasta_path::String, cds_coords::Array{Tuple{Int, Int}}, k::Int)
    reader = FASTA.Reader(open(fasta_path, "r"))
    record = first(reader)
    close(reader)

    seq_str = FASTA.sequence(String, record)
    cds_index = Dict{Tuple{Int, Int}, Set{DNAKmer{k}}}()

    for (start_idx, end_idx) in cds_coords
        start_idx = max(1, start_idx)
        end_idx = min(length(seq_str), end_idx)

        gene_kmers = Set{DNAKmer{k}}()
        if start_idx <= end_idx
            sub_seq_str = seq_str[start_idx:end_idx]
            for chunk in split(sub_seq_str, r"[^ACGTacgt]")
                if length(chunk) >= k
                    chunk_seq = LongDNA{4}(chunk)
                    union!(gene_kmers, CanonicalDNAMers{k}(chunk_seq))
                end
            end
        end
        cds_index[(start_idx, end_idx)] = gene_kmers
    end

    return cds_index
end

"""
    profile_cds(sink_fastq::String, cds_index::Dict, k::Int) -> Dict{Tuple{Int, Int}, Float64}

Compute fractional `k`-mer coverage across reads in `sink_fastq` for each coding sequence in `cds_index`.
"""
function profile_cds(sink_fastq::String, cds_index::Dict, k::Int)
    sink_kmers = Set{DNAKmer{k}}()
    reader = FASTQ.Reader(open(sink_fastq, "r"))

    try
        for record in reader
            kmers = extract_kmers(record, k)
            union!(sink_kmers, kmers)
        end
    finally
        close(reader)
    end

    results = Dict{Tuple{Int, Int}, Float64}()

    for (gene_coords, gene_kmers) in cds_index
        if isempty(gene_kmers)
            results[gene_coords] = 0.0
        else
            coverage = length(intersect(gene_kmers, sink_kmers)) / length(gene_kmers)
            results[gene_coords] = coverage
        end
    end

    return results
end

"""
    classify_cds_profile(cds_profile::Dict{Tuple{Int, Int}, Float64}; detection_threshold::Float64=0.05, intact_threshold::Float64=0.95)

Partition CDS coordinates from `cds_profile` into `(intact, fragmented, undetected)` vectors based on `detection_threshold` and `intact_threshold`.
"""
function classify_cds_profile(
    cds_profile::Dict{Tuple{Int, Int}, Float64};
    detection_threshold::Float64=0.05,
    intact_threshold::Float64=0.95
)
    intact = Tuple{Int, Int}[]
    fragmented = Tuple{Int, Int}[]
    undetected = Tuple{Int, Int}[]

    for (coords, cov) in cds_profile
        if cov >= intact_threshold
            push!(intact, coords)
        elseif cov > detection_threshold
            push!(fragmented, coords)
        else
            push!(undetected, coords)
        end
    end

    return (intact=intact, fragmented=fragmented, undetected=undetected)
end
