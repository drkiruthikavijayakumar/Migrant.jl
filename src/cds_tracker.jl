using FASTX
using BioSequences
using Kmers

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
