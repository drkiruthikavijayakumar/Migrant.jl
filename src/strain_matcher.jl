using FASTX
using BioSequences
using Kmers

function build_ncdna_footprint(fasta_path::String, coords::Array{Tuple{Int, Int}}, k::Int)
    reader = FASTA.Reader(open(fasta_path, "r"))
    # Assume 1 main chromosome/record
    record = first(reader)
    close(reader)
    
    seq = FASTA.sequence(LongDNA{4}, record)
    oral_kmers = Set{DNAKmer{k}}()
    
    for (start_idx, end_idx) in coords
        # Ensure indices are within bounds
        start_idx = max(1, start_idx)
        end_idx = min(length(seq), end_idx)
        
        if start_idx > end_idx
            continue
        end
        
        sub_seq = seq[start_idx:end_idx]
        for kmer in CanonicalDNAMers{k}(sub_seq)
            push!(oral_kmers, kmer)
        end
    end
    
    return oral_kmers
end

function calculate_kmer_containment(tumor_fastq::String, oral_kmers::Set, k::Int)
    tumor_kmers = Set{DNAKmer{k}}()
    reader = FASTQ.Reader(open(tumor_fastq, "r"))
    
    try
        for record in reader
            kmers = extract_kmers(record, k)
            union!(tumor_kmers, kmers)
        end
    finally
        close(reader)
    end
    
    if length(tumor_kmers) == 0
        return 0.0
    end
    
    containment = length(intersect(tumor_kmers, oral_kmers)) / length(tumor_kmers)
    return containment
end
