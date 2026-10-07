using FASTX
using BioSequences
using Kmers

function build_ncdna_footprint(fasta_path::String, coords::Array{Tuple{Int, Int}}, k::Int)
    reader = FASTA.Reader(open(fasta_path, "r"))
    # Assume 1 main chromosome/record
    record = first(reader)
    close(reader)
    
    seq_str = FASTA.sequence(String, record)
    source_kmers = Set{DNAKmer{k}}()
    
    for (start_idx, end_idx) in coords
        # Ensure indices are within bounds
        start_idx = max(1, start_idx)
        end_idx = min(length(seq_str), end_idx)
        
        if start_idx > end_idx
            continue
        end
        
        sub_seq_str = seq_str[start_idx:end_idx]
        for chunk in split(sub_seq_str, r"[^ACGTacgt]")
            if length(chunk) >= k
                chunk_seq = LongDNA{4}(chunk)
                union!(source_kmers, CanonicalDNAMers{k}(chunk_seq))
            end
        end
    end
    
    return source_kmers
end

function calculate_kmer_containment(sink_fastq::String, source_kmers::Set, k::Int)
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
    
    if length(sink_kmers) == 0
        return 0.0
    end
    
    containment = length(intersect(sink_kmers, source_kmers)) / length(sink_kmers)
    return containment
end
