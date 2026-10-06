using FASTX
using BioSequences
using Kmers

function extract_kmers(record::FASTQ.Record, k::Int)
    seq_str = FASTQ.sequence(String, record)
    kmers = Set{DNAKmer{k}}()
    for chunk in split(seq_str, r"[^ACGTacgt]")
        if length(chunk) >= k
            chunk_seq = LongDNA{4}(chunk)
            union!(kmers, CanonicalDNAMers{k}(chunk_seq))
        end
    end
    return kmers
end

function build_host_index(fasta_path::String, k::Int)
    reader = FASTA.Reader(open(fasta_path, "r"))
    host_index = Set{DNAKmer{k}}()
    try
        for record in reader
            seq_str = FASTA.sequence(String, record)
            for chunk in split(seq_str, r"[^ACGTacgt]")
                if length(chunk) >= k
                    chunk_seq = LongDNA{4}(chunk)
                    union!(host_index, CanonicalDNAMers{k}(chunk_seq))
                end
            end
        end
    finally
        close(reader)
    end
    return host_index
end

function filter_host_reads(fastq_in::String, fastq_out::String, host_index::Set, k::Int, threshold::Float64)
    reader = FASTQ.Reader(open(fastq_in, "r"))
    writer = FASTQ.Writer(open(fastq_out, "w"))
    try
        for record in reader
            read_kmers = extract_kmers(record, k)
            
            if isempty(read_kmers)
                write(writer, record)
                continue
            end
            
            overlap_frac = length(intersect(read_kmers, host_index)) / length(read_kmers)
            
            if overlap_frac < threshold
                write(writer, record)
            end
        end
    finally
        close(reader)
        close(writer)
    end
end
