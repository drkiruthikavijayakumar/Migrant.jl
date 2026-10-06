using FASTX
using BioSequences
using Kmers

function extract_kmers(record::FASTQ.Record, k::Int)
    seq = FASTQ.sequence(LongDNA{4}, record)
    return Set{DNAKmer{k}}(CanonicalDNAMers{k}(seq))
end

function build_host_index(fasta_path::String, k::Int)
    reader = FASTA.Reader(open(fasta_path, "r"))
    host_index = Set{DNAKmer{k}}()
    try
        for record in reader
            seq = FASTA.sequence(LongDNA{4}, record)
            for kmer in CanonicalDNAMers{k}(seq)
                push!(host_index, kmer)
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
