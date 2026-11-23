using FASTX
using BioSequences
using Kmers

"""
    extract_kmers(record::FASTQ.Record, k::Int) -> Set{DNAKmer{k}}

Extract all canonical `k`-mers from a `FASTQ.Record`, splitting across any non-ACGT ambiguous bases.
"""
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

"""
    build_target_index(fasta_path::String, k::Int) -> Set{DNAKmer{k}}

Build a canonical `k`-mer index from all sequences in a target reference FASTA file at `fasta_path`.
"""
function build_target_index(fasta_path::String, k::Int)
    reader = FASTA.Reader(open(fasta_path, "r"))
    target_index = Set{DNAKmer{k}}()
    try
        for record in reader
            seq_str = FASTA.sequence(String, record)
            for chunk in split(seq_str, r"[^ACGTacgt]")
                if length(chunk) >= k
                    chunk_seq = LongDNA{4}(chunk)
                    union!(target_index, CanonicalDNAMers{k}(chunk_seq))
                end
            end
        end
    finally
        close(reader)
    end
    return target_index
end

"""
    extract_target_reads(fastq_in::String, fastq_out::String, target_index::Set, k::Int, threshold::Float64)

Filter reads from `fastq_in` whose fraction of `k`-mers matching `target_index` is at least `threshold`, writing matching records to `fastq_out`.
"""
function extract_target_reads(fastq_in::String, fastq_out::String, target_index::Set, k::Int, threshold::Float64)
    reader = FASTQ.Reader(open(fastq_in, "r"))
    writer = FASTQ.Writer(open(fastq_out, "w"))
    try
        for record in reader
            read_kmers = extract_kmers(record, k)

            if isempty(read_kmers)
                continue
            end

            overlap_frac = length(intersect(read_kmers, target_index)) / length(read_kmers)

            if overlap_frac >= threshold
                write(writer, record)
            end
        end
    finally
        close(reader)
        close(writer)
    end
end
