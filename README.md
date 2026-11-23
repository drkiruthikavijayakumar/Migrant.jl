# Migrant.jl

[![CI](https://github.com/drkiruthikavijayakumar/Migrant.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/drkiruthikavijayakumar/Migrant.jl/actions/workflows/CI.yml)

## Overview
**Migrant.jl** is an alignment-free, memory-efficient Julia package for tracing microbial translocation via whole non-coding DNA (ncDNA) fingerprinting and functional coding sequence (CDS) profiling. It provides specialized tools for extracting target bacterial reads, building non-coding footprints, and quantifying gene-level coverage from metagenomic sequencing data.

## Installation
Install `Migrant.jl` using the Julia package manager:
```julia
using Pkg
Pkg.add("Migrant")
```

## Quickstart
Here is a self-contained example demonstrating target read extraction and translocation scoring:

```julia
using Migrant

# Create a small synthetic reference genome, GFF3 annotation, and FASTQ reads
workdir = mktempdir()
ref_fasta = joinpath(workdir, "ref.fasta")
ref_gff   = joinpath(workdir, "ref.gff3")
oral_fq   = joinpath(workdir, "oral.fastq")
tumor_fq  = joinpath(workdir, "tumor.fastq")

write(ref_fasta, ">chr1\n" * "ACGT"^15 * "TGCA"^15 * "\n")
write(ref_gff, "##gff-version 3\nchr1\tsrc\tCDS\t61\t120\t.\t+\t0\tID=cds1\n")
write(oral_fq,  "@r1\n" * "ACGT"^15 * "\n+\n" * "I"^60 * "\n")
write(tumor_fq, "@r2\n" * "ACGT"^15 * "\n+\n" * "I"^60 * "\n")

# 1. Build target k-mer index and filter bacterial reads
target_idx = build_target_index(ref_fasta, 15)
extract_target_reads(oral_fq, joinpath(workdir, "clean_oral.fastq"), target_idx, 15, 0.05)
extract_target_reads(tumor_fq, joinpath(workdir, "clean_tumor.fastq"), target_idx, 15, 0.05)

# 2. Run end-to-end migration and CDS functional profiling
res = analyze_migration(oral_fq, tumor_fq, ref_fasta, ref_gff; k = 15, target_threshold = 0.05)
println("Migration Containment Score: ", round(res.containment_score * 100, digits=2), "%")
```

Sample *Fusobacterium nucleatum* reference files (`examples/f_nucleatum.fna` and `examples/f_nucleatum.gff`) are also included in the `examples/` directory for testing on full bacterial genomes.
