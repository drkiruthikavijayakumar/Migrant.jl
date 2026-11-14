# migrant.jl

[![CI](https://github.com/username/migrant.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/username/migrant.jl/actions/workflows/CI.yml)

## Overview
**migrant.jl** is a high-performance Julia package for tracing tumor microbiome translocation via whole-non-coding DNA profiling and functional coding sequence (CDS) tracking. It provides specialized tools for extracting target bacterial sequences and non-coding footprints to efficiently isolate translocating bacterial signatures from sequencing data.

## Installation
You can install `migrant.jl` using the Julia package manager:
```julia
] add migrant
```

## Quickstart
Here is a basic example of tracking translocating DNA:

```julia
using Migrant

# 1. Build the target index and extract bacterial sequences
target_idx = build_target_index("data/bacterial_ref.fasta", 31)
extract_target_reads("data/oral.fastq", "clean_oral.fastq", target_idx, 31, 0.05)
extract_target_reads("data/tumor.fastq", "clean_tumor.fastq", target_idx, 31, 0.05)

# 2. Analyze the migration via the core engine
res = analyze_migration(
    "clean_oral.fastq",
    "clean_tumor.fastq",
    "data/bacterial_ref.fasta",
    "data/bacterial_annotations.gff3";
    k = 31,
    target_threshold = 0.05
)

println("Migration Containment Score: ", round(res.containment_score * 100, digits=2), "%")
```
