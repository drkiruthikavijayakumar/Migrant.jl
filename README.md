# migrant.jl

[![CI](https://github.com/username/migrant.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/username/migrant.jl/actions/workflows/CI.yml)

## Overview
**migrant.jl** is a high-performance Julia package for tracing tumor microbiome translocation via whole-non-coding DNA profiling and virulence tracking. It provides specialized tools for extracting host and non-coding sequence footprints to efficiently isolate translocating bacterial signatures from sequencing data.

## Installation
You can install `migrant.jl` using the Julia package manager:
```julia
] add migrant
```

## Quickstart
Here is a basic example of tracking translocating DNA:

```julia
using migrant

# 1. Build the host index and deplete host sequences
host_idx = build_host_index("data/host_ref.fasta", 31)
filter_host_reads("data/oral.fastq", "clean_oral.fastq", host_idx, 31, 0.1)
filter_host_reads("data/tumor.fastq", "clean_tumor.fastq", host_idx, 31, 0.1)

# 2. Analyze the migration via the core engine
res = analyze_migration(
    "clean_oral.fastq",
    "clean_tumor.fastq",
    "data/host_ref.fasta",
    "data/bacterial_ref.fasta",
    "data/bacterial_annotations.gff3";
    k = 31,
    host_threshold = 0.1
)

println("Migration Containment Score: ", round(res.containment_score * 100, digits=2), "%")
```
