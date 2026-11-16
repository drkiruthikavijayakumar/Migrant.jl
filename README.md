# Migrant.jl

[![CI](https://github.com/drkiruthikavijayakumar/Migrant.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/drkiruthikavijayakumar/Migrant.jl/actions/workflows/CI.yml)

## Overview
`Migrant.jl` is an alignment-free, memory-efficient Julia package for high-resolution microbial source tracking via whole-non-coding DNA (ncDNA) fingerprinting and functional coding sequence (CDS) profiling. By bypassing traditional host-depletion memory bottlenecks via positive target extraction, it isolates hyper-variable intergenic *k*-mers to verify clonal translocation while simultaneously profiling functional CDS retention and structural integrity.

## Installation
You can install `Migrant.jl` using the Julia package manager:

```julia
] add Migrant
```

## Quickstart
Here is a basic example of tracking translocating bacterial DNA from a source (e.g., saliva) to a sink (e.g., stool or tumor):

```julia
using Migrant

# 1. Build the target index and extract bacterial reads
target_idx = build_target_index("data/bacterial_ref.fasta", 31)
extract_target_reads("data/oral.fastq", "clean_oral.fastq", target_idx, 31, 0.05)
extract_target_reads("data/tumor.fastq", "clean_tumor.fastq", target_idx, 31, 0.05)

# 2. Analyze translocation via the core engine
res = analyze_migration(
    "clean_oral.fastq",
    "clean_tumor.fastq",
    "data/bacterial_ref.fasta",
    "data/bacterial_annotations.gff3";
    k=31,
    target_threshold=0.05
)

println("Strain Translocation Score (C): ", round(res.containment_score * 100, digits=2), "%")
```

## Citation
If you use `Migrant.jl` in your research, please cite the accompanying manuscript and software release (`CITATION.cff`).

## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
