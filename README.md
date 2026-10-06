# migrant.jl

[![CI](https://github.com/username/migrant.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/username/migrant.jl/actions/workflows/CI.yml)

## Overview
**migrant.jl** is an end-to-end clinical oncology pipeline written in Julia. It provides tools for tracking tumor microbiome migration using deep oral non-coding DNA (ncDNA) footprints, alongside comprehensive Virulence Profiling. The package uses a built-in Machine Learning engine to predict oncogenic (tumor tropism) risk directly from sequencing data.

## Installation
You can install `migrant.jl` using the Julia package manager:
```julia
] add migrant
```

## Quickstart
Here is a basic example of running the end-to-end Clinical Decision Support Tool (CDST) pipeline:

```julia
using migrant
using MLJ

# Load your trained MLJ machine (tropism model)
# mach = ... 

risk_score = run_cdst_pipeline(
    "Patient_001",
    "data/oral.fastq",
    "data/tumor.fastq",
    "data/host_ref.fasta",
    "data/bacterial_ref.fasta",
    "data/bacterial_annotations.gff3",
    mach,
    "patient_001_report.html";
    k = 31,
    host_threshold = 0.1
)

println("Tumor Tropism Risk: ", round(risk_score * 100, digits=2), "%")
```
