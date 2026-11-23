"""
    parse_cds_coordinates(gff_path::String) -> Vector{Tuple{Int, Int}}

Parse 1-based `(start, end)` genomic coordinates for all `CDS` features in a GFF3 file, sorted by start position.
"""
function parse_cds_coordinates(gff_path::String)
    cds_coords = Tuple{Int, Int}[]
    for line in eachline(gff_path)
        line = strip(line)
        if isempty(line) || startswith(line, "#")
            continue
        end

        parts = split(line, '\t')
        if length(parts) >= 5 && parts[3] == "CDS"
            start_idx = parse(Int, parts[4])
            end_idx = parse(Int, parts[5])
            push!(cds_coords, (start_idx, end_idx))
        end
    end
    sort!(cds_coords, by = x -> x[1])
    return cds_coords
end

"""
    extract_ncdna_regions(cds_coords::Array{Tuple{Int, Int}}, genome_length::Int) -> Vector{Tuple{Int, Int}}

Compute the complementary non-coding DNA (ncDNA) intervals across a genome of length `genome_length` given sorted `cds_coords`.
"""
function extract_ncdna_regions(cds_coords::Array{Tuple{Int, Int}}, genome_length::Int)
    ncdna = Tuple{Int, Int}[]
    current_pos = 1

    for (start_idx, end_idx) in cds_coords
        if current_pos < start_idx
            push!(ncdna, (current_pos, start_idx - 1))
        end
        current_pos = max(current_pos, end_idx + 1)
    end

    if current_pos <= genome_length
        push!(ncdna, (current_pos, genome_length))
    end

    return ncdna
end
