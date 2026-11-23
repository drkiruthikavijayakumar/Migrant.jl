using Documenter
using Migrant

makedocs(
    sitename = "Migrant.jl",
    modules = [Migrant],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ]
)
