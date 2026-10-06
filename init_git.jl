using Pkg
using LibGit2

repo_path = "."
repo = LibGit2.init(repo_path)
println("Initialized empty Git repository in ", repo_path)
