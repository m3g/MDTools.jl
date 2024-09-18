import Pkg
Pkg.add("Documenter")
using Documenter
push!(LOAD_PATH, "../")
push!(LOAD_PATH, "../src/")
using MolSimToolkit
makedocs(
    modules=[MolSimToolkit],
    sitename="MolSimToolkit.jl",
    doctest=false,
    pages=[
        "Home" => "index.md",
        "System setup" => "system_setup.md",
        "Molecular Minimum Distances" => "molecular_minimum_distances.md",
        "Block averages" => "block_averages.md",
        "Replica exchange" => "remd.md",
        "Structural alignment" => "procrustes.md",
        "Secondary structure" => "secondary_structures.md",
        "Miscelaneous functions" => "miscelaneous.md",
        "Developer zone" => "Developer.md",
        "Plotting style" => "plotting_style.md",
        "Experimental" => "Experimental.md",
        " ∘ Simulation Reweighting" => "Reweighting.md"
    ],
)
deploydocs(
    repo="github.com/m3g/MolSimToolkit.jl.git",
    target="build",
    branch="gh-pages",
    versions=["stable" => "v^", "v#.#"],
)
