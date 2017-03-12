using Documenter
using LibXC
makedocs(
         modules = [LibXC],
         clean = true,
         format = :html,
         sitename = "LibXC.jl",
         authors = "Mayeul d'Avezac",
         analytics = "UA-89508993-1",
         pages = Any[
                     "Home" => "index.md",
                    ]
        )

deploydocs(
    repo = "github.com/mdavezac/LibXC.jl.git",
    target = "build",
    julia = "release",
    deps = nothing,
    make = nothing,
)
