using Documenter, SimpleCrystals

makedocs(sitename = "SimpleCrystals.jl",
format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
), 
pages = [
    "Home"                      => "index.md"
]
)

deploydocs(
    repo="github.com/ejmeitz/SimpleCrystals.jl.git",
    push_preview = true,
    devbranch = "main"
)