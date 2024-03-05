using LevenbergMarquardtFit
using Documenter

DocMeta.setdocmeta!(LevenbergMarquardtFit, :DocTestSetup, :(using LevenbergMarquardtFit); recursive=true)

makedocs(;
    modules=[LevenbergMarquardtFit],
    authors="Alberto Mercurio",
    sitename="LevenbergMarquardtFit.jl",
    format=Documenter.HTML(;
        canonical="https://albertomercurio.github.io/LevenbergMarquardtFit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/albertomercurio/LevenbergMarquardtFit.jl",
    devbranch="main",
)
