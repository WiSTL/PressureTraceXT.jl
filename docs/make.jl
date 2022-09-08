using PressureTraceXT
using Documenter

DocMeta.setdocmeta!(PressureTraceXT, :DocTestSetup, :(using PressureTraceXT); recursive=true)

makedocs(;
    modules=[PressureTraceXT],
    authors="stillyslalom <aames@wisc.edu> and contributors",
    repo="https://github.com/stillyslalom/PressureTraceXT.jl/blob/{commit}{path}#{line}",
    sitename="PressureTraceXT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://stillyslalom.github.io/PressureTraceXT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/stillyslalom/PressureTraceXT.jl",
    devbranch="main",
)
