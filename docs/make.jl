using JEG
using Documenter

makedocs(
    sitename="Julia Electric Grid",
    modules = [JEG],
    pages = [
        "Welcome" => "index.md",
        "Guide" => [
            "Quick Start" => "quickstart.md",
            "The Nodeconstructor - Theory" => "NodeConstructor_Theory.md",
            "The Nodeconstructor - Application" => "NodeConstructor_Application.md",
            "Set up Classical Controllers" => "classical.md",
            "Contributing" => "dev.md"
        ],
        "Environment" => [
            "Configuring the Environment" => "Env_Create.md",
            "Interaction with the Environment" => "Env_Interaction.md"
        ],
        "Reinforcement Learning" => [
            "Reinforcement Learning using JEG" => "RL_Single_Agent.md",
            "Multicontroller" => "RL_Classical_Controllers_Merge.md",
            "Reinforcement Learning in Larger Grids " => "RL_Complex.md",
        ],
        "References" => "references.md",
    ],
    format = Documenter.HTML(
        mathengine = MathJax3(Dict(
            :loader => Dict("load" => ["[tex]/physics"]),
            :tex => Dict(
                "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
                "tags" => "ams",
                "packages" => ["base", "ams", "autoload", "physics"],
            ),
        )),
        assets = [
            #asset("https://cdn.plot.ly/plotly-2.3.0.min.js")
        ],
    ),
    )

deploydocs(
        repo = "github.com/upb-lea/JuliaElectricGrid.jl.git",
    )
