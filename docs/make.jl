using JEG
using Documenter

makedocs(
    sitename="Julia Electric Grid",
    modules = [JEG],
    pages = [
        "Welcome" => "index.md",
        "Guide" => [
            "Quick Start" => "quickstart.md",
            "Configuring the Environment" => "Env_Create.md",
            "Interaction with the Environment" => "Env_Interaction.md",
            "The Nodeconstructor - Theory" => "NodeConstructor_Theory.md",
            "The Nodeconstructor - Application" => "NodeConstructor_Application.md",
            "Set up Classical Controllers" => "classical.md",
            "Train an RL againt in the JEG framework - 1" => "RL_Single_Agent.md",
            "Train an RL againt interacting with a stable grid" => "RL_Classical_Controllers_Merge.md",
            "Train an RL againt in the JEG framework - 2" => "RL_Complex.md",
        ],
        "References" => "references.md",
    ]
    )

deploydocs(
        repo = "github.com/upb-lea/JuliaElectricGrid.jl.git",
    )
