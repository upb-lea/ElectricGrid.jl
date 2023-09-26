using ElectricGrid
using Documenter

makedocs(
    sitename="ElectricGrid.jl",
    modules = [ElectricGrid],
    pages = [
        "Welcome" => "index.md",
        "Quick Start" => "quickstart.md",
        "Environment" => [
            "Configuring the Environment" => "Env_Create.md",
            "Interaction with the Environment" => "Env_Interaction.md"
        ],
        "Classical Controllers" => [
            "Swing Mode" => "Classical_Controllers_Swing.md",
            "PQ Mode" => "Classical_Controllers_PQ.md",
            "Droop Controllers" => "Classical_Controllers_Droop.md",
            "Virtual Synchronous Generator" => "Classical_Controllers_VSG.md",
            "Auxiliaries" => "Auxiliaries_OU_process.md",
        ],
        "Reinforcement Learning" => [
            "Reinforcement Learning using ElectricGrid" => "RL_Single_Agent.md",
            "Multicontroller" => "RL_Classical_Controllers_Merge.md",
            "Reinforcement Learning in Larger Grids " => "RL_Complex.md",
        ],
        "Nodeconstructor" => [
            "The Nodeconstructor - Theory" => "NodeConstructor_Theory.md",
            "The Nodeconstructor - Application" => "NodeConstructor_Application.md"
        ],
        "Miscellanous" => [
            "Default Parameters" => "Default_Parameters.md"
            "Nonlinear Components" => "Nonlinear.md"
            "DC Link Models" => "DC_link_models.md"
            "GUI" => "Gui.md"
        ],
        "References" => "references.md",
        "Contributing" => "dev.md",
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
        repo = "github.com/upb-lea/ElectricGrid.jl.git",
    )
