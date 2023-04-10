using JEG
using Documenter

makedocs(
    sitename="Julia Electric Grid",
    modules = [JEG],
    pages = [
        "Welcome" => "index.md",
        "Guide" => [
            "Quick Start" => "quickstart.md",
            "Configuring the Environment" => "environment.md",
            "Set up Classical Controllers" => "classical.md",
            "Set up Agents" => "agents.md",
            "The Nodeconstructor" => "nodeconstructor.md",
        ],
        "References" => "references.md",
    ]
    )