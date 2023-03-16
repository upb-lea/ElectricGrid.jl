using JEG
using Documenter

makedocs(
    sitename="Julia Electric Grid",
    modules = [JEG],
    pages = [
        "Welcome" => "index.md",
        "Guide" => [
            "Quick Start" => "quickstart.md",
            "Configuring the Environment" => "quickstart.md",
            "Set up Classical Controllers" => "quickstart.md",
            "Set up Agents" => "quickstart.md",
            "The Nodeconstructor" => "quickstart.md",
        ],
        "References" => "references.md",
    ]
    )