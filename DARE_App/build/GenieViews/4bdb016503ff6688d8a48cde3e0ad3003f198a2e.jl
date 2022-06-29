# app\resources\ncconfigs\views\_ncconfig.jl.html 

function func_4bdb016503ff6688d8a48cde3e0ad3003f198a2e(;
    context = Genie.Renderer.vars(:context),
    ncconfig = Genie.Renderer.vars(:ncconfig),
    ncconfigs = Genie.Renderer.vars(:ncconfigs),
)

    [
        Genie.Renderer.Html.div(
            class = "container",
            style = "margin-top: 40px;",
            htmlsourceindent = "2",
        ) do
            [
                Genie.Renderer.Html.h3(htmlsourceindent = "3") do
                    [
                        ncconfig.name
                    ]
                end
                Genie.Renderer.Html.div(htmlsourceindent = "3") do
                    [
                        Genie.Renderer.Html.strong(htmlsourceindent = "4") do
                            [
                                """CM: """;
                            ]
                        end
                        ncconfig.cm
                    ]
                end
                Genie.Renderer.Html.div(htmlsourceindent = "3") do
                    [
                        Genie.Renderer.Html.strong(htmlsourceindent = "4") do
                            [
                                """Parameters: """;
                            ]
                        end
                        ncconfig.parameters
                    ]
                end
            ]
        end,
    ]
end
