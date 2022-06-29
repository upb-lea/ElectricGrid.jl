# app\resources\ncconfigs\views\index.jl.html 

function func_3287ab05ec4c35bc5b0dd6a45b70e9d533018484(;
    context = Genie.Renderer.vars(:context),
    ncconfigs = Genie.Renderer.vars(:ncconfigs),
)

    [
        Genie.Renderer.Html.h1(class = "display-1 text-center", htmlsourceindent = "2") do
            [
                """Watch tonight""";
            ]
        end
        Genie.Renderer.Html.div(
            class = "container",
            style = "margin-top: 40px;",
            htmlsourceindent = "2",
        ) do
            [
                Genie.Renderer.Html.form(
                    action = "$( linkto(:search_ncconfigs) )",
                    htmlsourceindent = "3",
                ) do
                    [
                        Genie.Renderer.Html.input(
                            name = "search_ncconfigs",
                            class = "form-control form-control-lg",
                            htmlsourceindent = "4",
                            placeholder = "Search for NC configs",
                            type = "search",
                        )
                    ]
                end;
            ]
        end
        if !isempty(ncconfigs)
            for_each(ncconfigs) do ncconfig
                partial(
                    joinpath(
                        Genie.config.path_resources,
                        "ncconfigs",
                        "views",
                        "_ncconfig.jl.html",
                    ),
                    ncconfig = ncconfig,
                )
            end
        else
            partial(
                joinpath(
                    Genie.config.path_resources,
                    "ncconfigs",
                    "views",
                    "_no_results.jl.html",
                ),
            )
        end
    ]
end
