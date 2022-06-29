# app\resources\ncconfigs\views\_no_results.jl.html 

function func_d4d6f7b98596f0a5f7adcf32dfcd7d3d639a7033(;
    context = Genie.Renderer.vars(:context),
    ncconfigs = Genie.Renderer.vars(:ncconfigs),
)

    [
        Genie.Renderer.Html.h4(class = "container", htmlsourceindent = "2") do
            [
                """
                    Sorry, no results were found for "$(params(:search_ncconfigs))"
                  """
            ]
        end,
    ]
end
