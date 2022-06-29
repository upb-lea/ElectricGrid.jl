# app\layouts\app.jl.html 

function func_4e4692c73ffd62174d1e1ae92682bcca8d9715d1(;
    context = Genie.Renderer.vars(:context),
    ncconfigs = Genie.Renderer.vars(:ncconfigs),
)

    [
        Genie.Renderer.Html.doctype()
        Genie.Renderer.Html.html(htmlsourceindent = "0", lang = "en") do
            [
                Genie.Renderer.Html.head(htmlsourceindent = "1") do
                    [
                        Genie.Renderer.Html.meta(charset = "utf-8", htmlsourceindent = "2")
                        Genie.Renderer.Html.title(htmlsourceindent = "2") do
                            [
                                """Genie :: The Highly Productive Julia Web Framework""";
                            ]
                        end
                        Genie.Renderer.Html.link(
                            crossorigin = "anonymous",
                            htmlsourceindent = "2",
                            href = "https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta2/dist/css/bootstrap.min.css",
                            rel = "stylesheet",
                            integrity = "sha384-BmbxuPwQa2lc/FVzBcNJ7UAyJxM6wuqIj61tLrc4wSX0szH/Ev+nYRRuWlolflfl",
                        )
                    ]
                end
                Genie.Renderer.Html.body(htmlsourceindent = "1") do
                    [
                        @yield
                    ]
                end
            ]
        end
    ]
end
