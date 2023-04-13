using Glob
using Random

function process_file(filepath::AbstractString)
    edit_latex_block = false
    edit_html_block = false
    edit_plotlyjs_block = false
    latex_dollar = false
    plotly_div_name = ""
    html_divs = 0
    original_lines = []
    output_lines = []
    prev_line = ""

    open(filepath) do file
        for line in eachline(file)
            push!(original_lines, line)

            #check for beginning of unencapsulated Latex blocks
            if occursin(r"\\begin{alig", line) || occursin(r"\\begin{equa", line)
                if !occursin("```math", prev_line)
                    edit_latex_block = true
                    push!(output_lines, "```math")
                end
            end

            if occursin("\$\$", line)
                if !latex_dollar
                    latex_dollar = true
                    if !occursin("```math", prev_line)
                        edit_latex_block = true
                        line = ""
                        push!(output_lines, "```math")
                    end
                else
                    latex_dollar = false
                    if edit_latex_block
                        line = ""
                        push!(output_lines, "```")
                        edit_latex_block = false
                    end
                end
            end

            #check for beginning of unencapsulated html blocks
            if occursin(r"<html", line)
                if !occursin("```@raw html", prev_line)
                    edit_html_block = true
                    push!(output_lines, "```@raw html")
                end
                html_divs = 1 + length(collect(eachmatch(r"<div", line)))
            end

            #check for beginning of unencapsulated div blocks
            if occursin(r"<div", line)
                if !occursin("```@raw html", prev_line) && html_divs == 0
                    edit_html_block = true
                    push!(output_lines, "```@raw html")
                end
                html_divs += length(collect(eachmatch(r"<div", line)))
            end

            #check for PlotlyJS plot beginning
            if occursin(r"webio-mountpoint", line) && occursin(r"<div", prev_line)
                edit_plotlyjs_block = true
                plotly_div_name = randstring(12)
                push!(output_lines, "id = $plotly_div_name > </div>")
                push!(output_lines, "<script>")
                push!(output_lines, "gd = '$plotly_div_name'")
                line = "require(['plotly'], function(plotly) {"
            elseif edit_plotlyjs_block
                html_divs -= length(collect(eachmatch(r"</div", line)))
                if occursin(r"Plotly.newPlot", line)
                    pline = "plotly." * match(r"newPlot.*?\)\;", line).match
                    push!(output_lines, replace(pline, r"\\\\\"" => "\""))
                    push!(output_lines, "});")
                    line = ""
                elseif html_divs > 0
                    line = ""
                else
                    push!(output_lines, "</script>")
                    if edit_html_block
                        push!(output_lines, "```")
                        edit_html_block = false
                    end
                    line = ""
                    edit_plotlyjs_block = false
                end
            end



            #remove julia logger artifacts
            regex_pattern = r"[^]*?m"
            line = replace(line, regex_pattern => "")

            #change figures/ to assets/
            line = replace(line, "figures/" => "assets/")

            #change  "") to )
            line = replace(line, " \"\")" => ")")


            push!(output_lines, line)


            #check for end of unencapsulated div blocks
            if occursin(r"</div", line)
                if html_divs > 0
                    html_divs -= length(collect(eachmatch(r"</div", line)))
                    if html_divs <= 0 && edit_html_block
                        push!(output_lines, "```")
                        edit_html_block = false
                    end
                end
            end

            #check for end of unencapsulated html blocks
            if occursin(r"</html", line)
                if html_divs > 0
                    html_divs -= 1 + length(collect(eachmatch(r"</div", line)))
                    if html_divs <= 0 && edit_html_block
                        push!(output_lines, "```")
                        edit_html_block = false
                    end
                end
            end

            #check for end of unencapsulated Latex blocks
            if occursin(r"\\end{alig", line) || occursin(r"\\end{equa", line)
                if edit_latex_block
                    push!(output_lines, "```")
                    edit_latex_block = false
                end
            end

            prev_line = line
        end
    end

    if output_lines != original_lines
        open(filepath, "w") do file
            for line in output_lines
                println(file, line)
            end
        end
    end
end

dirpath = string(@__DIR__) * "/src"
md_files = glob("*.md", dirpath)

for file in md_files
    process_file(file)
end
