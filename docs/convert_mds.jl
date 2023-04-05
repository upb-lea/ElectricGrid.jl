using Glob

function process_file(filepath::AbstractString)
    edit_latex_block = false
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


            #remove julia logger artifacts
            regex_pattern = r"[^]*?m"
            line = replace(line, regex_pattern => "")

            #change figures/ to assets/
            line = replace(line, "figures/" => "assets/")

            push!(output_lines, line)

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

dirpath = string(@__DIR__) * "./src"
md_files = glob("*.md", dirpath)

for file in md_files
    process_file(file)
end