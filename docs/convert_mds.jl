using Glob

function process_file(filepath::AbstractString)
    edit_block = false
    output_lines = []
    prev_line = ""

    open(filepath) do file
        for line in eachline(file)

            #check for beginning of unencapsulated Latex blocks
            if occursin(r"\\begin", line)
                if !occursin("```math", prev_line)
                    edit_block = true
                    push!(output_lines, "```math")
                end
            end


            #remove julia logger artifacts
            regex_pattern = r"[^]*?m"
            line = replace(line, regex_pattern => "")

            push!(output_lines, line)

            #check for end of unencapsulated Latex blocks
            if occursin(r"\\end", line)
                if edit_block
                    push!(output_lines, "```")
                    edit_block = false
                end
            end

            prev_line = line
        end
    end

    open(filepath, "w") do file
        for line in output_lines
            println(file, line)
        end
    end
end

dirpath = string(@__DIR__) * "./src"
md_files = glob("*.md", dirpath)

for file in md_files
    process_file(file)
end