
function create_folder_name(file_name::String)
    file_position = findfirst('.', file_name)
    last_folder_position = findlast('/', file_name)

    folder_name = file_name[begin:last_folder_position - 1]
    split_file_name = "split_" * file_name[last_folder_position + 1 : file_position - 1] * '_'
    return folder_name, split_file_name
end


function line_splitter(infile::IOStream, file_name::String, lines_per_file::Int)
    line_count::Int64 = 0
    file_count::Int64 = 0
    folder_name, split_file_name = create_folder_name(file_name)
    outfile::IOStream = open("$(folder_name)/head_$(split_file_name)$(file_count)", "w")

    for line in eachline(infile)
        write(outfile, line * '\n')
        line_count += 1
        if iszero(line_count % lines_per_file)
            close(outfile)
            file_count += 1
            outfile = open("$(folder_name)/head_$(split_file_name)$(file_count)", "w")
        end
    end
end



function main()
    if length(ARGS) != 2
        println("No arguments were provided\n\n" *
                "julia file_splitter.jl [LINE_COUNT] [FILE]\n\n" * 
                "<SPLIT BASED ON LINES>\n" * 
                "---------\n" *
                "LINECOUNT - the amount of lines to split in each file" * 
                "<FILE>\n" * 
                "Please provide a file with its absolute path as follows: /absolute/path/to/file/file.txt\n")
    elseif length(ARGS) == 2
        lines_per_file::Int = parse(Int, ARGS[1])
        infile::IOStream = open(ARGS[2])
        line_splitter(infile, ARGS[2], lines_per_file)
    end
end
main()