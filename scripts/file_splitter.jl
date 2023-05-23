#=
FILE SPLITTER
------------
This script splits up a mother file into multiple children files
When RAM usage is concerned in certain programs.
And that RAM usage can be split when the files themselves are split.
Though it should be stated that having anything less than 8GB is 
    not recommended to use a pipeline on since it is probably not enough
    RAM / memory to do the analyses to begin with.
    Since numerous assembly programs use terribly high numbers of RAM
    without even having mapped any of the reads.

Version: 1.0
Date: 2023-05-23
Author: Marc Wijnands
=#


function create_folder_name(file_name::String)
    #=
    Creates the folder name and creates the name of each of the children files
    in:
        the full name of the mother file
    out:
        folder_name: the name of the folder name of the children
        split_file_name: the name of each of the children files
    =#
    file_position = findfirst('.', file_name)
    last_folder_position = findlast('/', file_name)

    folder_name = file_name[begin:last_folder_position - 1]
    split_file_name = "split_" * file_name[last_folder_position + 1 : file_position - 1] * '_'
    return folder_name, split_file_name
end


function line_splitter(infile::IOStream, file_name::String, lines_per_file::Int)
    #=
    Splits a file based on line count
    in:
        the file to use splitting on
        the name of the file that will be split
        the line count which should be in each of the files
    out:
        <Linecount of mother file / line_count> files containg exactly <lines_per_file> lines
    =#
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