#=
Time summary script
---------
Summarizes the output of multiple tsv-files in a folder.
Into a single file containing all the times.

Version 1.1 - Added rounding for better readability

Date: 2023-06-27
Author: Marc Wijnands

NOTE:
This script is not timed, since it takes an insignificant amount of time
=#



function main()
    files = readdir(ARGS[1])
    prefix = ARGS[1]
    out_file = open(ARGS[2], "w")
    total_time = 0
    for i in files
        infile = open(prefix * i, "r")
        lines = readlines(infile)
        if length(lines) == 1
            filename = split(lines[begin], '.')[begin] * ".jl"
            file_time = round(parse(Float64, split(lines[begin], '\t')[2]), digits=3)
            write(out_file, "$(filename)\t$(file_time)s\n")
            total_time += parse(Float64, split(lines[begin], '\t')[end])
        end
        close(infile)
    end
    total_time = round(total_time, digits=3)
    write(out_file, "sum\t$(total_time)s\n")
    close(out_file)
end

main()