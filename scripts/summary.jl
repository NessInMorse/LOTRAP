#=
Time summary script
---------
Summarizes the output of multiple tsv-files in a folder.
Into a single file containing all the times.

Version 1.0 - Basic functionality with text

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
            write(out_file, "$(lines[begin])s\n")
            total_time += parse(Float64, split(lines[begin], '\t')[end])
        end
        close(infile)
    end
    write(out_file, "sum\t$(total_time)\n")
    close(out_file)
end

main()