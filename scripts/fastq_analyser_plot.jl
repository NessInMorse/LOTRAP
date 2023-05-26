#=
FASTQ ANALYSER PLOT
---------
An interactive plot analysis for the fastq-file:
    creating a boxplot for each of the positions 
    of the read where the count of qualities > 1000.

Version 0.1
Date: 2023-05-26
Author: Marc Wijnands
=#

function analyse_quality(quality_string::String, 
                         quality_per_position::Dict{Int, Vector{Int}},
                         count_per_position::Vector{Int})
    #=
    Analyses the quality string and updates 
        the quality per position, and
        the count per position
    in:
        quality_string: the line / string with the quality
        quality_per_position: A dictionary containing the qualities per position in vectors
        count_per_position: The count of each position in a vector
    =#
    qualities = [Int(i) - 33 for i in quality_string]
    for (i, quality) in enumerate(qualities)
        if !(i in keys(quality_per_position))
            quality_per_position[i] = zeros(Int, 100)
            push!(count_per_position, 0)
        end    
        quality_per_position[i][quality] += 1
        count_per_position[i] += 1
    end
    return quality_per_position, count_per_position
end

function readfastq(infile::IOStream)
    #=
    Reads a fastq file and retrieves the qualities & counts in two dictionaries
    in:
        infile: a file to read
    out:
        quality_per_position: a dictionary containing vectors of the qualities
        count_per_position: a dictionary containing the counts of each position


    =#
    quality_per_position::Dict{Int, Vector{Int}} = Dict()
    count_per_position::Vector{Int} = []
    for line in eachline(infile)
        readline(infile)
        readline(infile)
        line = readline(infile)
        quality_per_position, count_per_position = analyse_quality(line, 
                                                                   quality_per_position,
                                                                   count_per_position)
    end
    return quality_per_position, count_per_position
end

function main()
    if length(ARGS) == 2
        infile::IOStream = open(ARGS[1], "r")
        out_folder::String = ARGS[2]
        readfastq(infile)
    end
end
main()