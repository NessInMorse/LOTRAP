#=
FASTQ ANALYSER PLOT
---------
An interactive plot analysis for the fastq-file:
    creating a boxplot for each of the positions 
    of the read where the count of qualities > 1000.

Version 0.2 - added possibility to create plots

Date: 2023-05-26
Author: Marc Wijnands

Currently outputs to an example file.

NOTE:
Currently the program is running very slowly, since the plot is
creating A ton of lines in the process of making the plot.
=#
using PlotlyJS
# plot(a, x=:total_bill, kind="box")


function calculate_Qs(quality_per_position, count_per_position)
    #=
    Function that calculates the median, Q₀, Q₁, Q₃ & Q₄.
    And then creates the fastq plot showing the qualities
    in:
        quality_per_position: a dictionary containing the quality per position in
                a vector
        count_per_position: the counts per base in a vector
    out:
        A boxplot showing the spread of the quality per base.
    =#
    # count_per_position = count_per_position[1:50]
    box_per_position::Matrix = reshape([], (0, 0))
    
    for (i, count) in enumerate(count_per_position)
        if count > 1000
            running_sum = [0 for _ in quality_per_position[i]]
            for (j, sub_count) in enumerate(quality_per_position[i])
                if j > 1
                    running_sum[j] = running_sum[j-1] + sub_count
                else
                    running_sum[j] = sub_count
                end
            end
            println(running_sum)
            Q₀ = findfirst(x -> x > 0, running_sum)
            Q₄ = findfirst(x -> x == running_sum[end], running_sum)
            Q₂ = findfirst(x -> x >= (running_sum[end] ÷ 2), running_sum)
            Q₁ = findfirst(x -> x >= Q₂ ÷ 2, running_sum)
            Q₃ = findfirst(x -> x >= Q₂ + (Q₂ ÷ 2), running_sum)
            if Q₁ === nothing
                Q₁ = Q₂
            end
            if Q₃ === nothing
                Q₃ = Q₂
            end
            if isempty(box_per_position)
                box_per_position = reshape([Q₀, Q₁, Q₂, Q₃, Q₄], (5, 1))
            else
                box_per_position = hcat(box_per_position, reshape([Q₀, Q₁, Q₂, Q₃, Q₄], (5, 1)))
            end
        else
            break
        end
    end
    savefig(plot(box_per_position, x=:total_bill, kind="box"), "qualities.png")
    # println(box_per_position)
end


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
        quality_per_position, count_per_position = readfastq(infile)
        calculate_Qs(quality_per_position, count_per_position)
    end
end
main()