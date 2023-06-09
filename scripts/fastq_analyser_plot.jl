#=
FASTQ ANALYSER PLOT
---------
An interactive plot analysis for the fastq-file:
    creating a boxplot for each of the positions 
    of the read where the count of qualities > 1000.

Version 1.1 - Added timing of the script

Date: 2023-06-27
Author: Marc Wijnands

NOTE:
It takes around 5 minutes to plot a file of about 9Gb of data.

PLANNED:
A GC% Distribution barplot showing the distribution of the GC percentages of the reads

=#
try
    using PlotlyJS
catch
    using Pkg
    Pkg.add("PlotlyJS")
    using PlotlyJS
end


function plot_gc_distribution(gc_distribution::Vector{Int}, gc_file::String)
    trace = bar(x = 1:length(gc_distribution), y = gc_distribution, marker_color = "#3f167e")
    layout = Layout(title="GC-percentage distribution in sample")
    p = plot([trace], layout)
    savefig(p, "$(gc_file).html")
    savefig(p, "$(gc_file).png")
end


function calculate_Qs(quality_per_position, count_per_position, out_file)
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
            Q₀ = findfirst(x -> x > 0, running_sum)
            Q₄ = findfirst(x -> x == running_sum[end], running_sum)
            Q₂ = findfirst(x -> x >= (running_sum[end] ÷ 2), running_sum)
            Q₁ = findfirst(x -> x >= running_sum[Q₂] ÷ 2, running_sum)
            Q₃ = findfirst(x -> x >= ((running_sum[Q₂] + running_sum[Q₄]) ÷ 2), running_sum)
            if Q₁ === nothing
                Q₁ = Q₂
            end
            if Q₃ === nothing
                Q₃ = Q₂
            end
            if isempty(box_per_position)
                box_per_position = [Q₀ Q₁ Q₂ Q₃ Q₄]
            else
                box_per_position = vcat(box_per_position, [Q₀ Q₁ Q₂ Q₃ Q₄])
            end
        else
            break
        end
    end
    p = plot(box_per_position, 
                    Layout(title="Read quality of FASTQ-file",
                           yaxis_title="Phred Score",
                           xaxis_title="Position"))
    restyle!(p, 1, marker_color="rgba(255, 108, 108, 0.8)", name=:Q₀)
    restyle!(p, 2, marker_color="rgba(255, 243, 107, 0.8)", name=:Q₁)
    restyle!(p, 3, marker_color="rgba(119, 255, 107, 0.8)", name=:Q₂)
    restyle!(p, 4, marker_color="rgba(105, 255, 232, 0.8)", name=:Q₃)
    restyle!(p, 5, marker_color="rgba(151, 105, 255, 0.8)", name=:Q₄)
    savefig(p, "$(out_file).html")
    savefig(p, "$(out_file).png")
    # println(box_per_position)
end


function calculateGC(gc_distribution::Vector{Int}, read::String)::Vector{Int}
    #=
    Function that creates the GC-percentage of a read
        and adds this to the distribution vector
    in:
        a distribution of gc-percentages of the reads
        a read
    out:
        an updated distribution with the new GC-count of the read
    =#
    gc_count = sum([i == 'G' || i == 'C' for i in read])
    total = length(read)
    gc_distribution[Int(floor(((gc_count / total) * 100) + 0.5)) + 1] += 1
    return gc_distribution
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
        quality = (iszero(quality) * 1) + (!iszero(quality)) * quality
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
    gc_distribution::Vector{Int} = fill(0, 101)
    read::String = ""
    for line in eachline(infile)
        read = readline(infile)
        gc_distribution = calculateGC(gc_distribution, read)
        readline(infile)
        line = readline(infile)
        quality_per_position, count_per_position = analyse_quality(line, 
                                                                   quality_per_position,
                                                                   count_per_position)
    end
    return quality_per_position, count_per_position, gc_distribution
end

function main()
    infile::IOStream = open(ARGS[1], "r")
    out_file::String = ARGS[2]
    gc_file::String = ARGS[3]
    quality_per_position, count_per_position, gc_distribution = readfastq(infile)
    calculate_Qs(quality_per_position, count_per_position, out_file)
    plot_gc_distribution(gc_distribution, gc_file)
end

function time_main()
    time = @elapsed main()
    time_name = split(ARGS[4], '/')[end]
    outfile = open(ARGS[4], "w")
    write(outfile, "$(time_name)\t$(time)")
    close(outfile)
end

time_main()
