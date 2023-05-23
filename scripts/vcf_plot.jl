#=
VCF PLOT
--------
A script that creates an interactive plot of the data in the vcf-file
    containing all the SNP-mutation types per contig in a heatmap

VERSION: 0.1
Date: 2023-05-24
Author: Marc Wijnands

Added basic functionality to the script.
Is not yet able to produce plots
Was also not able to test functionality yet.

Planned to let it create plots using Plotly(JS)
=#

function returnIndex(base::SubString{String})
    #=
    Returns the index of a base
    A -> 1
    C -> 2
    G -> 3
    T -> 4
    =#
    if base == "A"
            return 1
    elseif base == "C"
            return 2
    elseif base == "G"
            return 3
    else
            return 4
    end
end


function readFile(infile::IOStream, titles::Vector{String})::Dict{String, Vector{Vector{Int}}}
    #=
    Reads a VCf-file and adds the counts to a Dictionary containing all the SNPs per contig.
    N and structural variants <*> are deleted.
    in:
            input IOStream
            the titles of the headers of the 2D vector
    out:
            The dictionary containing the SNP counts per contig
    =#
    all_items::Dict{String, Vector{Vector{Int}}} = Dict()
    
    items::Vector{SubString{String}} = []
    alt::Vector{SubString{String}} = []
    for line::String in eachline(infile)
            if !(line[1] == '#')
                    items = split(line, "\t")
                    contig = items[1]
                    if !(get(all_items, contig, true)) # if the contig can not be found in the contig, add said contig to dictionary
                        all_items[contig] = zeros(Int32, (4, 4))
                    end
                    ref = items[4]
                    alt = split(items[5], ",")
                    for alteration in alt
                            if alteration == "<*>" # I believe this to be a deletion
                                    alteration = ""
                            end
                            if ref != "N" # What could we possible do with aNy's
                                    if length(ref) == 1 && length(alteration) == 1
                                            WT = returnIndex(ref)
                                            MT = returnIndex(alteration)
                                            all_items[contig][((WT-1) * 4) + MT] += 1
                                    end
                            end
                    end

            end
    end
    close(infile)
    return scores
end


function main()
    if length(ARGS) == 2
        infile = open(ARGS[1], "r")
        outfile = open(ARGS[2], "w")
        titles::Vector{String} = ["A", "C", "G", "T"] # We remember this sequence
        scores = readFile(infile, titles)
        close(infile)
        close(outfile)
    end
end


main()