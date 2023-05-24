#=
VCF PLOT
--------
A script that creates an interactive plot of the data in the vcf-file
    containing all the SNP-mutation types per contig in a heatmap

VERSION: 0.2 - Added DataFrames package
Date: 2023-05-24
Author: Marc Wijnands

Added basic functionality to the script.
Is not yet able to produce plots
Was also not able to test functionality yet.

Planned to let it create plots using Plotly(JS)
=#
try
        using DataFrames
        using PlotlyJS
catch
        using Pkg
        Pkg.add("DataFrames")
        Pkg.add("PlotlyJS")
finally
        using DataFrames
        using PlotlyJS
end


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


function readFile(infile::IOStream, titles::Vector{String})::Dict{String, Matrix{Int}}
    #=
    Reads a VCf-file and adds the counts to a Dictionary containing all the SNPs per contig.
    N and structural variants <*> are deleted.
    in:
            input IOStream
            the titles of the headers of the 2D vector
    out:
            The dictionary containing the SNP counts per contig
    =#
    all_items::Dict{String, Matrix{Int}} = Dict() 
    items::Vector{SubString{String}} = []
    alt::Vector{SubString{String}} = []
    for line::String in eachline(infile)
            if !(line[1] == '#')
                    items = split(line, "\t")
                    contig = items[1]
                    if get(all_items, contig, 0) == 0 # if the contig can not be found in the contig, add said contig to dictionary
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
    return all_items
end


function main()
    if length(ARGS) == 2
        infile = open(ARGS[1], "r")
        # outfile = open(ARGS[2], "w")
        titles::Vector{String} = ["A", "C", "G", "T"] # We remember this sequence
        all_contigs = readFile(infile, titles)
        all_dfs::Dict{String, DataFrame} = Dict()
        for key in keys(all_contigs)
                df = DataFrame(all_contigs[key], :auto)
                rename!(df, combo)
                df[!, "combination"] = combo
                
                all_dfs[key] = df
        end
        
        # println(all_dfs)
        # println(length(all_dfs))
        # close(outfile)
    end
end


main()