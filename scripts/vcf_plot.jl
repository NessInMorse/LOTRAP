#=
VCF PLOT
--------
A script that creates an interactive plot of the data in the vcf-file
    containing all the SNP-mutation types per contig in a heatmap

VERSION: 1.0 - Added 16 heatmaps per batch
Date: 2023-05-24
Author: Marc Wijnands

Currently is able to create 16 heatmaps per file
and has the remainder in an extra file
Also removed the ticklabels since they do not add useful information.

currently unable to see the x-/y-axis titles
=#
try
        using PlotlyJS
catch
        using Pkg
        Pkg.add("PlotlyJS")
finally
        using PlotlyJS
end

function create_heatmaps(all_contigs::Dict{String, Matrix{Int}}, outfolder)
        #=
        Function that creates heatmaps of all the vcf-data per contig
        in:
                all_contigs: all SNP counts per contig
                outfolder: the output folder
        out:
                several heatmap plots containing the SNP counts per contig
                        and the severity of the counts
        =#
        heatmaps = [heatmap() for i in 1:16]
        contig_names = ["" for i in 1:16]
        plot_count::Int = 1
        contig_count::Int = 1
        for i in keys(all_contigs)
                layout = Layout(title="Directed SNP mutations in contig $(i)", 
                                xaxis_title = "Wild type",
                                yaxis_title = "Mutant type")
                data = all_contigs[i]
                hm = heatmap(
                        x = ["A", "C", "G", "T"],
                        y = ["A", "C", "G", "T"],
                        z = data,
                        colorbar_title = "Relative intensity",
                        colorscale = "Viridis",
                        colorbar_showticklabels = false,
                        xaxis_title = "FROM",
                        yaxis_title = "TO"
                )
                heatmaps[(contig_count-1) % 16 + 1] = hm
                contig_names[(contig_count-1) % 16 + 1] = i
                
                if iszero(contig_count % 16)
                        sub_plotties = make_subplots(rows=4, cols=4,
                                    specs=reshape([Spec() for i in 1:16], (4, 4)),
                                    subplot_titles=reshape(contig_names, (4, 4)),
                                    x_title = "FROM",
                                    y_title = "TO")

                        for (i, h) in enumerate(heatmaps)
                                c = ((i-1) % 4) + 1
                                r = ((i-1) ÷ 4) + 1
                                add_trace!(sub_plotties, h, row=r, col=c)
                        end
                        
                        relayout!(sub_plotties, title_text="VCF's of 16 contigs")
                        savefig(sub_plotties, outfolder * "heatmap_$(plot_count).png")
                        plot_count += 1
                        heatmaps = [heatmap() for i in 1:16]
                        contig_names = ["" for i in 1:16]
                end
                contig_count += 1
        end
        if !iszero(contig_count % 16)
                sub_plotties = make_subplots(rows=4, cols=4,
                specs=reshape([Spec() for i in 1:16], (4, 4)),
                subplot_titles=reshape(contig_names, (4, 4)))

                for (i, h) in enumerate(heatmaps)
                        c = ((i-1) % 4) + 1
                        r = ((i-1) ÷ 4) + 1
                        add_trace!(sub_plotties, h, row=r, col=c)
                end
                
                relayout!(sub_plotties, title_text="VCF's of 16 contigs")
                savefig(sub_plotties, outfolder * "heatmap_$(plot_count).png")
        end

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
        outfolder = ARGS[2]
        titles::Vector{String} = ["A", "C", "G", "T"] # We remember this sequence
        all_contigs = readFile(infile, titles)
        create_heatmaps(all_contigs, outfolder)
        
        # println(all_dfs)
        # println(length(all_dfs))
        # close(outfile)
    end
end


main()