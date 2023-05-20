using Plots

function plot_paf(paf_file)
    query_starts = []

    infile = open(paf_file, "r")

    for line in eachline(infile)
        fields = split(line, '\t')
        

    end
    close(infile)

    scatter(query_starts, target_starts, color="blue")
    scatter(query_ends, target_ends, color="red")
end

plot_paf("/home/ness/Programming/Julia/pipelines/input/head_enhydra_lutris.paf")