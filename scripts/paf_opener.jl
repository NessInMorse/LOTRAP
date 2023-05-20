

function open_paf(filename)
    infile = open(filename, "r")
    sequences::Vector{String} = []
    count_per_sequence::Dict{String, Int} = Dict()

    for line in eachline(infile)
        fields = split(line, '\t')
        
        count_per_sequence[fields[6]] = get(count_per_sequence, fields[6], 0) + 1
        if !(fields[6] in sequences)
            push!(sequences, fields[6])
        end
    end
    
    close(infile)
    return sequences, count_per_sequence
end