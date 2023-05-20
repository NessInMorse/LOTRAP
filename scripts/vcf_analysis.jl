
function writeFile(outname::String,
                   scores::String,
                   insertcount::Int64,
                   delcount::Int64,
                   indelscore::Float64)
        outfile = open(outname, "w")
        write(outfile, "Insert count: $insertcount\n" *
                       "Deletion count: $delcount\n" *
                       "INDEL-score; Insertion : deletion,  1 : $indelscore\n" *
                       "Mutations:\n $scores")
        close(outfile)
end

function returnPrintableGraph(scores::Vector{Vector{Int64}}, titles::Vector{String})        
        printstr::String = ""
        printstr = printstr * "\t\t\t\t" * "FROM\n\t" * "-" ^ 50 * "\n"
        printstr = printstr * "\t\t|\t" * join(titles, "\t") * "\n" * "\t\t|\t" * "-" ^ 34 * "\n"
        timely::String = ""
        for i in 1:length(scores)
                
                if i == length(scores) / 2
                        timely = "\tTO\t| $(titles[i])"
                else
                        timely = "\t\t| $(titles[i])"
                end
                for j in 1:length(scores)
                        timely =  timely * "\t" * join(scores[j][i], "")
                end
                timely = timely * "\n"
                printstr = printstr * timely
        end

        return printstr
end


function returnIndex(base::SubString{String})
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


function createGraph(titles::Vector{String})
        default::Vector{Int64} = [0 for i in 1:length(titles)]
        scores::Vector{Vector{Int64}} = []
        for i in titles
                push!(scores, copy(default))
        end

        return scores
end

function readFile(infile::IOStream, titles::Vector{String})
        delcount::Int64 = 0
        insertcount::Int64 = 0
        scores::Vector{} = createGraph(titles)
        items::Vector{SubString{String}} = []
        alt::Vector{SubString{String}} = []
        for line::String in eachline(infile)
                if !(line[1] == '#')
                        items = split(line, "\t")
                        ref = items[4]
                        alt = split(items[5], ",")
                        for alteration in alt
                                 if alteration == "<*>" # I believe this to be a deletion
                                        alteration = ""
                                end
                                if ref != "N" # What could we possible do with aNy's
                                        if length(ref) == 1 && length(alteration) == 1
                                                firstpos = returnIndex(ref)
                                                secondpos = returnIndex(alteration)
                                                scores[firstpos][secondpos] += 1
                                        elseif length(ref) < length(alteration)
                                                insertcount += 1
                                        elseif length(ref) > length(alteration)
                                                delcount += 1
                                        end
                                end
                        end

                end
        end
        close(infile)
        return scores, insertcount, delcount

end



function openFile(filename::String)
        infile = open("$filename", "r")
        return infile
end


function main()
        filename::String = ARGS[1]
        outname::String = ARGS[2]
        infile::IOStream = openFile(filename)

        titles::Vector{String} = ["A", "C", "G", "T"] # We remember this sequence
        scores, insertcount, delcount = readFile(infile, titles) # delcount, insertcount, indel_ratio, poss_mutations = 
        indelscore::Float64 = 1 / (insertcount / delcount)
        print_scores = returnPrintableGraph(scores, titles)
        
        writeFile(outname, print_scores, insertcount, delcount, indelscore)
end
if length(ARGS) == 2 && ARGS[1] != ARGS[2]
        main()
end