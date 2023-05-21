#=
VCF_analysis
----
This script does a simple VCf-analyses to check the SNP-mutations
and the INDEL mutations. To see the difference between the
formed asssembly from the reference and the reference itself.

Version: 1.0
Date: 2023-05-21
Author: Marc Wijnands

PLANNED FUNCTIONALITY:
A interactive heatmap and barplot showing the relationship
of the SNP mutations and the INDELscores respectively.
=#


function writeFile(outname::String,
                   scores::String,
                   insertcount::Int64,
                   delcount::Int64,
                   indelscore::Float64)
        #=
        A function that writes the outputfile
        in:
                - the name of the output file
                - the scores of the mutations
                - the insert count
                - the deletion count
                - the indelscore comparing insertions and deletions
        out:
                a file with the name <OUTNAME> with all the analysed data
        =# 
        outfile = open(outname, "w")
        write(outfile, "Insert count: $insertcount\n" *
                       "Deletion count: $delcount\n" *
                       "INDEL-score; Insertion : deletion,  1 : $indelscore\n" *
                       "Mutations:\n $scores")
        close(outfile)
end

function returnPrintableGraph(scores::Vector{Vector{Int64}}, titles::Vector{String})
        #=
        Returns a printable graph that is humanly readable
        in:     
                a 2D vector containing all the scores of the
                        SNPs
                The header titles of the 2D vector
        out:
                A string that can be printed to show 
                        the score in a visually pleasing matter
        =#
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


function createGraph(titles::Vector{String})
        #=
        Creates a 2D vector to hold all the scores
                of the mutations of the vcf
        in:
                a vector containing the titles of the headers of the
                        2D vector
        out:
                a 2D vector filled with 0
        =#
        default::Vector{Int64} = [0 for i in 1:length(titles)]
        scores::Vector{Vector{Int64}} = []
        for i in titles
                push!(scores, copy(default))
        end

        return scores
end

function readFile(infile::IOStream, titles::Vector{String})
        #=
        Reads a VCf-file and adds the counts to a 2D vector.
        Only counts the mutations with INDEL or a ACTG SNP.
        So N and structural variants <*> are deleted.
        in:
                input IOStream
                the titles of the headers of the 2D vector
        out:
                The 2D vector containing the counts
                the insertcount
                the deletion count
        =#
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
        #=
        Function that opens a file and returns the IOStream
        in:
                the name of the file
        out:
                the IOStream of the input file
        =#
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