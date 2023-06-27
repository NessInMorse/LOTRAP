#=
FASTQ analyser
------
A simple FASTQ-analyser which will show insight into general facts of the data.

Version 1.2 - Added timing of the script
Date 2023-06-27
Author: Marc Wijnands

PLANNED FUNCTIONALITY:
A boxplot per position visualising the quality of each base in the reads
=#


function calculateAverage(total, groups)
        #=
        Calculates the average of a total among the split groups
        in:
                a total
                the groupcount
        out:
                A float containing the average of the total & groups
        =#
        average::Float64 = total / groups
        return average
end

function calculateQuality(line::String)
        #=
        Returns the sum of the quality of the DNA 
        using the 33-Phred score.
        in:
                the line containing the quality in ASCII-encoding
        out:
                the sum of the quality of the quality string
        =#
        return sum([codepoint(i) - 33 for i in line])
end


function calculateGC(totalGC::Int64,
                     gcCountPerPos::Vector{Int64}, 
                     totalCountperpos::Vector{Int64},
                     line::String,
                     lineLen::Int64)
        #=
        Calculates the average GC and the GC per position
        in:
                totalGC: the current total GC to be updated
                gcCountPerPos: the GC count per position
                totalCountperpos: the total count of bases per position
                line: the DNA String
                lineLen: the length of the line
        out:
                updated gcCountPerPos: containing the count of GC per position
                updated totalCountperpos: containing the total count of bases per position
                updated totalGC: containing the total count of GC
        =#
        offset::Int64 = lineLen - length(gcCountPerPos)

        if length(gcCountPerPos) == 0
                gcCountPerPos = [0 for i in 1:lineLen]
                totalCountperpos = [1 for i in 1:lineLen]
        elseif offset > 0
                addon::Vector{Int64} = [0 for i in 1:offset]
                append!(gcCountPerPos, addon)
                append!(totalCountperpos, addon)
        end
        for i in 1:lineLen
                if line[i] == 'G' || line[i] == 'C'
                        gcCountPerPos[i] += 1
                        totalGC += 1
                end
                totalCountperpos[i] += 1
        end
        return gcCountPerPos, totalCountperpos, totalGC
end


function analyseFile(infile::IOStream)
        #=
        Analyses a fastq-file on quality, read length, GC-percentage, 
        GC percentage per position.
        in:
                the fastq-file to analyse
        out:
                the count of reads
                the average length of the reads
                the minimum length of the reads
                the maximum length of the reads
                the average GC-percentage of the reads
                the GC-average per position (in %)
                the average quality of the reads
        =#
        i::Int8 = 0 
        count::Int64 = 0

        lineLen::Int64 = 0
        totalLen::Int64 = 0
        minLen::Int64 = 0
        maxLen::Int64 = 0
        totalQuality::Int64 = 0

        totalGC::Int64 = 0
        gcCountPerPos::Vector{Int64} = Vector{Int}(undef, 0)
        totalCountperpos::Vector{Int64} = Vector{Int}(undef, 0)
        readline(infile)
        for line in eachline(infile)
                lineLen = length(line)
                totalLen += lineLen
                minLen = minimum((minLen, lineLen))
                maxLen = maximum((maxLen, lineLen))
                gcCountPerPos, totalCountperpos, totalGC = calculateGC(totalGC, gcCountPerPos, totalCountperpos,  line, lineLen)
                count += 1
                readline(infile)
                totalQuality += calculateQuality(readline(infile))


                readline(infile)
        end
        close(infile)
        averageLen::Float64 = round(calculateAverage(totalLen, count), digits=2)

        averageGC::String = "$(round(calculateAverage(totalGC, totalLen) * 100, digits=2))%"
        gcAvgPerPos::Vector{String} = ["$(round(calculateAverage(gcCountPerPos[item[1]],
                                          totalCountperpos[item[1]]) * 100, digits=2))%"
                                         for item in enumerate(totalCountperpos)]

        avgQuality::Float64 = totalQuality / totalLen
        return count, averageLen, minLen, maxLen, averageGC, gcAvgPerPos, avgQuality
end


function writeFile(filename::String,
                   count::Int64,
                   averageLen::Float64,
                   minLen::Int64,
                   maxLen::Int64,
                   averageGC::String,
                   gcAvgPerPos::Vector{String},
                   avgQuality::Float64)
        #=
        Function to write the data of the fastq-file
        to an output file
        in:
                filename: the filename of the output file
                count: the total counts of reads in the fastq-file
                averageLen: the average length of the reads
                minLen: the minimum length of the reads
                maxLen: the maximum length of the reads
                averageGC: the average GC-percentage in the reads
                gcAvgPerPos: the average GC-percentage per position in the reads
                AvgQuality: the average quality of every position in the read
        out:
                a file with all the data in a readable manner
        =#
        infile = open("$filename", "w")
        write(infile, "Read count:\t\t\t $count\n" *
                       "Average length of reads:\t $averageLen\n" *
                       "Minimum length of reads:\t $minLen\n" *
                       "Maximum length of reads:\t $maxLen\n" *
                       "Average quality of reads:\t $avgQuality\n" *
                       "Average GC%:\t\t\t $averageGC\n" *
                       "Average GC% per Position:\t $gcAvgPerPos")
        close(infile)
end


function openFile(filename::String)
        #=
        Opens a file and returns the IOStream of that file
        in:
                a filename of the input file
        out:
                an IOStream containing the file
        =#
        infile::IOStream = open("$filename", "r")
        return infile
end



function main()
        filename::String = ARGS[1]
        outfile::String = ARGS[2]
        infile::IOStream = openFile(filename)
        count::Int64, averageLen::Float64,
                minLen::Int64, maxLen::Int64,
                averageGC::String, gcAvgPerPos::Vector{String},
                avgQuality::Float64  = analyseFile(infile)
        writeFile(outfile, 
                  count, 
                  averageLen, 
                  minLen, 
                  maxLen, 
                  averageGC, 
                  gcAvgPerPos, 
                  avgQuality)
end


function time_main()
        time = @elapsed main()
        time_name = split(ARGS[3], '/')[end]
        outfile = open(ARGS[3], "w")
        write(outfile, "$(time_name)\t$(time)")
        close(outfile)
end

time_main()
