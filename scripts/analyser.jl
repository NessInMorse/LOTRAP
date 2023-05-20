function calculateAverage(total, groups)
        average::Float64 = total / groups
        return average
end

function calculateQuality(line::String)
        return sum([codepoint(i) - 33 for i in line])
end


function calculateGC(totalGC::Int64,
                     gcCountPerPos::Vector{Int64}, 
                     totalCountperpos::Vector{Int64},
                     line::String,
                     lineLen::Int64)
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

function checkIfLowest(minLen::Int64, lineLen::Int64)
        if minLen == 0 || lineLen < minLen
                return lineLen
        else
                return minLen
        end
end

function checkIfHighest(maxLen::Int64, lineLen::Int64) 
        if maxLen == 0 || lineLen > maxLen
                return lineLen
        else
                return maxLen
        end

end

function analyseFile(infile::IOStream)
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
                minLen = checkIfLowest(minLen, lineLen)
                maxLen = checkIfHighest(maxLen, lineLen)
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

@time main()