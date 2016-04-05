module ResonatorFit

import Base: complex


using DataArrays, DataFrames, Optim, JLD, PyPlot, Colors

export s2p, complex2dB, complex2Phase, dB2mag, s21Fit, savitsky_golay, s21FitAuto,
        bestFit, searchdir, s21FitPowerSweepJLD, s21FitTimeSweepJLD, plotAllTracesJLD,
        plotAllResultsPowerSweepJLD, plotAllResultsTimeSweepJLD, s21FitFunc, s21FitJLD,
        mergeJLDFileResults, mag2dB, feedOutArrayFromDataFrame, complexMapSNP,
        getQFromFitParams, fitGrade, cableAtten, ceifk, vnaPower2Photon, circleFit

function circleImag(phase, center, r)
    return r*cos(phase) + center[1] + 1im*(r*sin(phase) + center[2])
end

function phaseMatch(freq, s21, phase, circle)
    dist = zeros(length(freq),length(phase))
    phaseMatched = zeros(length(freq))
    for i in 1:length(freq)
        for j in 1:length(phase)
            dist[i,j] = real(s21[i] - circle[j])^2 + imag(s21[i] - circle[j])^2
        end
        phaseMatched[i]  = phase[indmin(dist[i,:])]
    end
    return phaseMatched
end

function bandwidthLine(x, f0idx, s21)
    return -real(s21[f0idx])/imag(s21[f0idx])*x
end

function f0Line(x, f0idx, s21)
    return imag(s21[f0idx])/real(s21[f0idx])*x
end

function lineIntersectIndex(x1, x2, y1, y2)
    dist = zeros(length(x2),length(x1))
    for i in 1:length(x2)
        for j in 1:length(x1)
            dist[i,j] = (x2[i] - x1[j])^2 + (y2[i] - y1[j])^2
        end
    end
    (min2, min1) = [mod(indmin(dist),length(x2)), convert(Int64,floor(indmin(dist)/length(x2)))]
    return min1, min2
end
# get quality factor for circle fit indicies
function circleFitQi(freq, f0idx, minPhaseN, minPhaseP)
    return freq[f0idx]/abs(freq[minPhaseP] - freq[minPhaseN])
end
# Qi extraction based off of arXiv:1201.3384
function circleFit(freq, s21; fitPts = 101, plotFit = true)
    phase = linspace(0, 2*pi, fitPts)
    f0idx = indmin(complex2dB(s21))
    s21inv = 1./s21
    xDim = maximum(real(s21inv)) - minimum(real(s21inv))
    yDim = maximum(imag(s21inv)) - minimum(imag(s21inv))
    r = maximum([xDim, yDim])/2
    center = [maximum(real(s21inv)) - xDim/2, maximum(imag(s21inv)) - yDim/2]
    s21inv = s21inv - (center[1] + 1im*center[2]) # Shift circle to orgin
    circle = circleImag(phase,[0,0] , r)
    phaseMatched = phaseMatch(freq, s21inv, phase, circle)
    cMatch = circleImag(phaseMatched, [0,0], r)

    xN = linspace(-r, 0, fitPts)
    xP = linspace(0, r, fitPts)
    bN = bandwidthLine(xN, f0idx, s21inv)
    bP = bandwidthLine(xP, f0idx, s21inv)
    fN = f0Line(xN, f0idx, s21inv)
    fP = f0Line(xP, f0idx, s21inv)

    minXNidx, minPhaseNidx = lineIntersectIndex(xN, real(cMatch), bN, imag(cMatch))
    minXPidx, minPhasePidx = lineIntersectIndex(xP, real(cMatch), bP, imag(cMatch))
    QiApprox = circleFitQi(freq, f0idx, minPhaseNidx, minPhasePidx)
    if plotFit
        figure()
        plot(real(s21inv), imag(s21inv))
        plot(real(cMatch), imag(cMatch))
        plot(xN, bN)
        plot(xP, bP)
        plot(xN, fN)
        plot(xP, fP)
        plot(real(s21inv[f0idx]), imag(s21inv[f0idx]), "ro")
        plot(real(s21inv[minPhaseNidx]), imag(s21inv[minPhaseNidx]), "go")
        plot(real(s21inv[minPhasePidx]), imag(s21inv[minPhasePidx]), "go")
        axis([-r; r; -r; r])
        title(string("f0 = ",freq[f0idx],"\nQi = ", QiApprox))
        xlabel("Re 1/S21")
        ylabel("Im 1/S21")
        # legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    end
    return QiApprox
end


"""
    searchdir(path,key)

Search for files matching key in path.\n
__Arg(s)__:\n
`path` - path to search\n
`key` - keyword for file(s) to search\n
__Output(s)__:\n
`files` - list of files in path matching key\n
__Example__:\n
```
juliaFiles = searchdir(path,".jl")
```
"""

function searchdir(path,key)
    return filter(x->contains(x,key), readdir(path))
end

"""
    s2p(file::AbstractString)

Imports arrays of frequency and s parameters from a .s2p file\n
__Arg(s)__:\n
`file` - Touchstone file with extension .s2p\n
__Output(s)__:\n
`frequency`, `S11`, `S21`, `S12`, `S22` - arrays of freqeuency and two port s paramters\n
__Example__:\n
```
f, s11, s21, s12, s22 = s2p("touchstone.s2p")
```
"""
function s2p(file::AbstractString)

    complexNames = [:frequency, :S11, :S21, :S12, :S22]
    data = complex(
        readtable(
        file,
        skipstart = 8,
        header = false,
        separator = ' '
        ),
        complexNames, complexMapSNP(2))
    return feedOutArrayFromDataFrame(data)
end

# function complex2dB(sparam::DataArray)
#     return DataArray(20*log10(abs(Array(sparam))))
# end

"""
    complex(x::DataArray, y::DataArray)

Complex function written for the DataArray type.\n
__Arg(s)__:\n
`x` - DataArray of real parts\n
`y` - DataArray of imaginary parts\n
__Output(s)__:\n
`output` - complex mapped DataArray\n
__Example__:\n
```
complexRand = complex(DataArray(rand(10)), DataArray(1im*rand(10)))
```
"""
function complex(x::DataArray, y::DataArray)
    return DataArray(complex(Array(x), Array(y)))
end

"""
    feedOutArrayFromDataFrame(data::DataFrame)

Reformats DataFrame as a sequence of column arrays output for function ouput.\n
__Arg(s)__:\n
`data` - DataArray of real parts\n
__Output(s)__:\n
`arrayOfArrays` - 1d array containing arrays extracted from DataFrame
"""
function feedOutArrayFromDataFrame(data::DataFrame)
    arrayOfArrays = []
    for i in 1:size(data, 2)
        push!(arrayOfArrays, Array(data[:,i]))
    end
    return arrayOfArrays
end

"""
    complexMapSNP(N)

Complex function map for sNp dataframes. A row with one element equal to one (sum = 1)
corresponds to doing nothing (freqeuency component). A row with two elements (one
of which is equal to one (real part), and the other which is equal to negative two (imaginary part))
so that sum = -1, corresponds to the combination of two seperate complex elements into
one complex number.\n
__Arg(s)__:\n
`N` - integer for a sNp (N port) measurement\n
__Output(s)__:\n
`map` - 2d array
"""
# complex function map for sNp dataframes
function complexMapSNP(N)
    map = round(Int,zeros(1 + N^2, 1 + 2*N^2))
    map[1, 1] = 1 # freqeuncy part
    for i in 2:size(map, 1)
        for j in 2:size(map, 2)
            if j == 2*(i-1)
                map[i, j] = 1 # s param real part
                map[i, j+1] = -2 # s param imag part
            end
        end
    end
    return map
end

"""
    complex(data::DataFrame, outNames::Array{Symbol,1},  cMap::Array{Int, 2})

Complex function written for the DataFrame type.\n
__Arg(s)__:\n
`data` - input of the DataFrame type\n
`outNames` - array of symbols corresponding output data column labels\n
`cMap` - complex mapping for sNp files, generated by ResonatorFit.complexMapSNP(N)\n
__Output(s)__:\n
`dataFrame` - complex mapped DataFrame\n
__Example__:\n
```
complexNames = [:frequency, :S11, :S21, :S12, :S22]
dfOut = complex(dfIn, complexNames, complexMapSNP(2))
```
"""
function complex(data::DataFrame, outNames::Array{Symbol,1},  cMap::Array{Int, 2})
    inNames = names(data)
    # Dimension checking between map and data
    size(cMap, 2) !== size(data, 2)?
    error(string("Dismension mismatch between map and input data. ",
    "Column length of input data for specified map should be ",
    size(cMap, 2), ". Current length is " , size(data, 2))): nothing
    # Dimension checking between map and output names
    size(cMap) !== (maximum(size(outNames)), maximum(size(inNames)))?
    error(string("Dismension mismatch between map and names. Name length for input data should be ",
    size(cMap, 2), ". Current value is " , maximum(size(inNames)) ,
    ". Name length for output data should be ", size(cMap, 1),
    ". Current value is " , maximum(size(outNames)), "."  )) : nothing

    dataFrame = DataFrame()
    nameAssign = [:real, :imag]

    for i in 1:size(cMap,1)
        # frequency part of map, no complex operation
        if countnz(cMap[i,:]) == 1 || countnz(cMap[i,:]) == 2
            if sum(sum(cMap[i,:], 2)) == 1
                for j in 1:size(cMap,2)
                    if cMap[i,j] == 1
                        dataFrame[outNames[i]] = data[inNames[j]]
                    end
                end
            # s parameter part of map, do complex operation
            elseif sum(sum(cMap[i,:], 2)) == -1
                for j in 1:size(cMap,2)
                    if cMap[i,j] == 1
                        nameAssign[1] = inNames[j]
                    elseif cMap[i,j] == -2
                        nameAssign[2] = inNames[j]
                    end
                 end
                dataFrame[outNames[i]] = complex(data[nameAssign[1]], data[nameAssign[2]])
            # unexpected map format
            else
                error("Cannot create a complex number from map! ",
                "Check that there is at least one integer with value 1 ",
                "(real part/no conversion) and optionally one integer with",
                " value -2 (imag part/ conversion) in each row.")
            end
        else
            error("There should only be one or two nonzero elements per row! ")
        end
    end
    return dataFrame
end

"""
    complex2dB(x)

Converts a complex number to a magnitube in dB.\n
__Arg(s)__:\n
`x` - complex number or array\n
__Output(s)__:\n
`db` - magnitude in dB\n
__Example__:\n
```
db = complex2dB([1 + 2im;1 - 2im])
```
"""
function complex2dB(x)
    return 20*log10(abs(x))
end

"""
    complex2Phase(x)

Converts a complex number to a phase relative to the x axis.\n
__Arg(s)__:\n
`x` - complex number or array\n
__Output(s)__:\n
`ph` - phase relative to the x axis\n
__Example__:\n
```
ph = complex2Phase([1 + 2im;1 - 2im])
```
"""
function complex2Phase(x)
    return atan2(imag(x), real(x))
end

"""
    dB2mag(db)

Converts from dB to magnitude.\n
__Arg(s)__:\n
`db` - dB input\n
__Output(s)__:\n
`mag` - magnitude output
"""
function dB2mag(db)
    return 10^(db/20)
end
"""
    mag2dB(mag)

Converts from magnitude to dB.\n
__Arg(s)__:\n
`mag` - mag input\n
__Output(s)__:\n
`db` - dB output
"""
function mag2dB(mag)
    return 20*log10(mag)
end

# decodes the fitCode as a smoothingFactor and a fitFreqOffset flag
function tryFit(fitCode)
    decoded = readdlm(IOBuffer(fitCode))
    sf = decoded[1]
    ffo = (decoded[2] != 0)
    return sf, ffo
end
# takes the best fit from the DataFrame output
function bestFit(fitResults::DataFrame)
    sort!(fitResults, cols = [:score])
    best = fitResults[1,:]
    # f0, Qi, score, grade = (float(fitResults[:f0][1]), float(fitResults[:Qi][1]), float(fitResults[:score][1]), fitResults[:grade][1])
    return best
end

function assignFromLabel(label::Symbol, dataframe::DataFrame)
    # $label is a symbol so eval($label) actually calls the variable associated
    # with label.
    eval(:($dataframe[$label] = eval($label)))
end

# Converts between the Julia color type and the traditional 3 tuple array
function colorArray2FloatArray{T<:Number}(colorArray::Array{ColorTypes.RGB{T},1})
    arr = zeros(size(colorArray, 1), 3)
    for (idx, color) in enumerate(colorArray)
        arr[idx, 1] = color.r
        arr[idx, 2] = color.g
        arr[idx, 3] = color.b
    end
    return arr
end

# Plots all the fitting traces from the s21FitJLD export file
function plotAllTracesJLD(jldFile)
    traces = load(jldFile)["traces"]
    results = load(jldFile)["results"]
    q = abs(Array(results[:Qi]))*1e6
    f0 = abs(Array(results[:f0]))
    for i in 1:maximum(size(traces))
        figure()
        plot(Array(traces[:freq][i])/1e9, complex2dB(Array(traces[:s21][1i])),"ro")
        # plot(Array(traces[:freq][i]), complex2dB(Array(traces[:s21smooth][1i])), "bo")
        plot(Array(traces[:freq][i])/1e9, complex2dB(Array(traces[:s21fit][1i])), "g", lw = 3)
        xlabel("Frequency [GHz]")
        ylabel("Magnitude [dB]")
        title(string("f0 = ", f0[i], "\nQi = ", q[i] ))
    end
end

# Plots all the fitting traces from the s21FitJLD export file
function mergeJLDFileResults(jldFiles; saveName = "")
    if saveName == ""
        saveName = "jldMergedResults.jld"
    end

    fileKeys = []
    for (idx, file) in enumerate(jldFiles)
        push!(fileKeys,keys(load(file)))
        for key in fileKeys[idx]
            cmd = string(key, "_", idx, " = load(\"", file, "\")[\"", key, "\"]")
            println(cmd)
            eval(parse(cmd))
            if idx > 1
                cmd = string(key, "_", idx, " = vcat(", key, "_", idx, ",", key, "_", idx - 1,")")
                println(cmd)
                eval(parse(cmd))
            end
        end
    end
    global saveArray # Julia made me do this. eval doesn't know about variables outside the for loop!
    saveArray = []
    for key in fileKeys[1]
        cmd = string("push!(saveArray,", "\"", key ,"\", " , key, "_", length(jldFiles) , ")")
        println(cmd)
        eval(parse(cmd))
    end
    save(saveName, saveArray...)
end

# Plots all resonators except ones specified by exlude (indexed by order of frequency)
function plotAllResultsPowerSweepJLD(jldFile, resonators; plotType = "photon", exclude = [])
    results = load(jldFile)["results"]
    fits = load(jldFile)["fits"]
    sort!(results, cols = (:power, :f0),rev = (true, false))
    resultsGroup = results
    numPowers = convert(Int64,maximum(size(resultsGroup))/resonators)
    f0Group = results[:f0][1:resonators]
    for i in 1:numPowers
        for j in 1:resonators
            resultsGroup[:f0][resonators*(i-1)+j] = f0Group[j]
        end
    end
    sort!(resultsGroup, cols = (:f0, :power),rev = (false, true))

    cs = colorArray2FloatArray(colormap("RdBu", size(results,1)))
    f = resultsGroup[:f0]
    p = resultsGroup[:power]
    if plotType == "photon"
        ph = zeros(length(f))
        temp = DataFrame()
        for i in 1:length(f)
            fit = fits[i]
            # println(string(i, " - ", vnaPower2Photon(f[i], p[i], fit.minimum)))
            ph[i] = vnaPower2Photon(f[i], p[i], fit.minimum)
            temp = vcat(temp, join(resultsGroup[i,:], DataFrame(photons = ph[i]) , kind = :cross))
        end
        resultsGroup = temp
        sort!(resultsGroup, cols = (:f0, :photons),rev = (false, true))

        f = resultsGroup[:f0]
        ph = resultsGroup[:photons]

        for i in 1:resonators
            if !(i in exclude)
                plot(p[(i-1)*numPowers + 1:i*numPowers], log10(ph[(i-1)*numPowers + 1:i*numPowers]), "o-")
                xlabel("VNA Power [dBm]")
                ylabel("Photons <n>")
            end
        end

    else

    end
    q = abs(Array(resultsGroup[:Qi]))*1e6
    figure()
    subplot(111, projection="3d")
    for i in 1:resonators
        if !(i in exclude)
            # plot(f[(i-1)*numPowers + 1:i*numPowers]/1e9, p[(i-1)*numPowers + 1:i*numPowers] ,q[(i-1)*numPowers + 1:i*numPowers]/1e6, "o-", label = string(f0Group[i]))
            if plotType == "photon"
                plot(f[(i-1)*numPowers + 1:i*numPowers]/1e9, log10(ph[(i-1)*numPowers + 1:i*numPowers]) ,q[(i-1)*numPowers + 1:i*numPowers]/1e6, "o-", label = string(f0Group[i]))
                ylabel("log10(Photon Number <n>)")
            else
                plot(f[(i-1)*numPowers + 1:i*numPowers]/1e9, p[(i-1)*numPowers + 1:i*numPowers] ,q[(i-1)*numPowers + 1:i*numPowers]/1e6, "o-", label = string(f0Group[i]))
                ylabel("Power [dBm]")
            end
        end
    end
    # scatter(f/1e9, p , zs = q/1e6, c = cs)
    xlabel("Frequency [GHz]")
    # ylabel("Power [dBm]")
    # ylabel("log10(Photon Number <n>)")
    zlabel("Internal Quality Factor [Mil]")

    figure()
    for i in 1:resonators
        if !(i in exclude)
            # semilogy(p[(i-1)*numPowers + 1:i*numPowers] ,q[(i-1)*numPowers + 1:i*numPowers], "o-", label = string(f0Group[i]))
            if plotType == "photon"
                loglog(ph[(i-1)*numPowers + 1:i*numPowers] ,q[(i-1)*numPowers + 1:i*numPowers], "o-", label = string(f0Group[i]))
                xlabel("Photon Number <n>")
            else
                semilogy(p[(i-1)*numPowers + 1:i*numPowers] ,q[(i-1)*numPowers + 1:i*numPowers], "o-", label = string(f0Group[i]))
                xlabel("Power [dBm]")
            end
        end
    end
    # xlabel("Power [dBm]")
    # xlabel("Photon Number <n>")
    ylabel("Internal Quality Factor")
    # legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=2, mode="expand", borderaxespad=0.)
    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    figure()
    # plot(f,q/1e6, "o")
    for i in 1:resonators
        if !(i in exclude)
            semilogy(f[(i-1)*numPowers + 1:i*numPowers]/1e9 ,q[(i-1)*numPowers + 1:i*numPowers], "o-", label = string(f0Group[i]))
        end
    end
    xlabel("Frequency [GHz]")
    ylabel("Internal Quality Factor")
end


# Plots all resonators except ones specified by exlude (indexed by order of frequency)
function plotAllResultsTimeSweepJLD(jldFile)
    results = load(jldFile)["results"]
    sort!(results, cols = (:timeHours),rev = (false))
    hours = results[:timeHours]
    Qi = results[:Qi]
    f0 = results[:f0]

    figure()
    plot(hours, Qi, "o-")

    xlabel("Time [Hours]")
    ylabel("Internal Quality Factor [Mil]")

    figure()
    plot(hours, f0, "o-")

    xlabel("Time [Hours]")
    ylabel("Resonance Frequency")
    # legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #            ncol=2, mode="expand", borderaxespad=0.)
end

# reads and fits resonators in the JLD format
function s21FitPowerSweepJLD(file; dataKey = "data", powerKey = "dbm", ifbKey = "ifs", attenVNA = 0,  saveResults = false, saveName = "", iterPerCode = 5000)
    saveName == "" ? saveName = string("fitResults_", basename(file)) : nothing
    res = load(file) # Z: is mapped to nyx
    powers = float(res[powerKey]) - attenVNA
    ifBandwidths = float(res[ifbKey])
    df = DataFrame()
    # powers = linspace(0, -80, 9)
    # codeList = ["0 1", "11 1"] # ,
    codeList = ["0 1"] #
    println(string(0, " of ", size(res[dataKey], 1), " processed"))
    resultIdx = []
    traceIdx = []
    fitIdx = 0
    for i in 1:size(res[dataKey], 1)
        for (idx, power) in enumerate(powers)
            dataset = res[dataKey][i,idx]
            freq = Array(dataset[:f])
            s21 = Array(dataset[:S21])
            fitResults, resultIdx, traceIdx, fitIdx = s21FitAuto(freq, s21, power = power, ifBandwidth = ifBandwidths[idx], codes = codeList, iterPerCode = iterPerCode)
            results = bestFit(fitResults)
            df = vcat(df, results)
        end
        println(string(i, " of ", size(res[dataKey], 1), " processed"))
    end
    sort!(df, cols = (:power, :f0),rev = (true, false))
    saveResults? save(saveName, "results", df[resultIdx], "traces", df[traceIdx], "fits", df[fitIdx]) : nothing
    return df[resultIdx]
end

# reads and fits resonators in the JLD format
function s21FitJLD(file; dataKey = "data", saveResults = false, saveName = "", iterPerCode = 2000)
    saveName == "" ? saveName = string("fitResults_", basename(file)) : nothing
    res = load(file) # Z: is mapped to nyx
    df = DataFrame()
    # powers = linspace(0, -80, 9)
    codeList = ["0 1", "11 1"] # ,
    # codeList = ["0 1"] #
    # println(string(0, " of ", size(res[dataKey], 1), " processed"))
    resultIdx = []
    traceIdx = []
    dataset = res[dataKey]
    freq = Array(dataset[:f])
    s21 = Array(dataset[:S21])
    fitResults, resultIdx, traceIdx, fitIdx = s21FitAuto(freq, s21, codes = codeList, iterPerCode = iterPerCode)
    results = bestFit(fitResults)
    df = vcat(df, results)
    sort!(df, cols = (:power, :f0),rev = (true, false))
    saveResults? save(saveName, "results", df[resultIdx], "traces", df[traceIdx], "fits", df[fitIdx]) : nothing
    return df[resultIdx]
end


# reads and fits resonators in the JLD format
function s21FitTimeSweepJLD(file; dataKey = "data", timeKey = "starttime",  saveResults = false, saveName = "", iterPerCode = 2000)
    saveName == "" ? saveName = string("fitResults_", basename(file)) : nothing
    res = load(file) # Z: is mapped to nyx
    df = DataFrame()
    timeHours = (res[timeKey] - res[timeKey][1])/3600
    codeList = ["0 1"] #, "11 1"] # ,
    println(string(0, " of ", size(res[dataKey], 1), " processed"))
    resultIdx = []
    traceIdx = []
    for (idx, hour) in enumerate(timeHours)
        dataset = res[dataKey][idx]
        freq = Array(dataset[:f])
        s21 = Array(dataset[:S21])
        fitResults, resultIdx, traceIdx, fitIdx = s21FitAuto(freq, s21, power = -40, ifBandwidth = 10, codes = codeList, iterPerCode = iterPerCode)
        results = bestFit(fitResults)
        results = join(results, DataFrame(timeHours = hour) , kind = :cross)
        df = vcat(df, results)
        println(string(idx, " of ", size(res[dataKey], 1), " processed"))
    end
    resultIdx = vcat(resultIdx, length(df))
    sort!(df, cols = (:timeHours),rev = (false))
    saveResults? save(saveName, "results", df[resultIdx], "traces", df[traceIdx], "fits", df[fitIdx]) : nothing
    return df[resultIdx]
end

macro assign(syms...)
    Expr(:block, [Expr(:(=), Expr(:ref, :fitTrails, QuoteNode(s)), s) for s in syms]...)
end

# Tries different fits given a code defining parameters for fit
function s21FitAuto(freq, s21; power = NA, ifBandwidth = NA, codes = ["0 1"], iterPerCode = 1000)
    """
    Automatically runs a set of ftting parameters and outputs a table with fit scoresnh n
    code[1] corresponds to the smoothingFactor parameter in the s21Fit function
    code[2] corresponds to the fitFreqOffset parameter in the s21Fit function
    """
    global fitTrails
    fitTrails = DataFrame()
    labels = [:code,:f0, :power, :ifBandwidth, :Qi, :dBOffset, :strayInductance, :dipDepth, :phaseOnResonance,
              :score, :grade, :results, :freq, :s21, :s21fit, :s21smooth] #:freqOffset,
    types = [AbstractString, Float64,  Float64,  Float64,  Float64,  Float64,  Float64,  Float64, Float64,
     Float64, AbstractString, Any, Array{Float64, 1}, Array{Complex{Float64},1},
     Array{Complex{Float64},1}, Array{Complex{Float64},1}]
    for (idx, label) in enumerate(labels)
        t = types[idx]
        fitTrails[label] = Array{t}(length(codes))
    end
    # warn(fitTrails)
    for i in 1:length(codes)
        # warn(codes[i])
        sf, ffo = tryFit(codes[i])
        results, score = s21Fit(freq, s21, iter = iterPerCode)
        grade = fitGrade(score, maximum(size(freq)))
        Qi = results.minimum[1]
        dBOffset = results.minimum[2]
        strayInductance = results.minimum[3]
        dipDepth = results.minimum[4]
        f0 = results.minimum[5]
        phaseOnResonance = results.minimum[6]
        code = codes[i]
        s21fit = s21FitFunc(freq, results.minimum)
        s21smooth = s21fit
        ## Gives error that whatever variable I'm using is undefined e.g. f0
        ## Probably a scope problem
        # for label in labels
        #     cmd = string("fitTrails[:", label,"] = ", label)
        #     eval(parse(cmd))
        #     # fitTrails[label] = eval(label)
        # end
        # @assign code, f0, power, ifBandwidth, Qi, dBOffset, strayInductance, dipDepth, phaseOnResonance, score, grade, results, freq, s21, s21fit, s21smooth
        # # Instead, I have to do it the stupid way
        fitTrails[symbol("code")][i] = code
        fitTrails[symbol("f0")][i] = f0
        fitTrails[symbol("power")][i] = power
        fitTrails[symbol("ifBandwidth")][i] = ifBandwidth
        fitTrails[symbol("Qi")][i] = Qi
        fitTrails[symbol("dBOffset")][i] = dBOffset
        fitTrails[symbol("strayInductance")][i] = strayInductance
        fitTrails[symbol("dipDepth")][i] = dipDepth
        fitTrails[symbol("phaseOnResonance")][i] = phaseOnResonance
        fitTrails[symbol("score")][i] = score
        fitTrails[symbol("grade")][i] = grade
        fitTrails[symbol("results")][i] = results
        fitTrails[symbol("freq")][i]  = freq
        fitTrails[symbol("s21")][i]  = s21
        fitTrails[symbol("s21fit")][i]  = s21fit
        fitTrails[symbol("s21smooth")][i]  = s21smooth
    end
    resultIdx = 2:11
    fitIdx = 12
    traceIdx = 13:16
    return fitTrails, resultIdx, traceIdx, fitIdx
end

"""
    s21Fit(freq, s21; pGuess = [0.5], smoothingFactor = 0, iter = 2000, funcTol = 1e-10 )

This function will fit a resonator's S21 transmission response to a fitting
model described in Megrant et. al. ArXiv 1201.3384 (2012). Minimizes
the square residuals using Nelder Mead optimizing provided in `Optim.jl`.\n
__Arg(s)__:\n
`freq` - real array of frequencies\n
`s21` - complex array of the S21 parameters\n
(`pGuess`) - array of length 6 for the initial fit parameter guess. The components
for fitting are as follows:\n
+ Qi/1e6\n
+ dBOffsetFromZero\n
+ strayInductance(nH)\n
+ dipDepth (magnitude)\n
+ f0(GHz)\n
+ phaseOnResonance\n
(`smoothingFactor`) - smooths the s21 input using Savitzky–Golay filter with a window size equal to the smoothing factor\n
(`iter`) - maximum number of iterations for convergence before optimization gives upgrade_scalar\n
(`funcTol`) - tolerance for the change in the fit function value between successive iterations\n
__Output(s)__:\n
`s21fit` - the fit s21 function given from the same freqeuency input and optimized parameters\n
`s21smoothed` - smoothed s21 function (if not, then original)\n
`fitResults` - Optim.jl fitting results where calling fitResults.minimum will give optimized paramters\n
`fitScore` - value of fitting function at optimized parameters\n
__Example__:\n
```
s21fit, s21smoothed, fitResults, fitScore = s21Fit(freq, s21)
```
"""
function s21Fit(freq, s21; pGuess = [0.5], iter = 2000, funcTol = 1e-10 ) # smoothingFactor = 0
    # # Save a copy of orginal s21 before smoothing
    # s21org = s21
    # # make smoothingFactor an odd integer
    # smoothingFactor = round(Int, smoothingFactor)
    # isodd(smoothingFactor) ? nothing : (smoothingFactor == 0? nothing : smoothingFactor += 1)
    # # smooth s21 if non-zero with polynomial of degree one (linear smoothing)
    # smoothingFactor > 0 ? s21 = savitsky_golay(s21, smoothingFactor, 1) : nothing
    """
     The six parameters are as follows:
    1. Qi/1e6
    2. dBOffsetFromZero
    3. strayInductance(nH)
    4. dipDepth (magnitude)
    5. f0(GHz)
    6. phaseOnResonance
    """
    pSize = 6
    # Make the parameter smaller if too big
    while maximum(size(pGuess, 1)) > pSize
        pop!(pGuess);
    end
    # Make the parameter bigger if too small
    while maximum(size(pGuess, 1)) < pSize
        push!(pGuess, 0)
    end
    # Separate into magnitude and phase
    s21dB = ResonatorFit.complex2dB(s21)
    s21Ph = ResonatorFit.complex2Phase(s21)
    # Get top value of flat area
    s1dBMax = s21dB[indmax(s21dB)]
    # Get resonance value
    s1dBMin = s21dB[indmin(s21dB)]
    # Index of resonance frequency
    f0Idx = indmin(s21dB)
    # Phase at reonance frequency
    # minPhase = s21Ph[f0Idx]  = p[6]
    # Resonance frequency
    pGuess[5] = freq[f0Idx] + pGuess[5]
    # Maximum value of trace to zero
    pGuess[2] = s1dBMax + pGuess[2]
    # Difference in magnitude from top of trace to the bottom of the dip
    pGuess[4] = dB2mag(pGuess[2] - s1dBMin) + pGuess[4]
    # The residual of the dB magnitude for the spectrum
    residualdB(freq, p) = complex2dB(s21) - complex2dB(s21FitFunc(freq,p))
    # The residual of the dB phase for the spectrum
    residualPh(freq, p) = complex2Phase(s21) - complex2Phase(s21FitFunc(freq,p))
    # # Errors created from smoothing the orginal data will count against fit
    # # Magnitude
    # residualSmoothdB = complex2dB(s21org) - complex2dB(s21)
    # # Phase
    # residualSmoothPh = complex2Phase(s21org) - complex2Phase(s21)
    # The magnitude of the combined residuals
    residualMag(freq, p) = sqrt(residualdB(freq, p).^2 + residualPh(freq, p).^2)
         #+ residualSmoothdB.^2 + residualSmoothPh.^2)
    # Least squares
    rsquared(p) = sum(abs(residualMag(freq, p)).^2)
    # Optimization
    fitResults = optimize(rsquared, pGuess, ftol = funcTol, iterations = iter)
    fitScore = rsquared(fitResults.minimum) # lower the better
    # minimumOutputKey = [Qi/1e6, dBOffsetFromZero, strayInductance(nH), dipDepthMag, f0, phaseOnResonance]
    return fitResults, fitScore #s21FitFunc(freq, fitResults.minimum),  s21,
end

"""
    s21FitFunc(freq, p)

This function defines the resonance behavior in a parameterized functional form.
The theory for this model can be found in section 2.3.3 and 2.4 of Ben Mazin's thesis\n
__Arg(s)__:\n
`freq` - frequency value or array in Hz\n
`p` -  array of length 6 for fit parameters. The components for fitting are as follows:\n
+ Qi/1e6\n
+ dBOffsetFromZero\n
+ strayInductance(nH)\n
+ dipDepth (magnitude)\n
+ f0(GHz)\n
+ phaseOnResonance\n
`Qi`, `strayInductance`, and `f0` were choosen to be represented in Millions, nanoHenries, and GHz
respectively to due to their scales (1e6, 1e-9, 1e9) relative to the other parameters.
__Output(s)__:\n
`s21val` - the complex value of s21 evaluated at `freq` given fit parameters\n
"""
function s21FitFunc(freq, p)
    # Define the fitting function in a parameterized functional form
    # The theory for this model can be found in section 2.3.3 and 2.4 of Ben Mazin's thesis

    # Difference in s21 between the dip and top of resonance curved,
    # normalized such that top of resonance curve is at 0 dB in units of
    # magnitude, i.e. NOT dB
    s21AdjMin(p) = p[4]
    # Offset between top of resonance curve and 0 dB
    totaloffset(p) = dB2mag(p[2])
    # Add some parasitic series indutance to the transmission line
    Z(freq, p) = 50 + 2im*pi*freq*p[3]/1E9
    # Define a scale factor based on this impedance scaling change
    scalefactor(freq, p) = 50./Z(freq, p)
    # Calcuate the coupling capacitance scaled by the new impedance
    Qc(freq, p) = scalefactor(freq, p)*s21AdjMin(p)*p[1]*1e6/(1 - s21AdjMin(p))
    # Calculate the loaded Q
    Q0(freq, p) = (p[1]*(1e6).*Qc(freq, p))./(p[1]*1e6+Qc(freq, p))
    # Some variable I arbitrarily defined
    k(freq, p) = 2im*Q0(freq, p).*(freq - p[5])/p[5]
    # The fitting function in all it's glory
    return (s21AdjMin(p) + k(freq, p)).*exp(1im*p[6])./(1 + k(freq, p))*totaloffset(p)
end

"""
    getQFromFitParams(f0, p)

This function defines the resonance behavior in a parameterized functional form.
The theory for this model can be found in section 2.3.3 and 2.4 of Ben Mazin's thesis.\n
__Arg(s)__:\n
`f0` - resonance frequency\n
`p` -  array of length 6 for fit parameters. The components for fitting are as follows:\n
+ Qi/1e6\n
+ dBOffsetFromZero\n
+ strayInductance(nH)\n
+ dipDepth (magnitude)\n
+ f0(GHz)\n
+ phaseOnResonance\n
__Output(s)__:\n
`Qc` - coupling Q\n
`Q0` - loaded Q\n
`Qi` - intrinsic Q/1e6, equal to `p[1]`
"""
function getQFromFitParams(f0, p)
    s21AdjMin(p) = p[4]
    # Offset between top of resonance curve and 0 dB
    totaloffset(p) = dB2mag(p[2])
    # Add some parasitic series indutance to the transmission line
    Z(f0, p) = 50 + 2im*pi*f0*p[3]/1e9
    # Define a scale factor based on this impedance scaling change
    scalefactor(f0, p) = 50./Z(f0, p)
    # Calcuate the coupling capacitance scaled by the new impedance
    Qc = abs(scalefactor(f0, p)*s21AdjMin(p)*p[1]*1e6/(1 - s21AdjMin(p)))
    # Calculate the loaded Q
    Q0 = (p[1]*(1e6).*Qc)./(p[1]*1e6 + Qc)
    Qi = p[1]*(1e6)
    return Qc, Q0, Qi
end

"""
    vnaPower2Photon(f0, vnaPower, p)


__Arg(s)__:\n
__Output(s)__:\n
"""
function vnaPower2Photon(f0, vnaPower, p)
    # # Attenuators at the output of VNA
    # vnaAtten = 0
    # # Cryo attenutators
    # cryoAtten = -60
    PACdB = powerAtChip(f0, vnaPower)
    # # Cable Attenutation
    # cableAtt = cableAtten(f0, cryoAtten = cryoAtten)
    # # Total power at chip in dB
    # PACdB = cableAtt + vnaPower + cryoAtten + vnaAtten
    Qc, Q0, Qi = getQFromFitParams(f0, p)
    energyInternal = resonatorEnergy(10, 5, f0, PACdB, Qi, Qc)
    plankConstant = 6.63e-34
    numberPhotons = energyInternal / (plankConstant * f0)
    return numberPhotons
end

"""
    vnaPower2Photon(f0, vnaPower, Qi, Qc)


__Arg(s)__:\n
__Output(s)__:\n
"""
function vnaPower2Photon(f0, vnaPower, Qi, Qc)
    # # Attenuators at the output of VNA
    # vnaAtten = 0
    # # Cryo attenutators
    # cryoAtten = -60
    PACdB = powerAtChip(f0, vnaPower)
    # # Cable Attenutation
    # cableAtt = cableAtten(f0, cryoAtten = cryoAtten)
    # # Total power at chip in dB
    # PACdB = cableAtt + vnaPower + cryoAtten + vnaAtten
    plankConstant = 6.63e-34
    energyInternal = resonatorEnergy(10, 5, f0, PACdB, Qi, Qc)
    numberPhotons = energyInternal / (plankConstant * f0)
    return numberPhotons
end

"""
    powerAtChip(f0, vnaPower; vnaAtten = 0,  cryoAtten = -60, hemt = 40, rtAmp = 26.5)


__Arg(s)__:\n
`f0` - resonance frequency\n
`vnaPower` - VNA output power in dBm\n
(`vnaAtten`) - attenuation at the output of the VNA\n
(`cryoAtten`) - cryogenic attenutation placed in wiring of fridge\n
__Output(s)__:\n
"""
function powerAtChip(f0, vnaPower; vnaAtten = 0,  cryoAtten = -60, hemt = 40, rtAmp = 26.5)
    # Attenuators at the output of VNA
    vnaAtten = 0
    # Cable Attenutation
    cableAtt = cableAtten(f0)
    # Total power at chip in dB
    PACdB = cableAtt + vnaPower + cryoAtten + vnaAtten
end
"""
    resonatorEnergy(widthCPW, gapCPW, f0, powerAtFeedline, Ql, Qc; dielSubstrate = 11.9)

__Arg(s)__:\n
`widthCPW` - center trace width for the resonator CPW geometry\n
`gapCPW` - center trace gap for the resonator CPW geometry\n
`f0` - resonance frequency of lambda/4 resonator\n
`powerAtFeedline` - effective input power at the feedline coupled to the resonator\n
`Qi` - resonator internal quality factor\n
`Qc` - resonator coupling quality factor\n
(`dielSubstrate`) - dieletric constant of device substrate\n
__Output(s)__:\n
`energyInternal` - internal energy (Joules) stored inside resonator at `f0`
"""
function resonatorEnergy(widthCPW, gapCPW, f0, powerAtFeedline, Qi, Qc; dielSubstrate = 11.9)
    Ql = (Qi.^(-1) + Qc.^(-1)).^(-1)
    dielConstApprox = (1 + dielSubstrate)/2 # Good approx for thick substrates
    k1 = widthCPW/(widthCPW+2*gapCPW)
    k2 = sqrt(1 - k1^2)
    vacPerm = 8.854187e-12
    capacitancePerUnitLength = 4*vacPerm*ceifk(k1, 24)/ceifk(k2, 24) #Rami Thesis Eq 3.3
    length = (3e8/sqrt(dielConstApprox)/f0/4)
    powerReadout = 10.^(powerAtFeedline/10)
    powerInternal = (2/pi)*(Ql.^2 / Qc)*powerReadout # Rami's Thesis eq 3.30
    powerInternaldBm = powerInternal/1000
    voltageResonator = 2*(50*powerInternaldBm).^.5
    energyInternal = .5*capacitancePerUnitLength*voltageResonator.^2*length
    return energyInternal
end

"""
    ceifk(x, n)

Complete elliptic integral of the first kind with an approximation of order 2(n -1).
Approximation to integral from: http://www.exstrom.com/math/elliptic/ellipint.html \n
__Arg(s)__:\n
`x` - elliptic integral argument, can be single value or array
`n` - order of approximation
__Output(s)__:\n
`k` - approximate value of the complete elliptic integral of the first kind at `x`
"""
function ceifk(x, n)
    s = 0
    for i in 0:n
        s += binomial(2*i, i).^2.*(x./4).^(2*i)
    end
    return length(x)>1 ? pi/2*sum(s,2) : pi/2*sum(s)
end

"""
    cableAtten(freq; hemt = 40, rtAmp = 26.5, cryoAtten = -60, measuredPts = [(3.5e9,-8.26),(8e9,-13.26)])

Linear frequency dependent approximation of the attenuation due to cable losses.\n
__Arg(s)__:\n
`freq` - frequency span for frequency dependent loss\n
`hemt` - gain in dB of the HEMT amplifier at 3K\n
`rtAmp` - gain in dB of the root temperature amplifier\n
`cryoAtten` - gain in dB of the cyro-attenuators (should be negative)\n
`measuredPts` - a two point frequency measurement to measure the slope of the attenuation in a two element tuple array format:\n
```
measuredPts = [(freqLow, s21dBValue), (freqHigh, s21dBValue)]
```
__Output(s)__:\n
"""
function cableAtten(freq; hemt = 40, rtAmp = 26.5, cryoAtten = -60, measuredPts = [(3.5e9,-8.26),(8e9,-13.26)])
    # Effect of amps and attens
    netGain = hemt + rtAmp + cryoAtten
    # Example: slope and offset obtained from:
    # -8.26 dB Transmission @ 3.5 GHz
    # -13.26 dB Transmission @ 8 GHz
    # cableAtten = slope*freq + offsetCorrection
    # VNA Measure = cableAtten + netGain
    # slope = - 5dB / 4.5 GHz
    slope = (measuredPts[2][2] - measuredPts[1][2])/(measuredPts[2][1] - measuredPts[1][1])
    offsetCorrection = measuredPts[1][2]- slope*measuredPts[1][1] - netGain
    cableOnlyAtten = slope*freq + offsetCorrection
    # close enough...
    oneWayAtten = cableOnlyAtten/2
end

"""
    fitGrade(fitScore, pts)

Give an arbitrary grade based on the fit score.\n
__Arg(s)__:\n
`fitScore` - square residual val from minimization function \n
`pts` -  length of fitting function, used so that the grade is normalized to the length of the fit\n
__Output(s)__:\n
`grade` - welcome back to school\n
"""
function fitGrade(fitScore, pts)
    if fitScore/pts < 0.0025
        "A"
    elseif fitScore/pts < 0.01
        "B"
    elseif fitScore/pts < 0.04
        "C"
    elseif fitScore/pts < 0.16
        "D"
    elseif fitScore/pts < 0.32
        "E"
    else
        "F"
    end
end

"""
    savitsky_golay(x::Array{Float64, 1}, windowSize::Int64, polyOrder::Int64; deriv::Int64=0)

Polynomial smoothing with the Savitsky Golay filters. The Savitzky-Golay is a
type of low-pass filter, particularly suited for smoothing noisy data.
The main idea behind this approach is to make for each point a least-square
fit with a polynomial of high order over a odd-sized window centered at
the point.\n
Orginal author BBN (https://github.com/BBN-Q/Qlab.jl). \n
Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf\n
Python Example: http://wiki.scipy.org/Cookbook/SavitzkyGolay\n
__Arg(s)__:\n
`x` - the values of the time history of the signal\n
`windowSize` -  size of smoothing window. Must be an odd integer number.\n
`polyOrder` - the order of the polynomial used in the filtering. Must be less then `window_size` - 1.\n
`deriv`- the order of the derivative to compute (default = 0 means only smoothing)\n
__Output(s)__:\n
`ys` - Array of the smoothed signal (or it's n-th derivative)\n
"""
function savitsky_golay(x::Array{Float64, 1}, windowSize::Int64, polyOrder::Int64; deriv::Int64=0)
#Some error checking
@assert isodd(windowSize) "Window size must be an odd integer."
@assert polyOrder < windowSize "Polynomial order must me less than window size."
halfWindow = round(Int,(windowSize-1)/2)
#Setup the S matrix of basis vectors.
S = zeros(windowSize, polyOrder+1)
for ct = 0:polyOrder
	S[:,ct+1] = collect(-halfWindow:halfWindow).^(ct)
end
#Compute the filter coefficients for all orders
#From the scipy code it seems pinv(S) and taking rows should be enough
G = S*pinv(S'*S)
#Slice out the derivative order we want
filterCoeffs = G[:,deriv+1] * factorial(deriv);
#Pad the signal with the endpoints and convolve with filter
paddedX = [x[1]*ones(halfWindow); x; x[end]*ones(halfWindow)]
y = conv(filterCoeffs[end:-1:1], paddedX)
#Return the valid midsection
return y[2*halfWindow+1:end-2*halfWindow]
end

"""
    savitsky_golay(x::Array{Complex{Float64},1}, windowSize::Int64, polyOrder::Int64)

Savitsky Golay filters for complex arrays.
__Arg(s)__:\n
`x` - the complex values of the time history of the signal\n
`windowSize` -  size of smoothing window. Must be an odd integer number.\n
`polyOrder` - the order of the polynomial used in the filtering. Must be less then `window_size` - 1.\n
__Output(s)__:\n
`ys` - complex array of the smoothed signal \n
"""
function savitsky_golay(x::Array{Complex{Float64},1}, windowSize::Int64, polyOrder::Int64)
    return savitsky_golay(real(x), windowSize, polyOrder) + 1im*savitsky_golay(imag(x), windowSize, polyOrder)
end


# Assign a varaible whos name is given by a string s to a value v
# function stringAsVarName(s::String,v::Any)
#          s=symbol(s)
#          @eval (($s) = ($v))
# end

# Automatically extrats variables from a dataframe. This function assigns it's
# column label as a variable to the label's column data
# Mar 24, 2016 THIS DOES NOT WORK
function assignVarFromName(dataframe::DataFrame, print = false)
    global foo = dataframe
    for name in names(foo)
        if print == true
            println(name)
        end
        eval(parse(string("global ", name)))
        eval(parse(string(name, " = foo[:",name,"]")))
        eval(parse(string("println(head(", name, "))")))
        # println(string("global ", name, " = foo[:",name,"]"))
    end

    foo = nothing
end

# Mar 24, 2016 THIS DOES NOT WORK
function dataframeReturnHelper(dataframe::DataFrame, returnType::DataType)
    global foo = dataframe
    returnArray = []
    for name in names(foo)
        println(name)
        println(returnType)
        # eval(parse(string("global ", name)))
        eval(parse(string(name, " = foo[:",name,"]")))
        # eval(parse(string("println(head(", name, "))")))
        # println(string("global ", name, " = foo[:",name,"]"))
        eval(parse(string(name, " = ", returnType, "(",name,")")))
        eval(parse(string("push!(returnArray,", name,")")))
    end
    foo = nothing
    println(returnArray)
    return returnArray
end


end # module
