module ResonatorFit

include("HelperFuncs.jl")
include("CircleFit.jl")

using DataArrays, DataFrames, Optim, JLD, PyPlot

export  s21Fit,
        s21FitAuto,
        # bestFit,
        s21FitPowerSweepJLD,
        s21FitTimeSweepJLD,
        plotAllTracesJLD,
        plotAllResultsPowerSweepJLD,
        plotAllResultsTimeSweepJLD,
        s21FitFunc,
        s21FitJLD,
        getQFromFitParams,
        fitGrade,
        cableAtten,
        vnaPower2Photon

# # decodes the fitCode as a smoothingFactor and a fitFreqOffset flag
# function tryFit(fitCode)
#     decoded = readdlm(IOBuffer(fitCode))
#     sf = decoded[1]
#     ffo = (decoded[2] != 0)
#     return sf, ffo
# end
# takes the best fit from the DataFrame output
# function bestFit(fitResults::DataFrame)
#     sort!(fitResults, cols = [:score])
#     best = fitResults[1,:]
#     # f0, Qi, score, grade = (float(fitResults[:f0][1]), float(fitResults[:Qi][1]), float(fitResults[:score][1]), fitResults[:grade][1])
#     return best
# end

# Plots all the fitting traces from the s21FitJLD export file
function plotAllTracesJLD(jldFile; plotType = "dB")
    traces = load(jldFile)["traces"]
    results = load(jldFile)["results"]
    q = abs(Array(results[:Qi]))
    f0 = abs(Array(results[:f0]))

    for i in 1:maximum(size(traces))
        figure()
        if plotType == "dB"
            plot(Array(traces[:freq][i])/1e9, complex2dB(Array(traces[:s21][i])),"ro")
            plot(Array(traces[:freq][i])/1e9, complex2dB(Array(traces[:s21fit][i])), "g", lw = 3)
            xlabel("Frequency [GHz]")
            ylabel("Magnitude [dB]")
        elseif plotType == "Phase"
            plot(Array(traces[:freq][i])/1e9, complex2Phase(Array(traces[:s21][i])),"ro")
            plot(Array(traces[:freq][i])/1e9, complex2Phase(Array(traces[:s21fit][i])), "g", lw = 3)
            xlabel("Frequency [GHz]")
            ylabel("Phase [rad]")
        elseif plotType == "ReIm"
            plot(real(Array(traces[:s21][i])), imag(Array(traces[:s21][i])),"ro")
            plot(real(Array(traces[:s21fit][i])), imag(Array(traces[:s21fit][i])), "g", lw = 3)
            xlabel("Real S21")
            ylabel("Imag S21")
        end
        title(string("f0 = ", f0[i], "\nQi = ", q[i] ))
    end
end

# Plots all resonators except ones specified by exlude (indexed by order of frequency)
function plotAllResultsPowerSweepJLD(jldFile; plotType = "photon", exclude = [])
    results = load(jldFile)["results"]
    sort!(results, cols = (:resonator, :power),rev = (false, true))
    resonators = results[:resonator]

    figure()
    subplot(111, projection="3d")
    plotType = "photon"
    for resonator in 1:length(resonators)
            res = extractResonator(results, resonator)
            yParam = res[:power]
            ylabel("Power [dBm]")
            if plotType == "photon"
                for (idx, power) in enumerate(yParam)
                    yParam[idx] = log10(vnaPower2Photon(res[:f0][idx], power, res[:Qi][idx], res[:Qc][idx]))
                end
                ylabel("log10(Photons)")
            end
            plot(Array(res[:f0]), Array(yParam), Array(res[:Qi]), "o-", label = string("Res#", resonator))
    end
    xlabel("Frequency [GHz]")
    zlabel("Internal Quality Factor")

    figure()
    for resonator in 1:length(resonators)
        res = extractResonator(results, resonator)
        yParam = res[:power]
        xlabel("Power [dBm]")
        if plotType == "photon"
            for (idx, power) in enumerate(res[:power])
                yParam[idx] = log10(vnaPower2Photon(res[:f0][idx], power, res[:Qi][idx], res[:Qc][idx]))
            end
            xlabel("log10(Photons)")
        end
        semilogy(Array(yParam), Array(res[:Qi]), "o-", label = string("Res#", resonator))
    end
    ylabel("Internal Quality Factor")

    figure()
    for resonator in 1:length(resonators)
        res = extractResonator(results, resonator)
        semilogy(Array(res[:f0]), Array(res[:Qi]), "o-", label = string("Res#", resonator))
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
function s21FitPowerSweepJLD(file; dataKey = "data", powerKey = "dbm", ifbKey = "ifs", attenVNA = 0,  saveResults = false, saveName = "", iterPerCode = 5000, doCircleFit = true)
    saveName == "" ? saveName = string("fitResults_", basename(file)) : nothing
    res = load(file) # Z: is mapped to nyx
    powers = float(res[powerKey]) - attenVNA
    ifBandwidths = float(res[ifbKey])
    df = DataFrame()
    # powers = linspace(0, -80, 9)
    # codeList = ["0 1", "11 1"] # ,
    # codeList = ["0 1"] #
    println(string(0, " of ", size(res[dataKey], 1), " processed"))
    resultIdx = []
    traceIdx = []
    fitIdx = []
    numRes = size(res[dataKey], 1)
    for i in 1:numRes
        for (idx, power) in enumerate(powers)
            dataset = res[dataKey][i,idx]
            freq = Array(dataset[:f])
            s21 = Array(dataset[:S21])
            fitResults, resultIdx, traceIdx, fitIdx = s21FitAuto(freq, s21, power = power, ifBandwidth = ifBandwidths[idx], iterPerCode = iterPerCode, doCircleFit = doCircleFit)
            # results = bestFit(fitResults)
            resLabel = DataFrame(resonator = i)
            results = hcat(resLabel, fitResults)
            df = vcat(df, results)
        end
        println(string(i, " of ", numRes, " processed"))
    end
    sort!(df, cols = (:power, :resonator),rev = (true, false))
    saveResults? save(saveName, "results", df[collect([1; resultIdx + 1])], "traces", df[traceIdx+1], "fits", df[fitIdx+1]) : nothing
    return df[collect([1; resultIdx + 1])]
end

# reads and fits resonators in the JLD format
function s21FitJLD(file; dataKey = "data", saveResults = false, saveName = "", iterPerCode = 2000, doCircleFit = true)
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
    fitResults, resultIdx, traceIdx, fitIdx = s21FitAuto(freq, s21, iterPerCode = iterPerCode, doCircleFit = doCircleFit)
    results = bestFit(fitResults)
    df = vcat(df, results)
    sort!(df, cols = (:power, :f0),rev = (true, false))
    saveResults? save(saveName, "results", df[resultIdx], "traces", df[traceIdx], "fits", df[fitIdx]) : nothing
    return df[resultIdx]
end


# reads and fits resonators in the JLD format
function s21FitTimeSweepJLD(file; dataKey = "data", timeKey = "starttime",  saveResults = false, saveName = "", iterPerCode = 2000, doCircleFit = true)
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
        fitResults, resultIdx, traceIdx, fitIdx = s21FitAuto(freq, s21, power = -40, ifBandwidth = 10, iterPerCode = iterPerCode, doCircleFit = doCircleFit)
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
    Expr(:block, [Expr(:(=), Expr(:ref, :fitTrials, QuoteNode(s)), s) for s in syms]...)
end

# Tries different fits given a code defining parameters for fit
function s21FitAuto(freq, s21; power = NA, ifBandwidth = NA, iterPerCode = 1000, doCircleFit = true)
    """
    Automatically extracts fit results into a table with fit scores
    """
    global fitTrials
    fit = DataFrame()
    labels = [:f0, :power, :ifBandwidth, :Qi, :Qc, :Q0, :dBOffset, :strayInductance, :dipDepth,
    :phaseOnResonance, :score, :grade, :optimResults, :guess, :freq, :s21, :s21fit]
    types = [Float64,  Float64,  Float64,  Float64,  Float64, Float64,  Float64,  Float64,  Float64,
    Float64, Float64, AbstractString, Optim.MultivariateOptimizationResults{Float64,1}, Array{Float64, 1},
    Array{Float64, 1}, Array{Complex{Float64},1}, Array{Complex{Float64},1}]
    for (idx, label) in enumerate(labels)
        t = types[idx]
        fit[label] = Array{t}(1)
    end
    optimResults, score, guess = s21Fit(freq, s21, iter = iterPerCode, doCircleFit = doCircleFit)
    grade = fitGrade(score, maximum(size(freq)))

    Qi = optimResults.minimum[1]
    dBOffset = optimResults.minimum[2]
    strayInductance = optimResults.minimum[3]
    dipDepth = optimResults.minimum[4]
    f0 = optimResults.minimum[5]
    phaseOnResonance = optimResults.minimum[6]
    Qc, Q0, Qi = getQFromFitParams(f0, optimResults.minimum)
    Qc = Qc*1e6
    Q0 = Q0*1e6
    Qi = Qi*1e6
    # code = codes[i]
    s21fit = s21FitFunc(freq, optimResults.minimum)

    ## Gives error that whatever variable I'm using is undefined e.g. f0
    ## Probably a scope problem
    # for label in labels
    #     cmd = string("fitTrials[:", label,"] = ", label)
    #     eval(parse(cmd))
    #     # fitTrials[label] = eval(label)
    # end
    # @assign code, f0, power, ifBandwidth, Qi, dBOffset, strayInductance, dipDepth, phaseOnResonance, score, grade, results, freq, s21, s21fit, s21smooth
    # # Instead, I have to do it the long way
    i = 1
    # Fit Results
    fit[symbol("f0")][i] = f0
    fit[symbol("power")][i] = power
    fit[symbol("ifBandwidth")][i] = ifBandwidth
    fit[symbol("Qi")][i] = Qi
    fit[symbol("Qc")][i] = Qc
    fit[symbol("Q0")][i] = Q0
    fit[symbol("dBOffset")][i]= dBOffset
    fit[symbol("strayInductance")][i] = strayInductance
    fit[symbol("dipDepth")][i] = dipDepth
    fit[symbol("phaseOnResonance")][i]= phaseOnResonance
    fit[symbol("score")][i] = score
    fit[symbol("grade")][i]= grade
    # Optim results output
    fit[symbol("optimResults")][i] = optimResults
    fit[symbol("guess")][i] = guess
    # Trace data
    fit[symbol("freq")][i] = freq
    fit[symbol("s21")][i] = s21
    fit[symbol("s21fit")][i] = s21fit

    resultIdx = 1:12
    optimIdx = 13:14
    traceIdx = 15:17
    return fit, resultIdx, traceIdx, optimIdx
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
function s21Fit(freq, s21; pGuess = [0.5], doCircleFit = true, iter = 4000, funcTol = 1e-12, weight = false) # smoothingFactor = 0
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

    if doCircleFit
        pGuess[1] = circleFit(freq, s21, plotFit = false)/1e6
    end

    # Make the parameter smaller if too big
    while maximum(size(pGuess, 1)) > pSize
        pop!(pGuess);
    end
    # Make the parameter bigger if too small
    while maximum(size(pGuess, 1)) < pSize
        push!(pGuess, 0)
    end
    # Separate into magnitude and phase
    s21dB = complex2dB(s21)
    s21Ph = complex2Phase(s21)
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
    # # Residuals for the real and imaginary parts
    # residualRe(freq, p) = real(s21) - real(s21FitFunc(freq,p))
    # residualIm(freq, p) = imag(s21) - imag(s21FitFunc(freq,p))

    # # The residual of the dB magnitude for the spectrum
    # residualdB(freq, p) = complex2dB(s21) - complex2dB(s21FitFunc(freq,p))
    # # The residual of the dB phase for the spectrum
    # residualPh(freq, p) = complex2Phase(s21) - complex2Phase(s21FitFunc(freq,p))


    # # The residual of the dB magnitude for the spectrum
    # residualdB(freq, p) = complex2dB(s21) - complex2dB(s21FitFunc2(freq,p)
    # # The residual of the dB phase for the spectrum
    # residualPh(freq, p) = complex2Phase(s21) - complex2Phase(s21FitFunc2(freq,p))

    # # Errors created from smoothing the orginal data will count against fit
    # # Magnitude
    # residualSmoothdB = complex2dB(s21org) - complex2dB(s21)
    # # Phase
    # residualSmoothPh = complex2Phase(s21org) - complex2Phase(s21)
    # # The magnitude of the combined residuals
    # residualMag(freq, p) = sqrt((complex2dB(s21) - complex2dB(s21FitFunc2(freq,p))).^2 + (complex2Phase(s21) - complex2Phase(s21FitFunc2(freq,p))).^2)
    #+ residualSmoothdB.^2 + residualSmoothPh.^2)

    # # The magnitude of the combined residuals
    # residualMag(freq, p) = sqrt(residualRe(freq, p).^2 + residualIm(freq, p).^2)
    if weight == true
        weights = -(s21dB - s1dBMax)/(s1dBMax - s1dBMin) + 0.5
        figure()
        plot(weights)
        title("fitting weight")
        rsquared(p) = sum((sqrt(weights.*(s21dB - complex2dB(s21FitFunc2(freq,p))).^2 + (s21Ph - complex2Phase(s21FitFunc2(freq,p))).^2)).^2)
    else
        rsquared(p) = sum((sqrt((s21dB - complex2dB(s21FitFunc2(freq,p))).^2 + (s21Ph - complex2Phase(s21FitFunc2(freq,p))).^2)).^2)

    end
    # figure()
    # plot(weights)
    # Least squares

    # Optimization
    fitResults = optimize(rsquared, pGuess, ftol = funcTol, iterations = iter)
    fitScore = rsquared(fitResults.minimum) # lower the better
    # minimumOutputKey = [Qi/1e6, dBOffsetFromZero, strayInductance(nH), dipDepthMag, f0, phaseOnResonance]
    return fitResults, fitScore, pGuess #s21FitFunc(freq, fitResults.minimum),  s21,
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


function s21FitFunc2(freq, p)
    # Define the fitting function in a parameterized functional form
    # The theory for this model can be found in section 2.3.3 and 2.4 of Ben Mazin's thesis

    # Difference in s21 between the dip and top of resonance curved,
    # normalized such that top of resonance curve is at 0 dB in units of
    # magnitude, i.e. NOT dB
    s21AdjMin = p[4]
    # Offset between top of resonance curve and 0 dB
    totaloffset = dB2mag(p[2])
    # Add some parasitic series indutance to the transmission line
    Z = 50 + 2im*pi*freq*p[3]/1E9
    # Define a scale factor based on this impedance scaling change
    scalefactor = 50./Z
    # Calcuate the coupling capacitance scaled by the new impedance
    Qc = scalefactor*s21AdjMin*p[1]*1e6/(1 - s21AdjMin)
    # Calculate the loaded Q
    Q0 = (p[1]*(1e6)*Qc)./(p[1]*1e6+Qc)
    # Some variable I arbitrarily defined
    k = 2im*Q0.*(freq - p[5])/p[5]
    # The fitting function in all it's glory
    return (s21AdjMin + k).*exp(1im*p[6])./(1 + k).*totaloffset
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
`Qc` - coupling Q in millions\n
`Q0` - loaded Q in millions\n
`Qi` - intrinsic Q in millions, equal to `p[1]`
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
    Qc = (abs(scalefactor(f0, p)*s21AdjMin(p)*p[1]*1e6/(1 - s21AdjMin(p))))/1e6
    # Calculate the loaded Q
    Q0 = ((p[1]*(1e6).*Qc*1e6)./(p[1]*1e6 + Qc*1e6))/1e6
    Qi = p[1]
    return abs(Qc), abs(Q0), abs(Qi)
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
    energyInternal = resonatorEnergy(10, 5, f0, PACdB, Qi*1e6, Qc*1e6)
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
    capacitancePerUnitLength = 4*vacPerm*dielConstApprox*ceifk(k1, 24)/ceifk(k2, 24) #Rami Thesis Eq 3.3
    length = (3e8/sqrt(dielConstApprox)/f0/4)
    powerReadout = 10.^(powerAtFeedline/10)
    powerInternal = (2/pi)*(Ql.^2 / Qc)*powerReadout # Rami's Thesis eq 3.30
    powerInternaldBm = powerInternal/1000
    voltageResonator = 2*(50*powerInternaldBm).^.5
    energyInternal = .5*capacitancePerUnitLength*voltageResonator.^2*length
    return energyInternal
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



end # module
