include("HelperFuncs.jl")

export circleFit

using PyPlot

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
function circleFit(freq, s21; fitPts = 0, plotFit = false)
    if fitPts == 0
        fitPts = length(freq)
    end
    phase = linspace(0, 2*pi, fitPts)
    f0idx = indmin(complex2dB(s21))
    s21inv = 1./s21

    # Find outer dimensions of circle and use to approximate center
    xDim = maximum(real(s21inv)) - minimum(real(s21inv))
    yDim = maximum(imag(s21inv)) - minimum(imag(s21inv))
    r = maximum([xDim, yDim])/2
    center = [maximum(real(s21inv)) - xDim/2, maximum(imag(s21inv)) - yDim/2]
    # Shift circle to orgin
    s21inv = s21inv - (center[1] + 1im*center[2])
    # It's easy to extract upper and lower bound on radius on a centered circle
    rMax = maximum(sqrt(real(s21inv).^2 + imag(s21inv).^2))
    rMin = minimum(sqrt(real(s21inv).^2 + imag(s21inv).^2))

    r = (rMax + rMin)/2

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
        plot(real(s21inv), imag(s21inv), "bo", markersize = 6)
        plot(real(cMatch), imag(cMatch), "ro", markersize = 6)
        plot(xN, bN, "g", lw = 3)
        plot(xP, bP, "g", lw = 3)
        plot(xN, fN, "m", lw = 3)
        plot(xP, fP, "m", lw = 3)
        plot(real(cMatch[f0idx]), imag(cMatch[f0idx]), "mo", markersize = 12)
        plot(real(cMatch[minPhaseNidx]), imag(cMatch[minPhaseNidx]), "go", markersize = 12)
        plot(real(cMatch[minPhasePidx]), imag(cMatch[minPhasePidx]), "go", markersize = 12)
        plot(real(circleImag(phase, [0,0], rMax)),  imag(circleImag(phase, [0,0], rMax)),"k--")
        plot(real(circleImag(phase, [0,0], rMin)),  imag(circleImag(phase, [0,0], rMin)),"k--")
        axis([-rMax; rMax; -rMax; rMax])
        title(string("f0 = ",freq[f0idx],"\nQi = ", QiApprox))
        xlabel("Re 1/S21")
        ylabel("Im 1/S21")
        # legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    end
    return QiApprox
end
