
<a id='ResonatorFit.jl-Documentation-1'></a>

# ResonatorFit.jl Documentation


<a id='Functions-1'></a>

## Functions

<a id='ResonatorFit.s2p' href='#ResonatorFit.s2p'>#</a>
**`ResonatorFit.s2p`** &mdash; *Function*.

---


```
s2p(file::AbstractString)
```

Imports arrays of frequency and s parameters from a .s2p file

__Arg(s)__:

`file` - Touchstone file with extension .s2p

__Output(s)__:

`frequency`, `S11`, `S21`, `S12`, `S22` - arrays of freqeuency and two port s paramters

__Example__:

```
f, s11, s21, s12, s22 = s2p("touchstone.s2p")
```

<a id='ResonatorFit.s21Fit' href='#ResonatorFit.s21Fit'>#</a>
**`ResonatorFit.s21Fit`** &mdash; *Function*.

---


```
s21Fit(freq, s21; pGuess = [0.5], smoothingFactor = 0, iter = 2000, funcTol = 1e-10 )
```

This function will fit a resonator's S21 transmission response to a fitting model described in Megrant et. al. ArXiv 1201.3384 (2012). Minimizes the square residuals using Nelder Mead optimizing provided in `Optim.jl`.

__Arg(s)__:

`freq` - real array of frequencies

`s21` - complex array of the S21 parameters

(`pGuess`) - array of length 6 for the initial fit parameter guess. The components for fitting are as follows:

  * Qi/1e6

  * dBOffsetFromZero

  * strayInductance(nH)

  * dipDepth (magnitude)

  * f0(GHz)

  * phaseOnResonance

(`smoothingFactor`) - smooths the s21 input using Savitzky–Golay filter with a window size equal to the smoothing factor

(`iter`) - maximum number of iterations for convergence before optimization gives upgrade_scalar

(`funcTol`) - tolerance for the change in the fit function value between successive iterations

__Output(s)__:

`s21fit` - the fit s21 function given from the same freqeuency input and optimized parameters

`s21smoothed` - smoothed s21 function (if not, then original)

`fitResults` - Optim.jl fitting results where calling fitResults.minimum will give optimized paramters

`fitScore` - value of fitting function at optimized parameters

__Example__:

```
s21fit, s21smoothed, fitResults, fitScore = s21Fit(freq, s21)
```

<a id='ResonatorFit.s21FitFunc' href='#ResonatorFit.s21FitFunc'>#</a>
**`ResonatorFit.s21FitFunc`** &mdash; *Function*.

---


```
s21FitFunc(freq, p)
```

This function defines the resonance behavior in a parameterized functional form. The theory for this model can be found in section 2.3.3 and 2.4 of Ben Mazin's thesis

__Arg(s)__:

`freq` - frequency value or array in Hz

`p` -  array of length 6 for fit parameters. The components for fitting are as follows:

  * Qi/1e6

  * dBOffsetFromZero

  * strayInductance(nH)

  * dipDepth (magnitude)

  * f0(GHz)

  * phaseOnResonance

`Qi`, `strayInductance`, and `f0` were choosen to be represented in Millions, nanoHenries, and GHz respectively to due to their scales (1e6, 1e-9, 1e9) relative to the other parameters. __Output(s)__:

`s21val` - the complex value of s21 evaluated at `freq` given fit parameters

<a id='ResonatorFit.getQFromFitParams' href='#ResonatorFit.getQFromFitParams'>#</a>
**`ResonatorFit.getQFromFitParams`** &mdash; *Function*.

---


```
getQFromFitParams(f0, p)
```

This function defines the resonance behavior in a parameterized functional form. The theory for this model can be found in section 2.3.3 and 2.4 of Ben Mazin's thesis.

__Arg(s)__:

`f0` - resonance frequency

`p` -  array of length 6 for fit parameters. The components for fitting are as follows:

  * Qi/1e6

  * dBOffsetFromZero

  * strayInductance(nH)

  * dipDepth (magnitude)

  * f0(GHz)

  * phaseOnResonance

__Output(s)__:

`Qc` - coupling Q

`Q0` - loaded Q

`Qi` - intrinsic Q/1e6, equal to `p[1]`

<a id='ResonatorFit.fitGrade' href='#ResonatorFit.fitGrade'>#</a>
**`ResonatorFit.fitGrade`** &mdash; *Function*.

---


```
fitGrade(fitScore, pts)
```

Give an arbitrary grade based on the fit score.

__Arg(s)__:

`fitScore` - square residual val from minimization function 

`pts` -  length of fitting function, used so that the grade is normalized to the length of the fit

__Output(s)__:

`grade` - welcome back to school


<a id='Utilities-1'></a>

## Utilities

<a id='ResonatorFit.searchdir' href='#ResonatorFit.searchdir'>#</a>
**`ResonatorFit.searchdir`** &mdash; *Function*.

---


```
searchdir(path,key)
```

Search for files matching key in path.

__Arg(s)__:

`path` - path to search

`key` - keyword for file(s) to search

__Output(s)__:

`files` - list of files in path matching key

__Example__:

```
juliaFiles = searchdir(path,".jl")
```

<a id='ResonatorFit.feedOutArrayFromDataFrame' href='#ResonatorFit.feedOutArrayFromDataFrame'>#</a>
**`ResonatorFit.feedOutArrayFromDataFrame`** &mdash; *Function*.

---


```
feedOutArrayFromDataFrame(data::DataFrame)
```

Reformats DataFrame as a sequence of column arrays output for function ouput.

__Arg(s)__:

`data` - DataArray of real parts

__Output(s)__:

`arrayOfArrays` - 1d array containing arrays extracted from DataFrame

<a id='ResonatorFit.complexMapSNP' href='#ResonatorFit.complexMapSNP'>#</a>
**`ResonatorFit.complexMapSNP`** &mdash; *Function*.

---


```
complexMapSNP(N)
```

Complex function map for sNp dataframes. A row with one element equal to one (sum = 1) corresponds to doing nothing (freqeuency component). A row with two elements (one of which is equal to one (real part), and the other which is equal to negative two (imaginary part)) so that sum = -1, corresponds to the combination of two seperate complex elements into one complex number.

__Arg(s)__:

`N` - integer for a sNp (N port) measurement

__Output(s)__:

`map` - 2d array

<a id='Base.complex' href='#Base.complex'>#</a>
**`Base.complex`** &mdash; *Function*.

---


```
complex(r, [i])
```

Convert real numbers or arrays to complex. `i` defaults to zero.

```
complex(data::DataFrame, outNames::Array{Symbol,1},  cMap::Array{Int, 2})
```

Complex function written for the DataFrame type.

__Arg(s)__:

`data` - input of the DataFrame type

`outNames` - array of symbols corresponding output data column labels

`cMap` - complex mapping for sNp files, generated by ResonatorFit.complexMapSNP(N)

__Output(s)__:

`dataFrame` - complex mapped DataFrame

__Example__:

```
complexNames = [:frequency, :S11, :S21, :S12, :S22]
dfOut = complex(dfIn, complexNames, complexMapSNP(2))
```

```
complex(x::DataArray, y::DataArray)
```

Complex function written for the DataArray type.

__Arg(s)__:

`x` - DataArray of real parts

`y` - DataArray of imaginary parts

__Output(s)__:

`output` - complex mapped DataArray

__Example__:

```
complexRand = complex(DataArray(rand(10)), DataArray(1im*rand(10)))
```

<a id='ResonatorFit.complex2dB' href='#ResonatorFit.complex2dB'>#</a>
**`ResonatorFit.complex2dB`** &mdash; *Function*.

---


```
complex2dB(x)
```

Converts a complex number to a magnitube in dB.

__Arg(s)__:

`x` - complex number or array

__Output(s)__:

`db` - magnitude in dB

__Example__:

```
db = complex2dB([1 + 2im;1 - 2im])
```

<a id='ResonatorFit.complex2Phase' href='#ResonatorFit.complex2Phase'>#</a>
**`ResonatorFit.complex2Phase`** &mdash; *Function*.

---


```
complex2Phase(x)
```

Converts a complex number to a phase relative to the x axis.

__Arg(s)__:

`x` - complex number or array

__Output(s)__:

`ph` - phase relative to the x axis

__Example__:

```
ph = complex2Phase([1 + 2im;1 - 2im])
```

<a id='ResonatorFit.dB2mag' href='#ResonatorFit.dB2mag'>#</a>
**`ResonatorFit.dB2mag`** &mdash; *Function*.

---


```
dB2mag(db)
```

Converts from dB to magnitude.

__Arg(s)__:

`db` - dB input

__Output(s)__:

`mag` - magnitude output

<a id='ResonatorFit.mag2dB' href='#ResonatorFit.mag2dB'>#</a>
**`ResonatorFit.mag2dB`** &mdash; *Function*.

---


```
mag2dB(mag)
```

Converts from magnitude to dB.

__Arg(s)__:

`mag` - mag input

__Output(s)__:

`db` - dB output

<a id='ResonatorFit.savitsky_golay' href='#ResonatorFit.savitsky_golay'>#</a>
**`ResonatorFit.savitsky_golay`** &mdash; *Function*.

---


```
savitsky_golay(x::Array{Complex{Float64},1}, windowSize::Int64, polyOrder::Int64)
```

Savitsky Golay filters for complex arrays. __Arg(s)__:

`x` - the complex values of the time history of the signal

`windowSize` -  size of smoothing window. Must be an odd integer number.

`polyOrder` - the order of the polynomial used in the filtering. Must be less then `window_size` - 1.

__Output(s)__:

`ys` - complex array of the smoothed signal 

```
savitsky_golay(x::Array{Float64, 1}, windowSize::Int64, polyOrder::Int64; deriv::Int64=0)
```

Polynomial smoothing with the Savitsky Golay filters. The Savitzky-Golay is a type of low-pass filter, particularly suited for smoothing noisy data. The main idea behind this approach is to make for each point a least-square fit with a polynomial of high order over a odd-sized window centered at the point.

Orginal author BBN (https://github.com/BBN-Q/Qlab.jl). 

Theory: http://www.ece.rutgers.edu/~orfanidi/intro2sp/orfanidis-i2sp.pdf

Python Example: http://wiki.scipy.org/Cookbook/SavitzkyGolay

__Arg(s)__:

`x` - the values of the time history of the signal

`windowSize` -  size of smoothing window. Must be an odd integer number.

`polyOrder` - the order of the polynomial used in the filtering. Must be less then `window_size` - 1.

`deriv`- the order of the derivative to compute (default = 0 means only smoothing)

__Output(s)__:

`ys` - Array of the smoothed signal (or it's n-th derivative)


<a id='Index-1'></a>

## Index

- [`Base.complex`](index.md#Base.complex)
- [`ResonatorFit.complex2Phase`](index.md#ResonatorFit.complex2Phase)
- [`ResonatorFit.complex2dB`](index.md#ResonatorFit.complex2dB)
- [`ResonatorFit.complexMapSNP`](index.md#ResonatorFit.complexMapSNP)
- [`ResonatorFit.dB2mag`](index.md#ResonatorFit.dB2mag)
- [`ResonatorFit.feedOutArrayFromDataFrame`](index.md#ResonatorFit.feedOutArrayFromDataFrame)
- [`ResonatorFit.fitGrade`](index.md#ResonatorFit.fitGrade)
- [`ResonatorFit.getQFromFitParams`](index.md#ResonatorFit.getQFromFitParams)
- [`ResonatorFit.mag2dB`](index.md#ResonatorFit.mag2dB)
- [`ResonatorFit.s21Fit`](index.md#ResonatorFit.s21Fit)
- [`ResonatorFit.s21FitFunc`](index.md#ResonatorFit.s21FitFunc)
- [`ResonatorFit.s2p`](index.md#ResonatorFit.s2p)
- [`ResonatorFit.savitsky_golay`](index.md#ResonatorFit.savitsky_golay)
- [`ResonatorFit.searchdir`](index.md#ResonatorFit.searchdir)
