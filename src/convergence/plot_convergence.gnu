filename = 'AnalyticResults.txt'
set terminal pdf colour
set output "AnalyticAnalysis.pdf"

interpolationTypes = "1 2 3"
dynamicDegrees = "1 2"
timeSteps = "0.100000E+00 0.500000E-01 0.100000E-01 0.500000E-02"

interpolationName = "LinearLagrange QuadraticLagrange CubicLagrange CubicHermite"

test = sprintf("%.6E",0.1)

selectData(it,dd,ts) = sprintf('< awk ''{if( ($1 == "%d" && $3 == "%d" && $6 == "%s")) print }'' %s',it,dd,ts,filename)

#set style data linespoints
set logscale x
set logscale y
set xlabel "Characteristic length, h"
set ylabel "Absolute RMS Error"
set title "RMS Error vs h, Linear timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,1,"0.100000E+00") using 5:7 title word(interpolationName,it)
set title "RMS Error vs h, Linear timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,1,"0.500000E-01") using 5:7 title word(interpolationName,it)
set title "RMS Error vs h, Linear timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,1,"0.100000E-01") using 5:7 title word(interpolationName,it)
set title "RMS Error vs h, Linear timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,1,"0.500000E-02") using 5:7 title word(interpolationName,it)

set xlabel "Characteristic length, h"
set ylabel "Integral Error"
set title "Integral vs h, Linear timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,1,"0.100000E+00") using 5:8 title word(interpolationName,it)
set title "Integral vs h, Linear timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,1,"0.500000E-01") using 5:8 title word(interpolationName,it)
set title "Integral vs h, Linear timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,1,"0.100000E-01") using 5:8 title word(interpolationName,it)
set title "Integral vs h, Linear timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,1,"0.500000E-02") using 5:8 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Number of DOFs"
set ylabel "Absolute RMS Error"
set title "RMS Error vs Number of DOFs, Linear timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,1,"0.100000E+00") using 4:7 title word(interpolationName,it)
set title "RMS Error vs Number of DOFs, Linear timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,1,"0.500000E-01") using 4:7 title word(interpolationName,it)
set title "RMS Error vs Number of DOFs, Linear timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,1,"0.100000E-01") using 4:7 title word(interpolationName,it)
set title "RMS Error vs Number of DOFs, Linear timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,1,"0.500000E-02") using 4:7 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Number of DOFs"
set ylabel "Integral Error"
set title "Integral Error vs Number of DOFs, Linear timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,1,"0.100000E+00") using 4:8 title word(interpolationName,it)
set title "Integral Error vs Number of DOFs, Linear timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,1,"0.500000E-01") using 4:8 title word(interpolationName,it)
set title "Integral Error vs Number of DOFs, Linear timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,1,"0.100000E-01") using 4:8 title word(interpolationName,it)
set title "Integral Error vs Number of DOFs, Linear timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,1,"0.500000E-02") using 4:8 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Characteristic length, h"
set ylabel "Absolute RMS Error"
set title "RMS Error vs h, Quadratic timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,2,"0.100000E+00") using 5:7 title word(interpolationName,it)
set title "RMS Error vs h, Quadratic timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,2,"0.500000E-01") using 5:7 title word(interpolationName,it)
set title "RMS Error vs h, Quadratic timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,2,"0.100000E-01") using 5:7 title word(interpolationName,it)
set title "RMS Error vs h, Quadratic timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,2,"0.500000E-02") using 5:7 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Characteristic length, h"
set ylabel "Integral Error"
set title "Integral Error vs h, Quadratic timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,2,"0.100000E+00") using 5:8 title word(interpolationName,it)
set title "Integral Error vs h, Quadratic timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,2,"0.500000E-01") using 5:8 title word(interpolationName,it)
set title "Integral Error vs h, Quadratic timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,2,"0.100000E-01") using 5:8 title word(interpolationName,it)
set title "Integral Error vs h, Quadratic timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,2,"0.500000E-02") using 5:8 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Number of DOFs"
set ylabel "Absolute RMS Error"
set title "RMS Error vs Number of DOFs, Quadratic timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,2,"0.100000E+00") using 4:7 title word(interpolationName,it)
set title "RMS Error vs Number of DOFs, Quadratic timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,2,"0.500000E-01") using 4:7 title word(interpolationName,it)
set title "RMS Error vs Number of DOFs, Quadratic timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,2,"0.100000E-01") using 4:7 title word(interpolationName,it)
set title "RMS Error vs Number of DOFs, Quadratic timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,2,"0.500000E-02") using 4:7 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Number of DOFs"
set ylabel "Integral Error"
set title "Integral Error vs Number of DOFs, Quadratic timestepping, ts = 0.1"
plot for [it=1:3] selectData(it,2,"0.100000E+00") using 4:8 title word(interpolationName,it)
set title "Integral Error vs Number of DOFs, Quadratic timestepping, ts = 0.05"
plot for [it=1:3] selectData(it,2,"0.500000E-01") using 4:8 title word(interpolationName,it)
set title "Integral Error vs Number of DOFs, Quadratic timestepping, ts = 0.01"
plot for [it=1:3] selectData(it,2,"0.100000E-01") using 4:8 title word(interpolationName,it)
set title "Integral Error vs Number of DOFs, Quadratic timestepping, ts = 0.005"
plot for [it=1:3] selectData(it,2,"0.500000E-02") using 4:8 title word(interpolationName,it)

set logscale x
set logscale y
set xlabel "Delta t"
set ylabel "Absolute RMS Error"
set title "RMS Error vs Delta t, Linear timestepping, Linear interpolation"
plot for [ts in timeSteps] selectData(1,1,ts) using 6:7 title ts
set title "RMS Error vs Delta t, Linear timestepping, Quadratic interpolation"
plot for [ts in timeSteps] selectData(2,1,ts) using 6:7 title ts
set title "RMS Error vs Delta t, Linear timestepping, Cubic interpolation"
plot for [ts in timeSteps] selectData(3,1,ts) using 6:7 title ts



#plot selectData("1","1","0.100000E+00"), selectData("2","1","0.100000E+00")
pause(1.0)