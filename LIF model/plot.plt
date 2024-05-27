
reset 
set terminal postscript eps enhanced size 12,4  color font 'Times-Roman, 50pt'
set output "LIF_NoReset_and_Reset.eps"

reset
set xlabel 'time (ms)' font 'Times-Roman, 50pt' offset 0, 0.5
set ylabel 'V (mV)' font 'Times-Roman, 50pt' offset 1.2,0
set xrange [0 : 500]
set yrange [-75 : -25]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror
set border 3
set xtics 250 nomirror
set ytics 25 nomirror
set multiplot layout 1,2

set title '(a)' offset -15,-0.5
set size 0.5, 1.05
plot 'LIF_NoReset.dat' u 1:2  with lines lc 'black' lw 8 notitle

set title '(b)' offset -15,-0.5
set size 0.5, 1.05
plot 'LIF.dat' u 1:2 with lines lc 'blue' lw 8 notitle

unset multiplot

set output
