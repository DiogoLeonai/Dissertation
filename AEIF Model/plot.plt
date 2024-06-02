
reset 
set terminal postscript eps enhanced size 12,5  color font 'Times-Roman, 50pt'
set output "ISI.eps"
reset
set xlabel 'time (ms)' font 'Times-Roman, 55pt' offset 0, 0.5
set ylabel 'V (mV)' font 'Times-Roman, 55pt' offset 1.2,0
set xrange [0 : 300]
set yrange [-80 : 0]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set xtics 75
set ytics 25
set border 3


plot 'AEIF_model.dat' u 1:2 with lines lc 'blue' lw 15 notitle



set output
