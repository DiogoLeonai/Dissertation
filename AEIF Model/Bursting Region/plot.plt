
reset 
set terminal postscript eps enhanced size 12,7  color font 'Times-Roman, 50pt'
set output "Bursting_Region.eps"
reset
set xlabel '{/Times-Italic V}_{reset} (mV)' font 'Times-Roman, 55pt' offset 0, 0.5
set ylabel '{/Times-Italic b} (pA)' font 'Times-Roman, 55pt' offset 1.2,0
set xrange [-60 : -45]
set yrange [0 : 240]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set cbtics 1
set cbtics add ("Burst" 1, "Spike" 0)
set xtics 5
set ytics 80
unset colorbox
set palette defined ( 0 "#f2dbbf", 1 "dark-red") maxcolors 2
plot 'Bursting_Region.dat' u 1:2:4 with image notitle


set output
