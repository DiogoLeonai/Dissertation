
reset 
set terminal postscript eps enhanced size 12,7  color font 'Times-Roman, 50pt'
set output "Bursting_Region.eps"
reset
set xlabel '{/Times-Italic V}_{reset} (mV)' font 'Times-Roman, 55pt' offset 0, 0.5
set ylabel '{/Times-Italic b} (pA)' font 'Times-Roman, 55pt' offset 1.2,0
set xrange [-60 : -45]
set yrange [0 : 350]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set cbtics 1
set cbtics add ("Adaptive spike" 1, "Tonic spike" 2, "Initial burst" 3, "Regular burst" 4, "Irregular burst" 5)
set xtics 5
set ytics 80
unset colorbox
set palette defined ( 1 "#f2dbbf", 2 "forest-green", 3 "orange", 4 "red", 5 "black") maxcolors 5
plot 'Bursting_Region.dat' u 1:2:3 with image notitle


set output
