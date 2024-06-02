
reset 
set terminal postscript eps enhanced size 14,8  color font 'Times-Roman, 50pt'
set output "Rasterplot_Z(t).eps"
set multiplot layout 2,2


########Ploting the raster-plot
reset
set xlabel 'time (ms)' font 'Times-Roman, 55pt' offset 0, 0.5
set ylabel 'i' font 'Times-Roman, 55pt' offset 2.3,0
set xrange [0 : 500]
set yrange [0 : 17324]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
#set xtics 200
#set xtics add ('0' 25000, '200' 25200, '400' 25400)

set ytics 8000
set border 3
set ytics add ('8x10^3' 8000, '16x10^3' 16000)
set origin -0.0, 0.0
#set title 'Z'
set size 0.55,0.5
plot 'Raster_R20.0_gexc_0.010.dat' every 10 u ($1 - 25200):2 with points pt 7 ps 0.3 lc 'blue' notitle

set origin 0.47, 0.0
set size 0.55,0.5
#set title 'CV = 0.863'
plot 'Raster_R70.0_gexc_0.050.dat' every 10 u ($1 - 25200):2 with points pt 7 ps 0.3 lc 'red' notitle

#########Plotting the Z(t)

reset 
set xlabel 'time' font 'Times-Roman, 55pt' offset 0, 0.3
set ylabel 'Z(t)' font 'Times-Roman, 55pt' offset 1.8,0
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set border 3
set xrange [0 : 500]

set origin -0.03, 0.6
set yrange [0: 1.1]
set size 0.52, 0.3
set ytics 0.25 
#set title '(a)' offset -15.3,10
plot 'Z(t)_R20.0_gexc_0.010.dat' u ($1 - 25200):2 with lines lw 9  lc 'blue' notitle

set origin 0.48, 0.6
set size 0.52 , 0.3
#set yrange [0.75: 1.1]
#set title '(b)' offset -15.3,10
plot 'Z(t)_R70.0_gexc_0.050.dat' u ($1 - 25200):2 with lines lw 9 lc 'red' notitle

unset multiplot
set output
