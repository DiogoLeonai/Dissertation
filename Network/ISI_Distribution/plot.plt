
reset 
set terminal postscript eps enhanced size 12,8  color font 'Times-Roman, 50pt'
set output "Potential_ISI_distribution.eps"
set multiplot layout 2,2


########Ploting the voltage traces
reset
set xlabel 'time (ms)' font 'Times-Roman, 55pt' offset 0, 0.5
set ylabel 'V (mV)' font 'Times-Roman, 55pt' offset 1.2,0
set xrange [25000 : 25400]
set yrange [-60 : -40]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set xtics 200
set xtics add ('0' 25000, '200' 25200, '400' 25400)
#set ytics 25
set border 3

set origin 0.0, 0.5
set title 'CV = 0.002'
plot 'Potential_R70.0_gexc0.0500.dat' u 1:2 with lines lc 'blue' lw 8 notitle

set origin 0.5, 0.5
set title 'CV = 0.863'
plot 'Potential_R70.0_gexc0.1500.dat' u 1:2 with lines lc 'red' lw 8 notitle

#########Plotting the ISI distributions

reset 
set xlabel 'ISI_i' font 'Times-Roman, 55pt' offset 0, 0.3
set ylabel 'P_{ISI}' font 'Times-Roman, 55pt' offset 1.8,0
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set border 3
set xrange [0:200]

set origin -0.03, 0.0
set yrange [0: 1.1]
set size 0.517 , 0.58
set ytics 0.25 
set title '(a)' offset -15.3,10
plot 'ISI_Histogram_R70.0_gexc0.0500.dat' u 1:2:(8) with boxes fill solid  lc 'blue' notitle

set origin 0.465, 0.0
set yrange [0: 0.15]
set size 0.535 , 0.58
set ytics 0.035
set title '(b)' offset -15.3,10
plot 'ISI_Histogram_R70.0_gexc0.1500.dat' u 1:2:(8) with boxes fill solid lc 'red' notitle

unset multiplot
set output
