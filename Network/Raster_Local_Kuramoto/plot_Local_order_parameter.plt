
reset 
set terminal postscript eps enhanced size 14,5  color font 'Times-Roman, 50pt'
set output "Local_Kuramoto.eps"
set multiplot layout 1,2

reset
set xlabel 'x ({/Symbol m}m)' font 'Times-Roman, 55pt' offset 0, 0.5
set ylabel 'y ({/Symbol m}m)' font 'Times-Roman, 55pt' offset 1.2,0
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250 nomirror
set ytics 250 nomirror
set cbtics 0.1 font 'Times-Roman, 40pt' nomirror
set cbrange [0.7: 1]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror 
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset colorbox
set size 0.49, 1
set title '(a)' offset  -16,0
plot 'LocalZ_R20.0_gexc_0.010.dat' u 1:2:3 w image notitle

set cblabel '~{Z}{0.9\_} _{/:Italic l,m}' font 'Times-Roman, 55pt' rotate by 0 offset -0.5,0
set colorbox
unset ylabel
set size 0.55, 1
set origin 0.45,0
set title '(b)' offset -16,0
plot 'LocalZ_R70.0_gexc_0.050.dat' u 1:2:3 w image notitle



unset multiplot
set output
