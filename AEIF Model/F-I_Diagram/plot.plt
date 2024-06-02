
reset 
set terminal postscript eps enhanced size 10,6  color font 'Times-Roman, 40pt'
#set terminal png transparent truecolor size 1500,1080 enhanced font 'Times, 32'
set output "FxI.eps"

reset
set xlabel 'I (pA)' font 'Times-Roman, 50pt' offset 0,0,0
set ylabel 'F (Hz)' font 'Times-Roman, 50pt'offset 0,0,0
set border lw 1
set border 3
set xrange [0 : 1250]
set yrange [-1 : 50]
set ytics 10 nomirror
set xtics 250 nomirror
set samples 20
set palette defined (0 "black", 1 "blue", 2 "magenta", 3 "brown", 4 "orange", 5 "red" )
unset colorbox

#set key left
set key at 1300, 23
set key reverse Left

#plot 'FxI_Diagram.dat' u 1:2:3 w lines lw 15 lc palette notitle

plot 'FxI_Diagram.dat' u 1:2:3 w lines lw 15 lc palette notitle, \
'FxI_Diagram.dat' u 1:($3==5 ? $2: 1/0) w points pt 9 ps 5 lc 'red' title "{/Symbol g} (35+ Hz)", \
'FxI_Diagram.dat' u 1:($3==4 ? $2: 1/0) w points pt 13 ps 5 lc 'orange' title "{/Symbol b} (13 - 35 Hz)", \
'FxI_Diagram.dat' u 1:($3==3 ? $2: 1/0) w points pt 7 ps 5 lc 'brown' title "{/Symbol a} (8 - 13 Hz)", \
'FxI_Diagram.dat' u 1:($3==2 ? $2: 1/0) w points pt 15 ps 5 lc 'magenta' title "{/Symbol q} (4 - 8 Hz)", \
'FxI_Diagram.dat' u 1:($3==1 ? $2: 1/0) w points pt 5 ps 5 lc 'blue' title "{/Symbol d} (0.5 - 4 Hz)", \
'FxI_Diagram.dat' u 1:($3==0 ? $2: 1/0) w lines lw 15 lc 'black' title "Not spiking"

set output
