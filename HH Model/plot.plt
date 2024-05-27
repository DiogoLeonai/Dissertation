
reset 
set terminal postscript eps enhanced size 14,8  color font 'Times-Roman, 50pt'
set output "HH_model.eps"
set multiplot layout 1,2

#####Current
reset
set xlabel 'time (ms)' font 'Times-Roman, 50pt' offset 0, 0.5
set ylabel 'I ({/Symbol m} A/cm^2)' font 'Times-Roman, 50pt' offset 1.2,0
set xrange [120 : 270]
set yrange [-1 : 10]
set ytics 5
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror
set border 2
set border lw 2
unset xlabel
unset xtics
set size 0.61,0.35
set origin -0.01, 0.7
set title '(a)' font 'Times-Roman, 50pt' offset -21, -0.005
plot 'HH_model.dat' u 1:2  with lines lc 'magenta' lw 14 notitle

#########Plotting the potential trace
reset
set xlabel 'time (ms)' font 'Times-Roman, 50pt' offset 0, 0.5
set ylabel 'V (mV)' font 'Times-Roman, 50pt' offset 2.3,0
set xrange [120 : 270]
set yrange [-30 : 120]
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror
set border lw 2
set border 3
set xtics 60 nomirror
set ytics 30 nomirror
set size 0.65,0.75
set origin -0.05, 0.0
set xtics add ('0' 120, '' 180, '75' 195, '' 240, '150' 270)
plot 'HH_model.dat' u 1:3  with lines lc 'black' lw 14 notitle


#########Plotting n, m and h variables

######n
reset
set ylabel 'n' font 'Times-Roman, 50pt' offset 2,0
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror
set xtics 100 nomirror
set ytics 0.5 nomirror
set border lw 2
set border 3
set yrange[0:1]
set yrange[0:1]
set xrange [120 : 270]
set size 0.45, 0.43
set origin 0.55, 0.62
set xtics add ('' 200, '' 300)
set title '(b)' font 'Times-Roman, 50pt' offset -15, -0.01
plot 'HH_model.dat' u 1:4  with lines lc 'blue' lw 14 notitle

######m
reset
set ylabel 'm' font 'Times-Roman, 50pt' offset 2,0
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror
set xtics 100 nomirror
set ytics 0.5 nomirror
set border lw 2
set border 3
set yrange[0:1]
set yrange[0:1]
set xrange [120 : 270]
set size 0.45, 0.35
set origin 0.55, 0.35
set xtics add ('' 200, '' 300)
plot 'HH_model.dat' u 1:5  with lines lc 'red' lw 14 notitle

#####h
reset
set ylabel 'h' font 'Times-Roman, 50pt' offset 2,0
set xtics font 'Times-Roman, 40pt' nomirror
set ytics font 'Times-Roman, 40pt' nomirror
set xtics 60 nomirror
set ytics 0.5 nomirror
set border lw 2
set border 3
set yrange[0:1]
set xrange [120 : 270]
set size 0.45, 0.4
set origin 0.55, -0.0
set xtics add ('0' 120, '' 180, '75' 195, '' 240, '150' 270)
set xlabel 'time (ms)' font 'Times-Roman, 50pt' offset 0, 0.5
plot 'HH_model.dat' u 1:6  with lines lc 'forest-green' lw 14 notitle

unset multiplot
set output
