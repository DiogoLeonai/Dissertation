reset 
#set terminal postscript eps enhanced size 12,6  color font 'Times-Roman, 50pt'
set terminal png transparent truecolor size 1500,1080 enhanced font 'Times, 32'
set output "FIG_FASE.png"
set multiplot layout 1, 6
############################################
###############R = 55 , g = 0.16###################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.24, 1.3
set origin -0.05, -0.2
set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
set palette defined (-80 "yellow",-60 "blue",0 "blue",20 "yellow")
set view 80, 47
set ticslevel 0
unset ztics
unset xlabel
unset ylabel
unset xtics
unset ytics
unset zlabel
unset cblabel
splot 'PHASE_R75.0_gexc0.128.dat' u 1:2:(($3+25068)/1000):100 w points pt 7 ps 1 lc palette notitle, \
    'PHASE_R75.0_gexc0.128.dat' u 1:2:(($3+25075)/1000):105 w points pt 7 ps 1 lc palette notitle, \
    'PHASE_R75.0_gexc0.128.dat' u 1:2:(($3+25085)/1000):110 w points pt 7 ps 1 lc palette notitle, \
    'PHASE_R75.0_gexc0.128.dat' u 1:2:(($3+25095)/1000):115 w points pt 7 ps 1 lc palette notitle, \
    'PHASE_R75.0_gexc0.128.dat' u 1:2:(($3+25105)/1000):120 w points pt 7 ps 1 lc palette notitle


############################################
###############R = 32 , g = 0.15###################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.24, 1.3
set origin 0.12, -0.2
set colorbox user origin 0.42, 0.38 size .02, .3
set palette defined (-80 "yellow",-60 "blue",0 "blue",20 "yellow")
set view 80, 47
set ticslevel 0
unset ztics
unset colorbox
unset xlabel
unset ylabel
unset zlabel
unset xtics
unset ytics
splot 'FASE_R32.0_gexc0.1500.dat' u 1:2:(($3+25068)/1000):85 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R32.0_gexc0.1500.dat' u 1:2:(($3+25075)/1000):95 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R32.0_gexc0.1500.dat' u 1:2:(($3+25085)/1000):105 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R32.0_gexc0.1500.dat' u 1:2:(($3+25095)/1000):115 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R32.0_gexc0.1500.dat' u 1:2:(($3+25105)/1000):125 w points pt 7 ps 1 lc palette notitle

############################################
###############R = 70 , g = 0.12###################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.24, 1.3
set origin 0.295, -0.2
set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
set palette defined (-80 "yellow",-60 "blue",0 "blue",20 "yellow")
set view 80, 47
set ticslevel 0
unset ztics
unset xlabel
unset ylabel
unset xtics
unset ytics
unset zlabel
unset colorbox
splot 'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25068)/1000):85 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25075)/1000):90 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25085)/1000):95 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25095)/1000):100 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25105)/1000):110 w points pt 7 ps 1 lc palette notitle

############################################
###############R = 50 , g = 0.24############
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.24, 1.3
set origin 0.47, -0.2
set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
set palette defined (-80 "yellow",-60 "blue",0 "blue",20 "yellow")
set view 80, 47
set ticslevel 0
unset ztics
unset xlabel
unset ylabel
unset xtics
unset ytics
unset zlabel
unset colorbox
splot 'FASE_R50.0_gexc0.2400.dat' u 1:2:(($3+25068)/1000):85 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R50.0_gexc0.2400.dat' u 1:2:(($3+25075)/1000):95 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R50.0_gexc0.2400.dat' u 1:2:(($3+25085)/1000):105 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R50.0_gexc0.2400.dat' u 1:2:(($3+25095)/1000):115 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R50.0_gexc0.2400.dat' u 1:2:(($3+25105)/1000):125 w points pt 7 ps 1 lc palette notitle

############################################
###############R = 75 , g = 0.14###################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.24, 1.3
set origin 0.645, -0.2
set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
set palette defined (-80 "yellow",-60 "blue",0 "blue",20 "yellow")
set view 80, 47
set ticslevel 0
unset ztics
unset xlabel
unset ylabel
unset xtics
unset ytics
unset zlabel
unset colorbox
splot 'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25068)/1000):70 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25075)/1000):90 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25085)/1000):105 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25095)/1000):120 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25105)/1000):130 w points pt 7 ps 1 lc palette notitle

############################################
###############R = 30 , g = 0.15###################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.24, 1.3
set origin 0.81, -0.2
set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
set palette defined (-80 "yellow",-60 "blue",0 "blue",20 "yellow")
set view 80, 47
set ticslevel 0
unset ztics
unset xlabel
unset ylabel
unset xtics
unset ytics
unset zlabel
unset colorbox
splot 'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25068)/1000):80 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25075)/1000):88 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25085)/1000):100 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25095)/1000):110 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25105)/1000):115 w points pt 7 ps 1 lc palette notitle






unset multiplot
set output
