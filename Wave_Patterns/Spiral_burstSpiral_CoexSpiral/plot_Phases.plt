reset 
#set terminal postscript eps enhanced size 12,6  color font 'Times-Roman, 50pt'
set terminal png transparent truecolor size 1500,1080 enhanced font 'Times, 28pt'
set output "FIG_FASE.png"
set multiplot layout 2, 3
############################################
#################################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.3, 0.8
set origin 0.05, 0.2
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
splot 'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25068)/1000):100 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25075)/1000):110 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25085)/1000):120 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25095)/1000):125 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R70.0_gexc0.1200.dat' u 1:2:(($3+25105)/1000):140 w points pt 7 ps 1 lc palette notitle

############################################
#################################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.3, 0.8
set origin 0.4, 0.2
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
splot 'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25068)/1000):85 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25075)/1000):95 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25085)/1000):105 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25095)/1000):115 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R75.0_gexc0.1400.dat' u 1:2:(($3+25105)/1000):125 w points pt 7 ps 1 lc palette notitle

############################################
##################################
###########################################
reset
set border lw 2
set xrange [0 : 1000]
set yrange [0 : 1000]
set zrange [25.170: 25.21]
set cbrange [0 : 6.28]
set ztics .01 nomirror
set cbtics ( 0, '2{/Symbol p}' 2*pi)
set size 0.3, 0.8
set origin 0.75, 0.2
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
splot 'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25068)/1000):85 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25075)/1000):90 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25085)/1000):95 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25095)/1000):100 w points pt 7 ps 1 lc palette notitle, \
    'FASE_R30.0_gexc0.1500.dat' u 1:2:(($3+25105)/1000):110 w points pt 7 ps 1 lc palette notitle


#################################
#################################
#################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set ylabel 'y ({/Symbol m}m)' font 'Times-Roman, 40pt'
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.4, 0.35
set origin 0 ,0.1
#set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
unset colorbox
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset cblabel
#set title '(a)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R70.0_gexc0.1200.dat' u 1:2:3 w image notitle

#################################
#################################
#################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
#set ylabel 'y ({/Symbol m}m)' font 'Times-Roman, 40pt'
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.4, 0.35
set origin 0.3 ,0.1
#set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
unset colorbox
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset cblabel
#set title '(a)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R75.0_gexc0.1400.dat' u 1:2:3 w image notitle

#################################
#################################
#################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
#set ylabel 'y ({/Symbol m}m)' font 'Times-Roman, 40pt'
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.4, 0.35
set origin 0.65 ,0.1
#set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
unset colorbox
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset cblabel
#set title '(a)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R30.0_gexc0.1500.dat' u 1:2:3 w image notitle





unset multiplot
set output
