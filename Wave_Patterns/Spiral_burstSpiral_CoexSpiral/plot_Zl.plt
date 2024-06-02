reset 
set terminal postscript eps enhanced size 14,6  color font 'Times-Roman, 32pt'
#set terminal png transparent truecolor size 1500,1080 enhanced font 'Times, 32'
set output "FIG_Zl.eps"
set multiplot layout 1, 6
############################################
###############R = 55 , g = 0.16###################
###########################################
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
set size 0.3, 0.53
set origin 0 ,0.52
#set colorbox user origin 0.23, 0.95 size .15, .02 horizontal
unset colorbox
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset cblabel
set title '(a)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R55.0_gexc0.1600.dat' u 1:2:3 w image notitle

############################################
###############R = 32 , g = 0.15###################
###########################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.3, 0.53
set origin 0.32, 0.52
set colorbox user origin 0.42, 0.38 size .02, .3
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset colorbox
set title '(b)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R32.0_gexc0.1500.dat' u 1:2:3 w image notitle

############################################
###############R = 70 , g = 0.12###################
###########################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.355, 0.53
set origin 0.64, 0.52
#set colorbox user origin 0.42, 0.38 size .02, .3
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
set cbtic offset -0.6, 0
set cblabel '~{Z}{0.9\_}_{/:Italic l,m}' font 'Times-Roman, 40pt' rotate by 0
set title '(c)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R70.0_gexc0.1200.dat' u 1:2:3 w image notitle

############################################
###############R = 50 , g = 0.24############
###########################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set ylabel 'y ({/Symbol m}m)' font 'Times-Roman, 40pt'
set xlabel 'x ({/Symbol m}m)' font 'Times-Roman, 40pt'
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.3, 0.56
set origin 0.0, 0
set colorbox user origin 0.42, 0.38 size .02, .3
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset colorbox
set title '(d)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R50.0_gexc0.2400.dat' u 1:2:3 w image notitle

############################################
###############R = 75 , g = 0.14###################
###########################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set xlabel 'x ({/Symbol m}m)' font 'Times-Roman, 40pt'
set cbtics 0.1
set size 0.3, 0.56
set origin 0.32, 0
set colorbox user origin 0.42, 0.38 size .02, .3
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
unset colorbox
set title '(e)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R75.0_gexc0.1400.dat' u 1:2:3 w image notitle

############################################
###############R = 30 , g = 0.15###################
###########################################
reset
set border lw 1
set xrange [0 : 976]
set yrange [0 : 976]
set xtics 250
set ytics 250
set xtics add ('1000' 976)
set ytics add ('1000' 976)
set xlabel 'x ({/Symbol m}m)' font 'Times-Roman, 40pt'
set cbrange[0.7: 1]
set cbtics 0.1
set size 0.355, 0.56
set origin 0.64, 0
set palette defined (0.7 "#000000", 0.8 "#08bdbd", 0.9 "#a30000", 1.0 "yellow")
set cblabel '~{Z}{0.9\_}_{/:Italic l,m}' font 'Times-Roman, 40pt' rotate by 0
set cbtic offset -0.6, 0
set title '(f)' font 'Times-Roman, 40pt' offset -17,-0.5
plot 'Local_Kuramoto_R30.0_gexc0.1500.dat' u 1:2:3 w image notitle






unset multiplot
set output
