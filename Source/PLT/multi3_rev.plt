reset
N_rot = 5
FILE='../000.dat'
load '../sch.inp'
#oad    'sch.ini'
x0G=0; y0G=0
hbar = h_bar_ * e_ * Bz_ / ( mp_ * v0_ )**2
#if($0==0) print "hbar = ", hbar, ", E0 = ", E0;
#print "$0=",int($0)
set term x11 font "Times, 12"
unset key
unset mouse
set style data li
#
EP = 2**5			# for plot wi points
EL = 1			# for plot wi lines
EM = 1 * EP	# for locus of momentum
set multiplot layout 2,3

set zeroaxis
#
set key top left
set xtics auto
#set xtics 0,0.25
set ytics 0,1
set ytics add("1/2" 0.5)
set ytics add("5/2" 2.5)

plot [0:N_rot][0:] FILE \
  	ev EL us ($1/2/pi):($6/hbar) wi li lt 1 lw 2 ti "Var. r",\
""	ev EL us ($1/2/pi):($4/hbar) wi li lt 4      ti "x", 0.5 wi li lt 0 noti, \
''	ev EL us ($1/2/pi):($5/hbar) wi li lt 3      ti "y", 2.5 wi li lt 0 noti
#pause -1
#
#set ytics auto
set ytics 0,0.5
#et ytics add ("1/8" 0.125)
set key top right
#et key below
plot [0:N_rot][0:1.6] FILE \
   ev EL us ($1/2/pi):($17/hbar) wi li lt 1 lw 2 ti "Var. P", \
'' ev EL us ($1/2/pi):($15/hbar) wi li lt 4      ti "P_x", \
'' ev EL us ($1/2/pi):($16/hbar) wi li lt 3      ti "P_y"
#q
#pause -1
set key bottom right
set ytics 0,0.5
#set ytics add ("5/8" 0.625)
#set ytics add ("5/4" 1.25)
plot [0:N_rot][0:1.6] FILE \
   ev EL us ($1/2/pi):($11/hbar) wi li      lw 2 ti "Var. mv", \
'' ev EL us ($1/2/pi):($9 /hbar) wi li lt 4      ti "mv_x", \
'' ev EL us ($1/2/pi):($10/hbar) wi li           ti "mv_y"
#q
#pause -1
set xtics auto
set ytics auto
set key top left
#set format y  "%3.0te%+T"
#q
set format "%g"
#et key center left
#set xtics 0,1e-14
set ytics nomirror
set y2tics auto
set xlabel 'Error in Paticle Conservation'
set ylabel 'Error in Energy'
set y2label 'Error in P_x'
plot FILE ev EL us 18:13 wi do noti, '' ev EL us 18:14 axes x1y2 wi do noti
unset xlabel
unset ylabel
unset y2label
set xtics auto
#q
unset y2tics
#pause -1
#set format y "%3.0te%+T"
set format "%g"
#et xlabel 'Change in x_G'
#et ylabel 'Change in y_G'
set key bottom left
plot FILE ev EL us 14:(q*Bz*$3-$7) wi po lt 1 lw 2 ti 'y_G'
#q
set xtics auto
set xtics 0.5
set ytics 0.5
set format "%g"
set key center
set xlabel "x or mv_x"
set ylabel "y or mv_y"
plot [:] FILE ev EM us  7: 8 wi po lt 4 pt 4 ti         "mv_x-mv_y", \
			FILE ev EL us  2: 3 wi li lt 1 lw 2 ti "Locus of Position"
         
#q
unset key
#q
unset multi
pause 6
#pause -1
reread

