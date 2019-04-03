reset
#unset mouse
#et term tek410x
unset key
set style data li
#
#set xrange [-8:8]
#set yrange [-6:6]
#
EV = 1
set multiplot layout 2,2

set zeroaxis
#
set key bottom left
plot [:] '../000.dat'  ev EV us 2:3 wi po pt 6 title "Locus of Position", \
         '../000.dat'  ev EV us 7:8 wi li lt 3 title "Momentum"
#
set key top left
plot [:] '../000.dat' \
	ev EV us ($1/2/pi):4 title "Variance in x", \
''	ev EV us ($1/2/pi):5 title             "y", \
''	ev EV us ($1/2/pi):6 title      "position"
#
set key top right
plot [:] '../000.dat' \
   ev EV us ($1/2/pi):11      ti "Variance in mv", \
'' ev EV us ($1/2/pi):17 lt 3 ti              "P"
#  ev EV us ($1/2/pi):9       ti "Variance in mv_x", \
#  ev EV us ($1/2/pi):10      ti "Variance in mv_y", \
#  ev EV us ($1/2/pi):15 lt 3 ti              "P_x"
#  ev EV us ($1/2/pi):16 lt 3 ti              "P_y"
#
set format y  "%3.0te%+T"
set format y2 "%3.0te%+T"
set y2tics
#
set key below
set ytics nomirror
#set xlabel 'time'
set ylabel 'Error in energy & particle'
set y2label 'Error in P_x'
plot [:] '../000.dat' ev EV us ($1/2/pi):13 wi li                ti 'Error in energy' , \
         '../000.dat' ev EV us ($1/2/pi):14 wi do lt 4 axes x1y2 ti 'P_x (right y-axis)',\
         '../000.dat' ev EV us ($1/2/pi):18 wi li lt 3           ti 'particle'
#
unset key
#
unset multi
pause 6
#pause -1
reread
