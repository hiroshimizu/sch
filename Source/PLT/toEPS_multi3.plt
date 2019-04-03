reset
#unset mouse
#et term tek410x
set term postscript eps enhanced color "TimesNewRoman" 10
set out "locus_info.eps"
unset key
set style data li
#
#set xrange [-8:8]
#set yrange [-6:6]
#
EV = 1
E2 = 3
set multiplot layout 2,2

set zeroaxis
#
set key bottom left
plot [:] '../000.dat'  ev EV us 2:3 wi po pt 6 lc 1 title "Locus of Position", \
         '../000.dat'  ev EV us 7:8 wi li lt 1 lc 3 title "Momentum"
#
set key top left
plot [:] '../000.dat' \
	ev EV us ($1/2/pi):4 lt 1 lc 1 title "Variance in {/TimesNewRoman-Italic x}", \
''	ev EV us ($1/2/pi):5 lt 1 lc 2 title             "{/TimesNewRoman-Italic y}", \
''	ev EV us ($1/2/pi):6 lt 1 lc 3 title    								   "position"
#
set key top right
plot [:] '../000.dat' \
   ev EV us ($1/2/pi):11 wi li lt 1        lc 1 ti "Variance in <{/TimesNewRoman-Italic mv}>", \
'' ev EV us ($1/2/pi):17 wi li lt 1        lc 3 ti              "<{/TimesNewRoman^Italic P}>", \
'' ev E2 us ($1/2/pi):9  wi po pt 1 ps 0.5 lc 1 ti "Variance in <{/TimesNewRoman-Italic mv_x}>", \
'' ev E2 us ($1/2/pi):10 wi po pt 2 ps 0.5 lc 1 ti "Variance in <{/TimesNewRoman-Italic mv_y}>", \
'' ev E2 us ($1/2/pi):15 wi po pt 6 ps 0.5 lc 3 ti "<{/TimesNewRoman-Italic P_x}>", \
'' ev E2 us ($1/2/pi):16 wi po pt 4 ps 0.5 lc 3 ti "<{/TimesNewRoman-Italic P_y}>"
#
set format y  "%3.0te%+T"
set format y2 "%3.0te%+T"
set y2tics
#
set key below
set ytics nomirror
#set xlabel 'time'
set ylabel 'Error in energy & particle'
set y2label 'Error in <{/TimesNewRoman-Italic P_x}>'
plot [:] '../000.dat' ev EV us ($1/2/pi):13 wi li                ti 'Error in energy' , \
         '../000.dat' ev EV us ($1/2/pi):14 wi do axes x1y2 ti '<{/TimesNewRoman-Italic P_x}> (right {/TimesNewRoman-Italic y}-axis)',\
         '../000.dat' ev EV us ($1/2/pi):18 wi li lt 3           ti 'particle'
#
unset key
#
unset multi
#pause 6
#pause -1
#reread
#set term wxt
set term x11
set out
