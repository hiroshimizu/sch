reset
#et term tek410x
unset key
set style data li
#
#set xrange [-8:8]
#set yrange [-6:6]
#
EV = 1
set multiplot layout 4,2
set title "position"
plot [:] '../000.dat'  ev EV us 2:3
#
set title "Deviation in Position"
plot [:] '../000.dat' ev EV us ($1/2/pi):6
#
set title "mv Momentum"
plot [:] '../000.dat'  ev EV us 7:8
#
set title "Deviation in mv Momentum"
plot [:] '../000.dat' ev EV us ($1/2/pi):9
#
set title "Energy_Error"
plot [:] '../000.dat' ev EV us ($1/2/pi):11 wi lp
#
set title "Px"
plot [:] '../000.dat' ev EV us ($1/2/pi):12 wi lp
#
set title "var_P"
plot [:] '../000.dat' ev EV us ($1/2/pi):13
#
unset multi
pause 6
#pause -1
reread
