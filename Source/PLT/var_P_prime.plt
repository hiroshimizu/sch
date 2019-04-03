reset
set term tek410x
#unset key
#set data style li
#
#set xrange [-8:8]
#set yrange [-6:6]
#
plot [:] '../000.dat' us ($1/2/pi):13 wi lp lt 5 ti 'var_P_prime'
