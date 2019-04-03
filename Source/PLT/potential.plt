reset
set pm3d map
set term gif
set out 'potential.gif'
set xrange [-39:39]
set yrange [-30:30]
splot '../potential.dat'
set out