# mvx, mvy, mvを近似的に再現する
set samples 200

f(x)= 0.15*cos(2*(2*pi*x)) + 0.45
g(x)=-0.25*cos(2*(2*pi*x)) + 0.55
set xrange [0:1.6]

plot	f(x)+g(x),\
		f(x) wi po pt 1 ps 0.5 lc 1,\
		g(x) wi po pt 2 ps 0.5 lc 1,\
		0.3 wi do, 0.9 wi do


