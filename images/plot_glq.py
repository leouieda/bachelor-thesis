import Gnuplot
import numpy as np
import pylab as pl
import glq
from plot_tesseroid_coordsys import tessSideGCirc, tessSidePar

gp = Gnuplot.Gnuplot(persist = 1)
output = 'tess_glq.ps'
#gp("set term postscript")
#gp('set output "%s"' % (output))
gp("set mapping spherical")
gp("set angle degrees")
#gp("set xrange [0:1]; set yrange [0:1]; set zrange [0:1]")
#gp('set format x ""')
#gp('set format y ""')
#gp('set format z ""')
gp('unset xtics; unset ytics; unset ztics')
gp('unset border')
gp('set size ratio 1')
gp('set xrange [0:1.2];set yrange [0:1.2];set zrange [0:1.2]')
gp("set ticslevel 0")
#gp('set xlabel "x"; set ylabel "y"; set zlabel "z"')
gp("set view 80,145")
gp("set parametric")
gp("set dummy u,v")
gp("set isosample 6")
gp("set urange [0:90]; set vrange [0:90]")

# Tesseroide
w = 18
e = 36
s = 36
n = 54
t = 1.25
b = 1
sides1 = [tessSideGCirc(s, n, w, t), tessSidePar(w, e, n, t), \
            tessSideGCirc(s, n, e, t), tessSidePar(w, e, s, t), \
            tessSideGCirc(s, n, e, b), tessSidePar(w, e, s, b), \
            [[w,s,t],[w,s,b]], [[e,n,t],[e,n,b]], [[e,s,t],[e,s,b]]]
sides2 = [tessSidePar(w, e, n, b), [[w,n,t],[w,n,b]], tessSideGCirc(s, n, w, b)]
gpdata = []
for line in sides1:
    gpdata.append(Gnuplot.PlotItems.Data(line, with_="l lt 1 lw 1.5", title=None))
for line in sides2:
    gpdata.append(Gnuplot.PlotItems.Data(line, with_="l lt 2 lw 1.5", title=None))

# make the abscissas
order = 3
absr = glq.Abscissas(order)
absr.scale(b, t)
abslon = glq.Abscissas(order)
abslon.scale(w, e)
abslat = glq.Abscissas(order)
abslat.scale(s, n)
abscissas = []
for lon in abslon._val:
    for lat in abslat._val:
        for r in absr._val:
            abscissas.append([lon,lat,r])
absdata = Gnuplot.PlotItems.Data(abscissas, with_="p pt 7 lt 7 ps 1", title=None)

gp.splot(absdata, \
        # Tesseroid
        *gpdata)
