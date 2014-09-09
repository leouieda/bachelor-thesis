"""
This is a script for plotting a sketch of the geocentric and local coordinate
systems unp.sing gnuplot.
"""

import Gnuplot
import numpy as np
from math import cos, sin, pi

d2r = pi/180.0


# The ROTAION and REFLECTION matrices used
# ##############################################################################
def R1(angle):
    """
    Returns in a matrix object the rotation matrix around the x axis.
    angle should be in degrees!
    """
    global d2r
    return np.matrix([[1.,0.,0.],[0.,np.cos(d2r*angle),np.sin(d2r*angle)],[0.,-np.sin(d2r*angle),np.cos(d2r*angle)]])

def R2(angle):
    """
    Returns in a matrix object the rotation matrix around the y axis.
    angle should be in degrees!
    """
    global d2r
    return np.matrix([[np.cos(d2r*angle),0.,-np.sin(d2r*angle)],[0.,1.,0.],[np.sin(d2r*angle),0.,np.cos(d2r*angle)]])

def R3(angle):
    """
    Returns in a matrix object the rotation matrix around the z axis.
    angle should be in degrees!
    """
    global d2r
    return np.matrix([[np.cos(d2r*angle),np.sin(d2r*angle),0.],\
                      [-np.sin(d2r*angle),np.cos(d2r*angle),0.],\
                      [0.,0.,1.]])

def P2():
    """
    Returns in a matrix object the reflection matrix of the y axis.
    """
    return np.matrix([[1.,0.,0.],[0.,-1.,0.],[0.,0.,1.]])
# ##############################################################################


def tessSideGCirc(s, n, lon, r, sent=1):
    res = []
    step = sent*0.2
    for lat in np.arange(s, n+step, step):
        res.append([lon, lat, r])
    return res

def tessSidePar(w, e, lat, r, sent=1):
    res = []
    step = sent*0.2
    for lon in np.arange(w, e+step, step):
        res.append([lon, lat, r])
    return res


def main():

    gp = Gnuplot.Gnuplot(persist = 1)

    output = 'tesseroid_coordsys_raw.ps'
    gp("set term postscript")
    gp('set output "%s"' % (output))
    gp("set mapping spherical")
    gp("set angle degrees")
    #gp("set xrange [0:1]; set yrange [0:1]; set zrange [0:1]")
    #gp('set format x ""')
    #gp('set format y ""')
    #gp('set format z ""')
    gp('unset xtics; unset ytics; unset ztics')
    gp('unset border')
    gp('set size ratio 1')
    gp('set xrange [-0.2:1.5];set yrange [-0.2:1.5];set zrange [-0.2:1.5]')
    gp("set ticslevel 0")
    #gp('set xlabel "x"; set ylabel "y"; set zlabel "z"')
    gp("set view 80,145")
    gp("set parametric")
    gp("set dummy u,v")
    gp("set isosample 6")
    gp("set urange [0:90]; set vrange [0:90]")

    # Make the plot items
    ############################################################################

    # The geocentric axis
    x_gc = Gnuplot.PlotItems.Data([[0,0,0, 1.2,0,0]], with_="vectors lt 1 lw 3", title=None)
    y_gc = Gnuplot.PlotItems.Data([[0,0,0, 0,1.2,0]], with_="vectors lt 1 lw 3", title=None)
    z_gc = Gnuplot.PlotItems.Data([[0,0,0, 0,0,1.2]], with_="vectors lt 1 lw 3", title=None)


    # The computation point
    lat = 36
    lon = 72
    r = 1
    comput_point = Gnuplot.PlotItems.Data([[lon,lat,r]], \
                                        with_="p pt 7 lt 7 ps 1", title=None)
    lines = [[lon,lat,r], [lon,0,r*cos(d2r*lat)], [0,0,0], [lon,lat,r]]
    cp_proj_lines = Gnuplot.PlotItems.Data(lines , with_="l lt 3 lw 1.5", title=None)

    # The local axis
    xlocal = np.matrix([0.3,0,0]).T
    ylocal = np.matrix([0,0.3,0]).T
    zlocal = np.matrix([0,0,0.4]).T
    # transform them to the geocentric system
    xlocal_gc = R3(180.0 - lon)*R2(90.0 - lat)*P2()*xlocal
    xlocal_gc = xlocal_gc.T.tolist()[0] # convert it to a list
    xlocal_vector_data = [[lon,lat,r, \
                          xlocal_gc[0], xlocal_gc[1], xlocal_gc[2]]]

    ylocal_gc = R3(180.0 - lon)*R2(90.0 - lat)*P2()*ylocal
    ylocal_gc = ylocal_gc.T.tolist()[0] # convert it to a list
    ylocal_vector_data = [[lon,lat,r, \
                          ylocal_gc[0], ylocal_gc[1], ylocal_gc[2]]]

    zlocal_gc = R3(180.0 - lon)*R2(90.0 - lat)*P2()*zlocal
    zlocal_gc = zlocal_gc.T.tolist()[0] # convert it to a list
    zlocal_vector_data = [[lon,lat,r, \
                          zlocal_gc[0], zlocal_gc[1], zlocal_gc[2]]]

    x_l_P = Gnuplot.PlotItems.Data(xlocal_vector_data, with_="vectors lt 1 lw 3", title=None)
    y_l_P = Gnuplot.PlotItems.Data(ylocal_vector_data, with_="vectors lt 1 lw 3", title=None)
    z_l_P = Gnuplot.PlotItems.Data(zlocal_vector_data, with_="vectors lt 1 lw 3", title=None)

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

    # Integration point Q
    lat_Q = (n+s)/2
    lon_Q = (e+w)/2
    r_Q = (t+b)/2
    int_point = Gnuplot.PlotItems.Data([[lon_Q,lat_Q,r_Q]], \
                                        with_="p pt 7 lt 7 ps 1", title=None)
    l_line = Gnuplot.PlotItems.Data([[lon_Q,lat_Q,r_Q],[lon,lat,r]] , with_="l lt 6 lw 1.5", title=None)
    ################################################################################

    # Make the plots
    #lines = [[w,0,0],[w,0,1]] + tessSideGCirc(0, 90, w, 1, 1) + \
                #tessSideGCirc(90, 0, e, 1, -1) + [[e,0,1],[e,0,0]]
    #guidelineslon = Gnuplot.PlotItems.Data(lines, with_="l lt 4 lw 1", title=None)
    #guidelineslat = Gnuplot.PlotItems.Data(tessSidePar(0, 90, 0, 1), with_="l lt 4 lw 1", title=None)
    gridlines = "cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lt 4 lw 1 notitle"

    # Plot the geocentric system
    gp.splot(gridlines, \
             #guidelineslon, \
             #guidelineslat, \
             x_gc, y_gc, z_gc, \
    # The computation point P
             comput_point,\
             cp_proj_lines,\
    # The local system P
             x_l_P, y_l_P, z_l_P, \
    # Tesseroid
             int_point, l_line, \
    # Ponto de integracao Q
             *gpdata)


if __name__ == '__main__':
    main()