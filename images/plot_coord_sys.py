"""
This is a script for plotting a sketch of the geocentric and local coordinate
systems unp.sing gnuplot.
"""

import Gnuplot
import numpy as np
from math import cos, sin, pi, asin, atan2, sqrt

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

    output = 'coordinate_sys_raw.ps'
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
    gp('set xrange [-0.2:2];set yrange [-0.2:2];set zrange [-0.2:2]')
    gp("set ticslevel 0")
    #gp('set xlabel "x"; set ylabel "y"; set zlabel "z"')
    gp("set view 70,140")
    gp("set parametric")
    gp("set dummy u,v")
    gp("set isosample 5")
    gp("set urange [0:90]; set vrange [0:90]")

    # Make the plot items
    ############################################################################

    global d2r
    r2d = 180.0/pi
    
    # The geocentric axis
    x_gc = Gnuplot.PlotItems.Data([[0,0,0, 1.5,0,0]], with_="vectors lt 1 lw 3", title=None)
    y_gc = Gnuplot.PlotItems.Data([[0,0,0, 0,1.5,0]], with_="vectors lt 1 lw 3", title=None)
    z_gc = Gnuplot.PlotItems.Data([[0,0,0, 0,0,1.5]], with_="vectors lt 1 lw 3", title=None)
    

    # The computation point
    lat = 45
    lon = 75
    r = 2
    comput_point = Gnuplot.PlotItems.Data([[lon,lat,r]], \
                                        with_="p pt 7 lt 7 ps 1", title=None)
    #r_vector_data = [[0, 0, 0,\
                #r*cos(d2r*lat)*cos(d2r*lon), r*cos(d2r*lat)*sin(d2r*lon), r*sin(d2r*lat)]]
    #r_vector = Gnuplot.PlotItems.Data(r_vector_data, with_="vectors lt 7", title=None)
    lines = [[lon,lat,r], [lon,0,r*cos(d2r*lat)], [0,0,0], [lon,lat,r]]
    cp_proj_lines = Gnuplot.PlotItems.Data(lines , with_="l lt 3 lw 1.5", title=None)

    # The local axis
    xlocal = np.matrix([0.3,0,0]).T
    ylocal = np.matrix([0,0.3,0]).T
    zlocal = np.matrix([0,0,0.5]).T
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

    x_l_P = Gnuplot.PlotItems.Data(xlocal_vector_data, with_="vectors lt 1 lw 2", title=None)
    y_l_P = Gnuplot.PlotItems.Data(ylocal_vector_data, with_="vectors lt 1 lw 2", title=None)
    z_l_P = Gnuplot.PlotItems.Data(zlocal_vector_data, with_="vectors lt 1 lw 2", title=None)

    # Prisma
    lon_Q = 15
    lat_Q = 45
    r_Q = 2
    w = -0.15
    e = 0.15
    s = -0.15
    n = 0.15
    t = 0
    b = -0.5
    sides1_l = [[[s,w,t],[n,w,t]],[[n,w,t],[n,e,t]],[[n,e,t],[s,e,t]],[[s,e,t],[s,w,t]], \
                [[s,w,b],[s,e,b]],[[s,e,b],[n,e,b]], \
                [[s,w,t],[s,w,b]],[[s,e,t],[s,e,b]],[[n,e,t],[n,e,b]]]
    sides2_l = [[[n,w,t],[n,w,b]],[[n,w,b],[n,e,b]],[[n,w,b],[s,w,b]]]

    # converte the sides to global system
    e_prisma = [cos(d2r*lat_Q)*cos(d2r*lon_Q), \
                cos(d2r*lat_Q)*sin(d2r*lon_Q), \
                sin(d2r*lat_Q)]
    e_prisma = r_Q*np.matrix(e_prisma).T
    sides1 = []
    for side in sides1_l:
        side_g = []
        for point in side:
            point = np.matrix(point).T
            point_g = R3(180.0 - lon_Q)*R2(90.0 - lat_Q)*P2()*point + e_prisma
            point_g = point_g.T.tolist()[0] # convert it to a list
            # Converte pra esfericas            
            xG = point_g[0]
            yG = point_g[1]
            zG = point_g[2]
            r = sqrt(xG**2 + yG**2 + zG**2)
            lat = r2d*asin(zG/r)
            lon = r2d*atan2(yG, xG)
            point_g = [lon,lat,r]
            side_g.append(point_g)
        sides1.append(side_g)
    sides2 = []
    for side in sides2_l:
        side_g = []
        for point in side:
            point = np.matrix(point).T
            point_g = R3(180.0 - lon_Q)*R2(90.0 - lat_Q)*P2()*point + e_prisma
            point_g = point_g.T.tolist()[0] # convert it to a list
            # Converte pra esfericas
            xG = point_g[0]
            yG = point_g[1]
            zG = point_g[2]
            r = sqrt(xG**2 + yG**2 + zG**2)
            lat = r2d*asin(zG/r)
            lon = r2d*atan2(yG, xG)
            point_g = [lon,lat,r]
            side_g.append(point_g)
        sides2.append(side_g)

        
    gpdata = []
    for line in sides1:
        gpdata.append(Gnuplot.PlotItems.Data(line, with_="l lt 1 lw 1.5", title=None))
    for line in sides2:
        gpdata.append(Gnuplot.PlotItems.Data(line, with_="l lt 2 lw 1.5", title=None))

    # Integration point Q
    int_point = Gnuplot.PlotItems.Data([[lon_Q,lat_Q,r_Q]], \
                                        with_="p pt 7 lt 7 ps 1", title=None)
    lines = [[lon_Q,lat_Q,r_Q], [lon_Q,0,r_Q*cos(d2r*lat_Q)], [0,0,0], [lon_Q,lat_Q,r_Q]]
    ip_proj_lines = Gnuplot.PlotItems.Data(lines , with_="l lt 3 lw 1.5", title=None)
    #l_line = Gnuplot.PlotItems.Data([[lon_Q,lat_Q,r_Q],[lon,lat,r]] , with_="l lt 6 lw 2", title=None)
    
    # The local axis
    xlocal = np.matrix([0.3,0,0]).T
    ylocal = np.matrix([0,0.3,0]).T
    zlocal = np.matrix([0,0,0.4]).T
    # transform them to the geocentric system
    xlocal_gc = R3(180.0 - lon_Q)*R2(90.0 - lat_Q)*P2()*xlocal
    xlocal_gc = xlocal_gc.T.tolist()[0] # convert it to a list
    xlocal_vector_data = [[lon_Q,lat_Q,r_Q, \
                          xlocal_gc[0], xlocal_gc[1], xlocal_gc[2]]]

    ylocal_gc = R3(180.0 - lon_Q)*R2(90.0 - lat_Q)*P2()*ylocal
    ylocal_gc = ylocal_gc.T.tolist()[0] # convert it to a list
    ylocal_vector_data = [[lon_Q,lat_Q,r_Q, \
                          ylocal_gc[0], ylocal_gc[1], ylocal_gc[2]]]

    zlocal_gc = R3(180.0 - lon_Q)*R2(90.0 - lat_Q)*P2()*zlocal
    zlocal_gc = zlocal_gc.T.tolist()[0] # convert it to a list
    zlocal_vector_data = [[lon_Q,lat_Q,r_Q, \
                          zlocal_gc[0], zlocal_gc[1], zlocal_gc[2]]]

    x_l_Q = Gnuplot.PlotItems.Data(xlocal_vector_data, with_="vectors lt 1 lw 2", title=None)
    y_l_Q = Gnuplot.PlotItems.Data(ylocal_vector_data, with_="vectors lt 1 lw 2", title=None)
    z_l_Q = Gnuplot.PlotItems.Data(zlocal_vector_data, with_="vectors lt 1 lw 2", title=None)
    ################################################################################

    # Make the plots
    #gridlines = "cos(u)*cos(v),cos(u)*sin(v),sin(u) w l lt 4 lw 1 notitle"
    lines = [[w,0,0],[w,0,1]] + tessSideGCirc(0, 90, w, 1, 1) + \
                tessSideGCirc(90, 0, e, 1, -1) + [[e,0,1],[e,0,0]]
    guidelineslon = Gnuplot.PlotItems.Data(lines, with_="l lt 4 lw 1", title=None)
    
    # Plot the geocentric system
    gp.splot(#guidelineslon, \
             x_gc, y_gc, z_gc, \
    # The computation point P
             comput_point,\
             cp_proj_lines,\
    # The local system P
             x_l_P, y_l_P, z_l_P, \
    # The local system Q
             x_l_Q, y_l_Q, z_l_Q, \
    # Ponto de integracao Q
             int_point, \
             ip_proj_lines, \
    # Prisma
             *gpdata)
    

if __name__ == '__main__':
    main()