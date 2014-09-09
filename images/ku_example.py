# -*- coding: utf-8 -*-
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.ticker import MaxNLocator
import tesseroids_lib.tesseroidgravity as tg
import tesseroids_lib.tesseroid as ts

def main():

    lonO = 0
    latO = 65
    gdlon = 30
    gdlat = 30
    #bm = Basemap(projection='merc', \
                 #llcrnrlon=lonO-0.5*gdlon, llcrnrlat=latO-0.5*gdlat, \
                 #urcrnrlon=lonO+0.5*gdlon, urcrnrlat=latO+0.5*gdlat, \
                 #lon_0=lonO, lat_0=latO, resolution='c', area_thresh=10000)
    #bm = Basemap(projection='ortho', \
                 #lon_0=lonO, lat_0=latO, resolution='c', area_thresh=10000)
    bm = Basemap(projection='aeqd', \
                 llcrnrlon=lonO-0.5*gdlon, llcrnrlat=latO-0.5*gdlat, \
                 urcrnrlon=lonO+0.5*gdlon, urcrnrlat=latO+0.5*gdlat, \
                 lon_0=lonO, lat_0=latO, resolution='c', area_thresh=10000)
    #bm = Basemap(projection='laea', lat_ts=latO, \
                 #llcrnrlon=lonO-0.5*gdlon, llcrnrlat=latO-0.5*gdlat, \
                 #urcrnrlon=lonO+0.5*gdlon, urcrnrlat=latO+0.5*gdlat, \
                 #lon_0=lonO, lat_0=latO, resolution='c', area_thresh=10000)

    dlon = 5
    dlat = 5
    dz = 30000
    tess = ts.Tesseroid(lonO-0.5*delon, lonO+0.5*delon, \
                        latO-0.5*delat, latO+0.5*delat, \
                        0, dz, 1000)
    model = [tess]

    # Make a grid
    gw = lonO-0.5*gdlon
    ge = lonO+0.5*gdlon
    steplon = (ge - gw)/100.
    gs = latO-0.5*gdlat
    gn = latO+0.5*gdlat
    steplat = (gn - gs)/100.
    print "Gerando grid: %g/%g/%g/%g" % (gw, ge, gs, gn)
    globalGridLons = np.arange(gw, ge+steplon, steplon)
    globalGridLats = np.arange(gs, gn+steplat, steplat)
    h = 250000.
    # Cria grids com os lons e lats
    glons, glats = bm(*pl.meshgrid(globalGridLons, globalGridLats))


if __name__ == "__main__":
    main()