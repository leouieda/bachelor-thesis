#!/usr/bin/env python
"""
MapIntegrand:
    Map how the integral kernel varies with respect to the integration point for
    a given calculation point.
"""
################################################################################
# Created on 7-May-2009 12:24 PM
# Last modified by $Author: root $
__author__ = 'Leonardo Uieda (leouieda@gmail.com)'
__version__ = '$Revision: 50 $'
__date__ = '$Date: 2009-05-22 15:32:11 -0300 (Fri, 22 May 2009) $'
################################################################################

import os
import sys
import pylab as pl
import numpy as np
import math

import glq
import tesseroid as ts
import tesseroidgravity as tg

def mapint(lon, lat, height, order, folder):
    """
    Map the integrand for the calculation point (lon,lat,height).
    """

    pl.figure(figsize=(15,10))
    pl.suptitle(r"Lon: $%g^\circ$   Lat: $%g^\circ$   Altitude: $%g\ km$" % \
                (lon, lat, height/1000.0), fontsize=16)

    # Convert lon, lat to radians
    deg2rad = math.pi/180.0
    rad2deg = 180.0/math.pi
    lon = deg2rad*lon
    lat = deg2rad*lat

    # Create the abscissas and weights
    abslon = glq.Abscissas(order)
    abslat = glq.Abscissas(order)
    absr = glq.Abscissas(order)
    wlon = glq.Weights(abslon)
    wlat = glq.Weights(abslat)
    wr = glq.Weights(absr)

    # Create a calculator class for the integrands of the GGT
    calculators = []
    names = []
    tessgxx = tg.TesseroidGxx(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgxx)
    names.append("Gxx")
    tessgxy = tg.TesseroidGxy(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgxy)
    names.append("Gxy")
    tessgxz = tg.TesseroidGxz(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgxz)
    names.append("Gxz")
    tessgyy = tg.TesseroidGyy(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgyy)
    names.append("Gyy")
    tessgyz = tg.TesseroidGyz(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgyz)
    names.append("Gyz")
    tessgzz = tg.TesseroidGzz(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgzz)
    names.append("Gzz")

    # Create a tesseroid
    tess = ts.Tesseroid(-0.5, 0.5, 44.5, 45.5, 0, 10000, 200, '1')
    r1 = tessgxx.R - tess['bottom']
    r2 = tessgxx.R - tess['top']
    r = tessgxx.R + height

    # For each calculator, map the integrand and plot
    index = [1,2,3,5,6,9]
    intlons = np.arange(tess['w'], tess['e']+0.01, 0.01)
    intlats = np.arange(tess['s'], tess['n']+0.01, 0.01)
    intlonsrad = deg2rad*intlons
    intlatsrad = deg2rad*intlats
    lonsgrid, latsgrid = pl.meshgrid(intlons, intlats)
    for calc, name, i in zip(*[calculators, names, index]):
        # Map the integrand
        integrand = []
        for latl in intlatsrad:
            tmp = []
            for lonl in intlonsrad:
                try:
                    tmp.append(calc.kernel(r, lon, lat, r1, r2, lonl, latl))
                except tg.SingularityError:
                    try:
                        tmp.append(tmp[-1])
                    except:
                        tmp.append(0)
            integrand.append(tmp)

        # Plot the integrand
        pl.subplot(3,3,i)
        #pl.title(name)
        #pl.xlabel('Longitude')
        #pl.ylabel('Latitude')
        pc = pl.pcolor(lonsgrid, latsgrid, integrand, cmap=pl.cm.jet)
        pl.colorbar(orientation='vertical', format='%2g')#,aspect=50)
        abslon.scale(tess['w'], tess['e'])
        abslat.scale(tess['s'], tess['n'])
        for alon in abslon._val:
            pl.plot([alon]*len(abslat._val), abslat._val, 'k+')
        pl.xlim(tess['w'], tess['e'])
        pl.ylim(tess['s'], tess['n'])

    filename = folder + os.path.sep + 'lon%g-lat%g-h%g-o%d.png' % \
                 (lon*rad2deg, lat*rad2deg, height/1000.0, order)
    pl.savefig(filename, fmt='png')
    pl.close()


if __name__ == '__main__':


    heights = [0, 1000, 10000, 20000, 30000, 50000, 100000, 150000, 250000]
    lons = [-20, -15, -10, -5, -4, -3, -2, -1.9, -1.8, -1.7, \
            -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, \
            -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, \
            0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, \
            1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 10, 15, 20]
    lats = [25, 30, 35, 40, 41, 42, 43, 43.1, 43.2, 43.3, 43.4, 43.5, \
            43.6, 43.7, 43.8, 43.9, 44, 44.1, 44.2, 44.3, 44.4, 44.5, \
            44.6, 44.7, 44.8, 44.9, 45, 45.1, 45.2, 45.3, 45.4, 45.5, 45.6, \
            45.7, 45.8, 45.9, 46, 46.1, 46.2, 46.3, 46.4, 46.5, 46.6, \
            46.7, 46.8, 46.9, 47, 48, 49, 50, 55, 60, 65]
    order = 5
    for height in heights:
        print "Height: %g km" % (height/1000.0)
        print "  Creating dir" 
        folder = 'height-%g' % (height/1000.0)
        # Check if folder exists
        if not os.path.exists(folder):
            os.mkdir(folder)

        # Draw the map of the field
        print "  Drawing maps"
        drawmaps(height, folder, lons, lats)

        #print "  Latitude profile"
        #subfolder = 'lat-profile'
        ## Check if it exists
        #if not os.path.exists(folder+os.path.sep+subfolder):
            #os.mkdir(folder+os.path.sep+subfolder)
        #lat = 45
        #for lon in lons:
            #print "    lat: %g  lon: %g" % (lat, lon)
            #mapint(lon, lat, height, order, folder+os.path.sep+subfolder)

        #print "  Longitude profile"
        #subfolder = 'lon-profile'
        ## Check if it exists
        #if not os.path.exists(folder+os.path.sep+subfolder):
            #os.mkdir(folder+os.path.sep+subfolder)
        #lon = 0
        #for lat in lats:
            #print "    lon: %g  lat: %g" % (lon, lat)
            #mapint(lon, lat, height, order, folder+os.path.sep+subfolder)
    