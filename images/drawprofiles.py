#!/usr/bin/env python
"""
DrawProfiles:
    Draws profiles from the maps.
"""
################################################################################
# Created on 25-May-2009 15:41 PM
# Last modified by $Author: root $
__author__ = 'Leonardo Uieda (leouieda@gmail.com)'
__version__ = '$Revision: 51 $'
__date__ = '$Date: 2009-05-25 17:25:14 -0300 (Mon, 25 May 2009) $'
################################################################################

import os
import sys
import pylab as pl
import numpy as np
import math
import time
import logging
from matplotlib.font_manager import FontProperties

import glq
import tesseroid as ts
import tesseroidgravity as tg

from drawmaps import breaktess


def drawprofile(model, plon, plat, height, orders, folder):

    deg2rad = math.pi/180.0
    rad2deg = 180.0/math.pi

    # Create a logger for the prog
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # Create console handler and set level to debug
    ch = logging.StreamHandler(strm=sys.stderr)
    # Create formatter
    formatter = logging.Formatter("\n%(name)-20s %(levelname)-8s: %(message)s")
    # Add formatter to ch
    ch.setFormatter(formatter)
    ch.setLevel(logging.WARNING)
    # Add ch to logger
    logger.addHandler(ch)

    # Create the profiles
    w = -3
    e = 3
    dlon = (e-w)/100.0
    s = 42
    n = 48
    dlat = (n-s)/100.0
    lonsplon = np.arange(w, e + dlon, dlon)
    latsplon = [plat]*len(lonsplon)
    hsplon = [height]*len(lonsplon)
    latsplat = np.arange(s, n + dlat, dlat)
    lonsplat = [plon]*len(latsplat)
    hsplat = [height]*len(latsplat)

    names = []
    names.append("Gxx")
    names.append("Gxy")
    names.append("Gxz")
    names.append("Gyy")
    names.append("Gyz")
    names.append("Gzz")

    titles = []
    titles.append(r"$g_{xx}$")
    titles.append(r"$g_{xy}$")
    titles.append(r"$g_{xz}$")
    titles.append(r"$g_{yy}$")
    titles.append(r"$g_{yz}$")
    titles.append(r"$g_{zz}$")
    

    # Create a calculator class for the integrands of the GGT
    # Create the abscissas and weights
    abslon = glq.Abscissas(2)
    abslat = glq.Abscissas(2)
    absr = glq.Abscissas(2)
    wlon = glq.Weights(abslon)
    wlat = glq.Weights(abslat)
    wr = glq.Weights(absr)
    calculators = []
    tessgxx = tg.TesseroidGxx(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgxx)
    tessgxy = tg.TesseroidGxy(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgxy)
    tessgxz = tg.TesseroidGxz(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgxz)
    tessgyy = tg.TesseroidGyy(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgyy)
    tessgyz = tg.TesseroidGyz(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgyz)
    tessgzz = tg.TesseroidGzz(abslon, wlon, abslat, wlat, absr, wr)
    calculators.append(tessgzz)

    # The legend font
    font= FontProperties(size='small')

    print "Calculating LON profile at %f km" % (height/1000.0)
    # Create a figure
    fig = pl.figure(figsize=(7,10))
    #pl.figure()
    pl.suptitle(r"Perfil ao longo da latitude $%g^\circ$ a $%g$ km" \
                % (plat, height/1000.0), \
                fontsize=14)
    pl.subplots_adjust(hspace=0.5, wspace=0.31)

    # For each calculator, calculate the field and plot the map
    #index = [1,2,3,5,6,9]
    index = [1,2,3,4,5,6]
    totaltime = 0
    for calc, name, i, title in zip(*[calculators, names, index, titles]):
        print "  %s" % (name)
        subplt = pl.subplot(3, 2, i)
        pl.title(title, fontsize=18)
        pl.xlabel(r'Longitude ($ ^\circ$)', fontsize=10)

        for order in orders:
            print "    O: %d" % (order)

            # Create the abscissas and weights
            calc.ablon = glq.Abscissas(order)
            calc.ablat = glq.Abscissas(order)
            calc.abr = glq.Abscissas(order)
            calc.wlon = glq.Weights(calc.ablon)
            calc.wlat = glq.Weights(calc.ablat)
            calc.wr = glq.Weights(calc.abr)

            fieldfile = folder + os.path.sep + 'profile-lon-%s-h%g-o%d.txt' % \
                                                (name, height/1000.0, order)
            if os.path.exists(fieldfile):
            #if False:
                print "    Loading from file..."
                field = pl.load(fieldfile)
            else:
                # Calculate the field
                start = time.clock()
                field = calc.calculate(model, lonsplon, latsplon, hsplon)
                end = time.clock()
                deltat = end - start
                totaltime += deltat
                print "    Time it took (s): %f\n" % (deltat)
                # Save the data
                pl.save(fieldfile, field)

            pl.plot(lonsplon, field, label='%d'%(order))

        pl.legend(shadow=True, prop=font, loc='lower right')
        pl.xlim(w,e)

    filename = folder + os.path.sep + 'profile-lon-h%g.png' % \
                (height/1000.0)
    pl.savefig(filename, fmt='png')
    print "Total time (s): %f\n" % (totaltime)
    #pl.show()
    #pl.close()


    print "Calculating LAT profile at %f km" % (height/1000.0)
    # Create a figure
    fig = pl.figure(figsize=(7,10))
    #pl.figure()
    pl.suptitle(r"Perfil ao longo da longitude $%g^\circ$ a $%g$ km" \
                % (plon, height/1000.0), \
                fontsize=14)
    pl.subplots_adjust(hspace=0.5, wspace=0.31)

    # For each calculator, calculate the field and plot the map
    index = [1,2,3,4,5,6]
    totaltime = 0
    for calc, name, i, title in zip(*[calculators, names, index, titles]):
        print "  %s" % (name)
        subplt = pl.subplot(3, 2, i)
        pl.title(title, fontsize=18)
        pl.xlabel(r'Latitude ($^\circ$)', fontsize=10)

        for order in orders:
            print "    O: %d" % (order)

            # Create the abscissas and weights
            calc.ablon = glq.Abscissas(order)
            calc.ablat = glq.Abscissas(order)
            calc.abr = glq.Abscissas(order)
            calc.wlon = glq.Weights(calc.ablon)
            calc.wlat = glq.Weights(calc.ablat)
            calc.wr = glq.Weights(calc.abr)

            fieldfile = folder + os.path.sep + 'profile-lat-%s-h%g-o%d.txt' % \
                                                (name, height/1000.0, order)
            if os.path.exists(fieldfile):
            #if False:
                print "    Loading from file..."
                field = pl.load(fieldfile)
            else:
                # Calculate the field
                start = time.clock()
                field = calc.calculate(model, lonsplat, latsplat, hsplat)
                end = time.clock()
                deltat = end - start
                totaltime += deltat
                print "    Time it took (s): %f\n" % (deltat)
                # Save the data
                pl.save(fieldfile, field)

            pl.plot(latsplat, field, label='%d'%(order))

        pl.legend(shadow=True, prop=font, loc='lower right')
        pl.xlim(s,n)

    filename = folder + os.path.sep + 'profile-lat-h%g.png' % \
                (height/1000.0)
    pl.savefig(filename, fmt='png')
    print "Total time (s): %f\n" % (totaltime)
    #pl.show()
    pl.close()


if __name__ == '__main__':

    # Create the model
    tess = ts.Tesseroid(-0.5, 0.5, 44.5, 45.5, 0, 10000, 200, '1')
    nlon = 1
    nlat = 1
    model = breaktess(tess, nlon, nlat, 1)

    lon = 0.25
    lat = 45.2

    #height = (10)*1000
    #orders = [13,14,15,16,17,18,19,20,21,22]
    height = (250)*1000
    orders = [2,3,5,6,7,8,9]
    folder = 'model-%dx%d' % (nlon, nlat)

    # Check if folder exists
    if not os.path.exists(folder):
        print "  Creating dir"
        os.mkdir(folder)

    # Draw the map of the field
    drawprofile(model, lon, lat, height, orders, folder)