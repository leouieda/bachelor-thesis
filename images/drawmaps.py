#!/usr/bin/env python
"""
DrawMaps:
    Make some sample maps at varying heights and GLQ orders to see when it fails
"""
################################################################################
# Created on 12-May-2009 15:41 PM
# Last modified by $Author: leo $
__author__ = 'Leonardo Uieda (leouieda@gmail.com)'
__version__ = '$Revision: 55 $'
__date__ = '$Date: 2009-07-01 22:04:03 -0300 (Wed, 01 Jul 2009) $'
################################################################################

import os
import sys
import pylab as pl
import numpy as np
import math
import time
import logging

import glq
import tesseroid as ts
import tesseroidgravity as tg

def breaktess(tess, nlon, nlat, nr):
    """
    Break the tesseroid into nlon x nlat x nr pieces.
    """
    dlon = (tess['e'] - tess['w'])/(nlon)
    dlat = (tess['n'] - tess['s'])/(nlat)
    dr = (tess['bottom'] - tess['top'])/(nr)

    ws = np.arange(tess['w'], tess['e'], dlon)
    es = np.arange(tess['w']+dlon, tess['e']+dlon, dlon)
    ss = np.arange(tess['s'], tess['n'], dlat)
    ns = np.arange(tess['s']+dlat, tess['n']+dlat, dlat)
    tops = np.arange(tess['top'], tess['bottom'], dr)
    bottoms = np.arange(tess['top']+dr, tess['bottom']+dr, dr)

    density = tess['density']
    model = []
    i = 1
    for top, bottom in zip(*[tops, bottoms]):
        for w, e in zip(*[ws, es]):
            for s, n in zip(*[ss, ns]):
                tess = ts.Tesseroid(w, e, s, n, top, bottom, density, '%d' % (i))
                model.append(tess)
                i += 1

    return model


def drawmaps(model, lons, lats, height, order, folder):
    """
    Draw the tensor maps at height and save them in folder.
    """

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
    ch.setLevel(logging.DEBUG)
    # Add ch to logger
    logger.addHandler(ch)

    # Create a calculation grid
    lonslist = []
    latslist = []
    heights = []
    # Iterate over the lats
    for lat in lats:
        # Iterate over the lons
        for lon in lons:
            lonslist.append(lon)
            latslist.append(lat)
            heights.append(height)
    # Mesh a grid with lons, lats
    glons, glats = pl.meshgrid(lons, lats)

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

    # Create the abscissas and weights
    abslon = glq.Abscissas(order)
    abslat = glq.Abscissas(order)
    absr = glq.Abscissas(order)
    wlon = glq.Weights(abslon)
    wlat = glq.Weights(abslat)
    wr = glq.Weights(absr)
      
    # Create a calculator class for the integrands of the GGT
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

    # Calculate the average distances between nodes
    print "\nCalculating distance between nodes:"
    for tess in model:
    #for tess in [model[0]]:
        print "  Tesseroid %s:" % (tess['tag'])
        abslon.scale(tess['w'], tess['e'])
        abslat.scale(tess['s'], tess['n'])
        R = tessgxx.R - (tess['bottom'] - tess['top'])/2.0
        londists = []
        latdists = []
        for i in range(0, len(abslon)-1):
            londists.append(abs(R*deg2rad*abslon[i+1] - \
                            R*deg2rad*abslon[i]))
            latdists.append(abs(R*deg2rad*abslat[i+1] - \
                            R*deg2rad*abslat[i]))
        londists = np.array(londists)
        latdists = np.array(latdists)
        print "    Longitude nodes:"
        print "      Min: %f km" % (londists.min()/1000.0)
        print "      Max: %f km" % (londists.max()/1000.0)
        print "      Mean: %f km" % (londists.mean()/1000.0)
        print "    Latitude nodes:"
        print "      Min: %f km" % (latdists.min()/1000.0)
        print "      Max: %f km" % (latdists.max()/1000.0)
        print "      Mean: %f km\n" % (latdists.mean()/1000.0)

    print "Calculating GGT at %f km" % (height/1000.0)
    # Create a figure
    pl.figure(figsize=(10,6))
    #pl.figure()
    pl.suptitle("Gradiente da Gravidade a %g km (Ordem:%d)" \
                % (height/1000.0, order), \
                fontsize=14)
    pl.subplots_adjust(hspace=0.05, wspace=0.31)

    # For each calculator, calculate the field and plot the map
    #index = [1,2,3,5,6,9]
    index = [1,2,3,4,5,6]
    ggt = {}
    totaltime = 0
    for calc, name, i, title in zip(*[calculators, names, index, titles]):
        print "  %s" % (name)

        fieldfile = folder + os.path.sep + '%s-h%g-o%d.txt' % \
                                            (name, height/1000.0, order)
        if os.path.exists(fieldfile):
        #if False:
            print "    Loading from file..."
            fieldgrd = pl.load(fieldfile)
        else:
            # Calculate the field
            start = time.clock()
            field = calc.calculate(model, lonslist, latslist, heights)
            end = time.clock()
            deltat = end - start
            totaltime += deltat
            print "    Time it took (s): %f\n" % (deltat)

            # Convert the field to a matrix
            fieldlist = []
            for j in range(len(lons),len(field)+1,len(lons)):
                fieldlist.append(field[j-len(lons):j])
            # Make a Numpy array
            fieldgrd = np.array(fieldlist)

            # Save the data
            pl.save(fieldfile, fieldgrd)
        ggt[name] = fieldgrd
        # Plot the integrand
        #pl.subplot(3,3,i)
        pl.subplot(2, 3, i, aspect='equal')
        pl.title(title, fontsize=18)
        #pl.xlabel('Longitude')
        #pl.ylabel('Latitude')
        pc = pl.pcolor(glons, glats, fieldgrd, cmap=pl.cm.jet)
        pl.colorbar(orientation='vertical', format='%g', shrink=0.73)

        ## Plot the tesseroids
        #for tess in model:
            #tesslons = [tess['w'], tess['e'], tess['e'], tess['w'], tess['w']]
            #tesslats = [tess['s'], tess['s'], tess['n'], tess['n'], tess['s']]
            #pl.plot(tesslons, tesslats, 'w-', linewidth=1)
            ## Plot the abscissas
            #abslon.scale(tess['w'], tess['e'])
            #abslat.scale(tess['s'], tess['n'])
            #for alon in abslon._val:
                #pl.plot([alon]*len(abslat._val), abslat._val, 'k+')

        pl.xlim(lons[0], lons[-1])
        pl.ylim(lats[0], lats[-1])

    filename = folder + os.path.sep + 'ggt-h%g-o%d.png' % \
                (height/1000.0, order)
    pl.savefig(filename, fmt='png')
    print "Total time (s): %f\n" % (totaltime)

    #print "Calculating invariants...\n"
    #pl.figure()
    #pl.suptitle("Invariantes a %g km (Ordem:%d)" \
                #% (height/1000.0, order))
    #pl.subplots_adjust(hspace=0.3, wspace=0.31)

    ## I0
    #pl.subplot(2, 2, 1, aspect='equal')
    #pl.title(r"$I_0$")
    #trace = ggt['Gxx'] + ggt['Gyy'] + ggt['Gzz']
    #pl.pcolor(glons, glats, trace, cmap=pl.cm.jet)
    #pl.colorbar(orientation='vertical', format='%g')
    #pl.xlim(lons[0], lons[-1])
    #pl.ylim(lats[0], lats[-1])

    ## I1   
    #pl.subplot(2, 2, 2, aspect='equal')
    #pl.title(r"$I_1$")
    #I1 = ggt['Gxx']*ggt['Gyy'] + ggt['Gyy']*ggt['Gzz'] + \
         #ggt['Gxx']*ggt['Gzz'] - ggt['Gxy']**2 - \
         #ggt['Gyz']**2 - ggt['Gxz']**2
    #pl.pcolor(glons, glats, I1, cmap=pl.cm.jet)
    #pl.colorbar(orientation='vertical', format='%g')
    #pl.xlim(lons[0], lons[-1])
    #pl.ylim(lats[0], lats[-1])

    ## I2
    #pl.subplot(2, 2, 3, aspect='equal')
    #pl.title(r"$I_2$")
    #I2 = ggt['Gxx']*(ggt['Gyy']*ggt['Gzz'] - ggt['Gyz']**2) + \
         #ggt['Gxy']*(ggt['Gyz']*ggt['Gxz'] - ggt['Gxy']*ggt['Gzz']) + \
         #ggt['Gxz']*(ggt['Gxy']*ggt['Gyz'] - ggt['Gxz']*ggt['Gyy'])
    #pl.pcolor(glons, glats, I2, cmap=pl.cm.jet)
    #pl.colorbar(orientation='vertical', format='%g')
    #pl.xlim(lons[0], lons[-1])
    #pl.ylim(lats[0], lats[-1])

    ## I
    #pl.subplot(2, 2, 4, aspect='equal')
    #pl.title(r"$I$")
    #I = -((I2/2.0)**2)/((I1/3.0)**3)
    #pl.pcolor(glons, glats, I, cmap=pl.cm.jet, vmin=0, vmax=1)
    #pl.colorbar(orientation='vertical', format='%g')
    #pl.xlim(lons[0], lons[-1])
    #pl.ylim(lats[0], lats[-1])

    #filename = folder + os.path.sep + 'invariants-h%g-o%d.png' % \
                #(height/1000.0, order)
    #pl.savefig(filename, fmt='png')
    #pl.show()
    pl.close()

if __name__ == '__main__':

    # Create the model
    tess = ts.Tesseroid(-0.5, 0.5, 44.5, 45.5, 0, 10000, 200, '1')
    nlon = int(sys.argv[3])
    nlat = int(sys.argv[4])
    model = breaktess(tess, nlon, nlat, 1)

    # Computaion grid
    w = -3
    e = 3
    dlon = (e-w)/100.0
    s = 42
    n = 48
    dlat = (n-s)/100.0
    lons = np.arange(w, e + dlon, dlon)
    lats = np.arange(s, n + dlat, dlat)

    height = float(sys.argv[1])*1000
    order = int(sys.argv[2])
    folder = 'model-%dx%d' % (nlon, nlat)

    # Check if folder exists
    if not os.path.exists(folder):
        print "  Creating dir"
        os.mkdir(folder)

    # Draw the map of the field
    drawmaps(model, lons, lats, height, order, folder)