#!/usr/bin/env python

"""
Copyright (c) 2013, LAMP development team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the LAMP development team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import sys, math, optparse
import flower.flower_svg as svg
import flower.flower_readfile as readfile
from flower.flower_readfile import *
from optparse import OptionParser

__author__ = "Takayuki ITOH"


# definition of main function
def main():

    if len(sys.argv) < 2:
        print 'Usage: # python %s resultfilename (target(csv)filename) (value(tab)filename)' % sys.argv[0]
        quit()

    # default values of sizes of petals
    EWIDTH = 120
    EHEIGHT = 50
    MINWIDTH = 120
    MINHEIGHT = 40

    # shift from the origin(=upper-left end)
    SHIFTX = 400
    SHIFTY = 400

    # coeffcient for the adjustment of petal position
    PETALPOSCOEF = 0.8

    # coefficient for the adjustment of torus size
    TORUSSIZECOEF = 1.2

    # coefficients for size calculation
    sizecoef = 0.3

    # add command line options
    parser = OptionParser()
    parser.add_option("--ewidth", "--EWIDTH", dest="EWIDTH", default="120")
    parser.add_option("--eheight", "--EHEIGHT", dest="EHEIGHT", default="50")
    parser.add_option("--minwidth", "--MINWIDTH", dest="MINWIDTH", default="120")
    parser.add_option("--minheight", "--MINHEIGHT", dest="MINHEIGHT", default="40")
    parser.add_option("--shiftx", "--SHIFTX", dest="SHIFTX", default="400")
    parser.add_option("--shifty", "--SHIFTY", dest="SHIFTY", default="400")
    parser.add_option("--petalposcoef", "--PETALPOSCOEF", dest="PETALPOSCOEF", default="0.8")
    parser.add_option("--torussizecoef", "--TORUSSIZECOEF", dest="TORUSSIZECOEF", default="1.2")
    (options, args) = parser.parse_args()
    EWIDTH = int(options.EWIDTH)
    EHEIGHT = int(options.EHEIGHT)
    MINWIDTH = int(options.MINWIDTH)
    MINHEIGHT = int(options.MINHEIGHT)
    SHIFTX = int(options.SHIFTX)
    SHIFTY = int(options.SHIFTY)
    PETALPOSCOEF = float(options.PETALPOSCOEF)
    TORUSSIZECOEF = float(options.TORUSSIZECOEF)

    # read the result file
    csvfname = ''
    tabfname = ''
    if len(sys.argv) > 3:
        cfvfname = sys.argv[2]
        tabfname = sys.argv[3]
    significance = readfile.readResult(sys.argv[1], csvfname, tabfname)
    colorvalMin = -math.log(1.0 / significance)

    # for each combination
    for j in range(len(combiName)):

        # reset the end position values
        xmin = 1.0e+30
        xmax = -1.0e+30
        ymin = 1.0e+30
        ymax = -1.0e+30

        # determine the scale of drawing
        nameset = combiName[j]
        for k in range(len(nameset)):
            
            # search for the motif and specify its ID
            nameid = -1
            for i in range(len(motifName)):
                if(motifName[i].find(nameset[k]) >= 0):
                    nameid = i
                    break

            # calculate the end position values
            sizeval = math.log(motifNgenes[nameid]) * sizecoef
            rad = k * math.pi * 2 / len(combiName[j]) - math.pi * 0.5
            sizex = EWIDTH * sizeval
            sizey = EHEIGHT * sizeval
            shiftx = (math.cos(rad) * (sizex - EHEIGHT * 2)) + SHIFTX
            shifty = (math.sin(rad) * (sizex - EHEIGHT * 2)) + SHIFTY
            if(xmax < shiftx + sizex):
                xmax = shiftx + sizex
            if(xmin > shiftx - sizex):
                xmin = shiftx - sizex
            if(ymax < shifty + sizey):
                ymax = shifty + sizey
            if(ymin > shifty - sizey):
                ymin = shifty - sizey

        # calculate the scaling factor
        scalex = (xmax - xmin) / (SHIFTX * 3)
        scaley = (ymax - ymin) / (SHIFTY * 3)

        # open the SVG file
        nameset = combiName[j]
        svgfilename = sys.argv[1] + '-flower'
        #for k in range(len(nameset)):
        #    svgfilename += ('-' + nameset[k])
        svgfilename += str(combiRank[j])
        svgfilename += '.svg'
        svgfile = svg.openFile(svgfilename)

        # for each motif in the combination
        for k in range(len(nameset)):
            
            # search for the motif and specify its ID
            nameid = -1
            for i in range(len(motifName)):
                if(motifName[i].find(nameset[k]) >= 0):
                    nameid = i
                    break

            # calculate the size and color of the ellipsoid
            sizeval = math.log(motifNgenes[nameid]) * sizecoef
            if(motifApvalue[nameid] > significance):
                colorval = -math.log(motifApvalue[nameid] / significance)
                if(colorval < colorvalMin):
                    colorval = colorvalMin
            if(motifApvalue[nameid] <= significance):
                colorval = motifApvalue[nameid] / significance

            # initialize variables before calculating the shape of the ellipsoid
            rad = k * math.pi * 2 / len(combiName[j]) - math.pi * 0.5
            sizex = EWIDTH * sizeval * scalex
            sizey = EHEIGHT * sizeval * scaley
            if(sizex < MINWIDTH):
                sizex = MINWIDTH
            if(sizey < MINHEIGHT):
                sizey = MINHEIGHT
            shiftx = (math.cos(rad) * (sizex - EHEIGHT * PETALPOSCOEF)) + SHIFTX
            shifty = (math.sin(rad) * (sizex - EHEIGHT * PETALPOSCOEF)) + SHIFTY
            annox = (math.cos(rad) * (sizex * 2  - EHEIGHT * PETALPOSCOEF)) + SHIFTX
            annoy = (math.sin(rad) * (sizex * 2  - EHEIGHT * PETALPOSCOEF)) + SHIFTY

            # draw a motif
            svg.drawMotif(sizex, sizey, shiftx, shifty, rad, colorval, svgfile)
            svg.annotateMotif(motifName[nameid], motifApvalue[nameid],
                              annox, annoy, svgfile)

        # draw the combination
        colorval = combiApvalue[j] / significance
        size = EHEIGHT * scalex * TORUSSIZECOEF
        if(size < MINHEIGHT):
            size = MINHEIGHT
        svg.drawMotif(size, size, SHIFTX, SHIFTY, 0.0, colorval, svgfile)
        svg.annotateMotif("", combiApvalue[j], SHIFTX-30, SHIFTY-5, svgfile)

        # close the SVG file
        svg.closeFile(svgfile)

    
# call the main function
if __name__ == '__main__':
    main()

