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

import math

__author__ = "Takayuki ITOH"


# open the SVG file
def openFile(filename):
    svgfile = open(filename, 'w')
    
    svgfile.write("<!DOCTYPE html>\n")
    #svgfile.write("<html>\n")
    #svgfile.write(" <body>\n")
    svgfile.write("  <svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n")

    return svgfile


# close the SVG file
def closeFile(svgfile):
    svgfile.write("  </svg>\n")
    #svgfile.write(" </body>\n")
    #svgfile.write("</html>\n")
    svgfile.close()


# draw a loof for a motif
def drawMotif(sizex, sizey, shiftx, shifty, rot, colorval, svgfile):
    deg = rot * 180.0 / math.pi
    if(colorval < 0.0):
        opacity = 0.5 + (colorval + 1.0) * 0.2
        if(opacity < 0.1):
            opacity = 0.1
        if(opacity > 1.0):
            opacity = 1.0
        b = -colorval * 50.0
        if(b > 192):
            b = 192
        if(b < 0):
            b = 0
        svgfile.write("    <ellipse stroke=\"#888\" fill=\"rgb(255,255," + str(int(b)) + ")\" fill-opacity=\"" + str(opacity) + "\" cx=\"" + str(shiftx) + "\" cy=\"" + str(shifty) + "\" rx=\"" + str(sizex) + "\" ry=\"" + str(sizey) + "\" transform=\"rotate(" + str(deg) + "," + str(shiftx) + "," + str(shifty) + ")\" />\n")
    else:
        g = int(colorval * 128)
        if(g < 0):
            g = 0
        if(g > 128):
            g = 128
        opacity = (1.25 - colorval) / 0.8
        if(opacity < 0.4):
            opacity = 0.4
        if(opacity > 1.0):
            opacity = 1.0
        svgfile.write("    <ellipse stroke=\"#888\" fill=\"rgb(255," + str(g) + ",0)\" fill-opacity=\"" + str(opacity) + "\" cx=\"" + str(shiftx) + "\" cy=\"" + str(shifty) + "\" rx=\"" + str(sizex) + "\" ry=\"" + str(sizey) + "\" transform=\"rotate(" + str(deg) + "," + str(shiftx) + "," + str(shifty) + ")\" />\n")

#annotate a motif
def annotateMotif(name, pvalue, x, y, svgfile):
    LINEHEIGHT = 13
    if(pvalue > 1.0):
        vword = '>1'
    else:
        vword = str.format('{0:.3g}', pvalue) 
#        vword = str.format('{0:.3e}', pvalue) 
    buf = "    <text x=\"" + str(x) + "\" y=\"" + str(y) + "\">" + name + "</text>\n"
    svgfile.write(buf)
    buf = "    <text x=\"" + str(x) + "\" y=\"" + str(y + LINEHEIGHT) + "\">" + vword + "</text>\n"
    svgfile.write(buf)
