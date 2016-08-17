/*
 * Copyright (c) 2013, LAMP development team
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the LAMP development team nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAMP DEVELOPMENT TEAM BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Flower_svg.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <set>
#include <math.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace boost;

/**
 * draw a loof for a motif
 * @param sizex width of petal
 * @param sizey height of petal
 * @param shiftx shift x from the origin
 * @param shifty shift y from the origin
 * @param rot rotation
 * @param colorval color
 * @param svgfile output stream
 */
void Flower_svg::drawMotif(int sizex, int sizey, int shiftx, int shifty, double rot, double colorval, ofstream *svgfile) {
	double deg = rot * 180.0 / pi;
	double opacity;
	if (colorval < 0.0) {
		opacity = 0.5 + (colorval + 1.0) * 0.2;
		if (opacity < 0.1) {
			opacity = 0.1;
		}
		if (opacity > 1.0) {
			opacity = 1.0;
		}
		int b = int(-colorval * 50.0);
		if (b > 192) {
			b = 192;
		}
		if (b < 0) {
			b = 0;
		}
		*svgfile << "    <ellipse stroke=\"#888\" fill=\"rgb(255,255," << b << ")\" fill-opacity=\"" << opacity
				 << "\" cx=\"" << shiftx << "\" cy=\"" << shifty << "\" rx=\"" << sizex
				 << "\" ry=\"" << sizey << "\" transform=\"rotate(" << deg << "," << shiftx << "," << shifty << ")\" />\n" << endl;
	} else {
		int g = int(colorval * 128);
		if (g < 0) {
			g = 0;
		}
		if (g > 128) {
			g = 128;
		}
		opacity = (1.25 - colorval) / 0.8;
		if (opacity < 0.4) {
			opacity = 0.4;
		}
		if (opacity > 1.0) {
			opacity = 1.0;
		}
		*svgfile << "    <ellipse stroke=\"#888\" fill=\"rgb(255," << g << ",0)\" fill-opacity=\"" << opacity
				 << "\" cx=\"" << shiftx << "\" cy=\"" << shifty << "\" rx=\"" << sizex
				 << "\" ry=\"" << sizey << "\" transform=\"rotate(" << deg << "," << shiftx << "," << shifty << ")\" />\n" << endl;
	}

}

/**
 * annotate a motif
 * @param name motif name
 * @param pvalue p-value
 * @param x Location x
 * @param y Location y
 * @param svgfile output stream
 */
void Flower_svg::annotateMotif(string name, double pvalue, double x, double y, ofstream *svgfile) {
	int LINEHEIGHT = 13;
	string vword;
	if (pvalue > 1.0) {
		vword = ">1";
	} else {
		vword = boost::str(boost::format("%.4f") % pvalue);
	}
	*svgfile << "    <text x=\"" << x << "\" y=\"" << y << "\">" << name << "</text>\n" << endl;
	*svgfile << "    <text x=\"" << x << "\" y=\"" << (y + LINEHEIGHT) << "\">" << vword << "</text>\n" << endl;
}
