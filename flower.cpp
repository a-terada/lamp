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

#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>
#include <boost/program_options.hpp>
#include <sstream>
#include <fstream>
#include <list>
#include "flowers/Flower_readfile.h"
#include "flowers/Flower_svg.h"

using namespace std;
namespace po = boost::program_options;
using namespace boost;

/** @file */

/**
 * Draw flower diagram
 * @param argc
 * @param argv
 * @return 
 */
int main(int argc, char** argv) {
	if (argc < 2) {
		string mssg = "Usage: flower resultfilename (target(csv)filename) (value(tab)filename)";
		cout << mssg << endl;
		return 0;
	}

	// default values of sizes of petals
	int EWIDTH = 120;
	int EHEIGHT = 50;
	int MINWIDTH = 120;
	int MINHEIGHT = 40;
	double sizeval;
	double rad;
	int sizex;
	int sizey;
	int shiftx;
	int shifty;

	// shift from the origin(=upper-left end)
	int SHIFTX = 400;
	int SHIFTY = 400;

	// coeffcient for the adjustment of petal position
	double PETALPOSCOEF = 0.8;

	// coefficient for the adjustment of torus size
	double TORUSSIZECOEF = 1.2;

	// coefficients for size calculation
	double sizecoef = 0.3;

	po::options_description options1("opt1");
	options1.add_options()
			("--ewidth", po::value<int>(&EWIDTH), "width of petals")
			("--eheight", po::value<int>(&EHEIGHT), "height of petals")
			("--minwidth", po::value<int>(&MINWIDTH), "minimum value of width of petals")
			("--minheight", po::value<int>(&MINHEIGHT), "minimum value of height of petals")
			("--shiftx", po::value<int>(&SHIFTX), "shift X from the origin(=upper-left end)")
			("--shifty", po::value<int>(&SHIFTY), "shift Y from the origin(=upper-left end)")
			("--petalposcoef", po::value<double>(&PETALPOSCOEF), "coeffcient for the adjustment of petal position")
			("--torussizecoef", po::value<double>(&TORUSSIZECOEF), "coefficient for the adjustment of torus size")
			;
	po::variables_map values;
	po::store(parse_command_line(argc, argv, options1), values);

	// read the result file
	string csvfname;
	string tabfname;

	if (argc > 3) {
		csvfname = argv[2];
		tabfname = argv[3];
	}
	double significance = -1.0;
	Flower_readfile fread;
	string resfname = argv[1];
	try {
		significance = fread.readResult(resfname, csvfname, tabfname);
	} catch (string &e) {
		cerr << e << endl;
		return 1;
	} catch (...) {
		cerr << "Error: An unexpected error occurred." << endl;
		return 1;
	}
	double colorvalMin = -log(1.0 / significance);

	// for each combination
	Flower_svg svg;
	for (int j = 0; j < (int)fread.combiName.size(); j++) {
		// reset the end position values
		double xmin = 1.0e+30;
		double xmax = -1.0e+30;
		double ymin = 1.0e+30;
		double ymax = -1.0e+30;

		// determine the scale of drawing
		vector<string> nameset = fread.combiName[j];
		for (int k = 0; k < (int)nameset.size(); k++) {
			// search for the motif and specify its ID
			int nameid = -1;
			for (int i = 0; i < (int)fread.motifName.size(); i++) {
				if (fread.motifName[i].find(nameset[k]) >= 0) {
					nameid = j;
					break;
				}
			}
			// calculate the end position values
			sizeval = log(fread.motifNgenes[nameid]) * sizecoef;
			rad = (double)k * svg.pi * 2 / (double)fread.combiName[j].size() - svg.pi * 0.5;
			sizex = EWIDTH * sizeval;
			sizey = EHEIGHT * sizeval;
			shiftx = (cos(rad) * (sizex - EHEIGHT * 2)) + SHIFTX;
			shifty = (sin(rad) * (sizex - EHEIGHT * 2)) + SHIFTY;
			if (xmax < shiftx + sizex) {
				xmax = shiftx + sizex;
			}
			if (xmin > shiftx - sizex) {
				xmin = shiftx - sizex;
			}
			if (ymax < shifty + sizey) {
				ymax = shifty + sizey;
			}
			if (ymin > shifty - sizey) {
				ymin = shifty - sizey;
			}
		}

		// calculate the scaling factor
		double scalex = (xmax - xmin) / (SHIFTX * 3);
		double scaley = (ymax - ymin) / (SHIFTY * 3);

		// open the SVG file
		nameset = fread.combiName[j];
		string svgfilename = resfname + "-flower";
		svgfilename += to_string(int(fread.combiRank[j]));
		svgfilename += ".svg";

		cout << svgfilename << endl;

		ofstream svgfile(svgfilename);
		svgfile << "<!DOCTYPE html>\n" << endl;
		svgfile << "  <svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" << endl;

		// for each motif in the combination
		for (int k = 0; k < (int)nameset.size(); k++) {
			// search for the motif and specify its ID
			int nameid = -1;
			for (int i = 0; i < (int)fread.motifName.size(); i++) {
				if (fread.motifName[i] == nameset[k]) {
					nameset[k];
					fread.motifName[i];
					nameid = i;
					break;
				}
			}
			// calculate the size and color of the ellipsoid
			double colorval;
			sizeval = log(fread.motifNgenes[nameid]) * sizecoef;
			if (fread.motifApvalue[nameid] > significance) {
				colorval = -log(fread.motifApvalue[nameid] / significance);
				if (colorval < colorvalMin) {
					colorval = colorvalMin;
				}
			}
			if (fread.motifApvalue[nameid] <= significance) {
				colorval = fread.motifApvalue[nameid] / significance;
			}
			// initialize variables before calculating the shape of the ellipsoid
			rad = (double)k * svg.pi * 2 / (double)fread.combiName[j].size() - svg.pi * 0.5;
			sizex = EWIDTH * sizeval * scalex;
			sizey = EHEIGHT * sizeval * scaley;
			if (sizex < MINWIDTH) {
				sizex = MINWIDTH;
			}
			if (sizey < MINHEIGHT) {
				sizey = MINHEIGHT;
			}
			shiftx = (cos(rad) * (sizex - EHEIGHT * PETALPOSCOEF)) + SHIFTX;
			shifty = (sin(rad) * (sizex - EHEIGHT * PETALPOSCOEF)) + SHIFTY;
			double annox = (cos(rad) * (sizex * 2 - EHEIGHT * PETALPOSCOEF)) + SHIFTX;
			double annoy = (sin(rad) * (sizex * 2 - EHEIGHT * PETALPOSCOEF)) + SHIFTY;

			// draw a motif
			svg.drawMotif(sizex, sizey, shiftx, shifty, rad, colorval, &svgfile);
			string name = fread.motifName[nameid];
			double pvalue = fread.motifApvalue[nameid];
			svg.annotateMotif(name, pvalue, annox, annoy, &svgfile);

			// draw the combination
			colorval = fread.combiApvalue[j] / significance;
			int size = EHEIGHT * scalex * TORUSSIZECOEF;
			if (size < MINHEIGHT) {
				size = MINHEIGHT;
			}
			double apvalue = fread.combiApvalue[j];
			svg.drawMotif(size, size, SHIFTX, SHIFTY, 0.0, colorval, &svgfile);
			svg.annotateMotif("", apvalue, (SHIFTX - 30), (SHIFTY - 5), &svgfile);
		}
		//# close the SVG file
		svgfile << "  </svg>\n" << endl;
	}
	return 0;
}

