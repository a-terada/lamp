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
#include <ctime>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

#include "LampCore.h"

using namespace boost::program_options;

/** @file */
/**
 * Run multiple testing correction.
 * This script need transaction file, expression-file and significance-level.
 * transaction-file: The file includes associations between TFs and genes.
 *     Each line indicates a gene.
 *     If gene is targeted by the TF, then value is 1, otherwise 0.
 * expression-file: Each line indicates a gene. The column1 is gene name.
 *     If gene has the feature, the column2 is 1. The other case 0.
 * significance-level: The statistical significance threshold.
 * @author Terada 26, June, 2011
 * @author Miura 31, 12, 2015 Convert Python program to C++ code.
 */
int main(int argc, char** argv) {
	int term_code = 1;
	// command line options
	std::string pvalue_procedure, max_comb_str, log_filename, alternative_str;
	
	options_description opt("options");
	opt.add_options()
		("help,h",																		"print help message.")
		("pvalue,p",	value<std::string>(&pvalue_procedure),							"Choose the p-value calculation procedure from 'fiehser' (Fisher's exact test), 'chi' (Chi-square test) or 'u_test' (Mann-Whitney's U-test)")
		("max_comb",	value<std::string>(&max_comb_str)->default_value("all"),		"Set the maximum size of combination to be tested.")
		("log_file,e",	value<std::string>(&log_filename),								"The file name to output log.")
		("alternative",	value<std::string>(&alternative_str)->default_value("greater"),	"Indicate which alternative hypothesis is used. Select \"greater\", \"less\" or \"two.sided\"\n, and the default is \"greater\".")
	;
	try {
		variables_map values;
		auto const parsing_result = parse_command_line(argc, argv, opt);
		store(parsing_result, values);
		notify(values);
		
		std::vector<std::string> args = collect_unrecognized(parsing_result.options, include_positional);
		
		if( values.count("help") ) {
			std::cerr << "usage: lamp [options] transaction_file value_file significance_probability" << std::endl;
			std::cerr << opt << std::endl;
			exit(0);
		}
		// check arguments
		if (args.size() != 3) {
			throw std::string("Error: input [target-file], [expression-file] and [significance-level].");
		}
		
		int max_comb = -1;
		std::transform(max_comb_str.begin(), max_comb_str.end(), max_comb_str.begin(), ::tolower);
		if (!max_comb_str.compare("all") == 0) {
			if (std::all_of(max_comb_str.cbegin(), max_comb_str.cend(), ::isdigit)) {
				max_comb = std::stoi(max_comb_str);
			} else {
				throw std::string("Error: max_comb must be an integer value or all.");
			}
		}
		
		// check the file exist.
		std::string transaction_file = args[0];
		std::ifstream tfile(transaction_file);
		if (!tfile.good()) {
			throw std::string("IOError: No such file: \'" + transaction_file + "\'");
		}
		tfile.close();
		std::string flag_file = args[1];
		std::ifstream vfile(flag_file);
		if (!vfile.good()) {
			throw std::string("IOError: No such file: \'" + flag_file + "\'");
		}
		vfile.close();
		
		double threshold = -1;
		try {
			threshold = std::stod(args[2]);
		} catch (std::exception &ex) {}
		if ((threshold < 0) || (threshold > 1)) {
			throw std::string("Error: significance probability must be a float value from 0.0 to 1.0.");
		}
		
		// check the value of alternative hypothesis
		int alternative = 0;
		if (alternative_str.compare("greater") == 0) {
			alternative = 1;
		}
		else if (alternative_str.compare("less") == 0) {
			alternative = -1;
		}
		else if (alternative_str.compare("two.sided") == 0) {
			alternative = 0;
		}
		else {
			throw std::string("Error: \"alternative\" should be one of {\"greater\", \"less\", \"two.sided\"}");
		}
	
		// change log file
		std::string log_file;
		if (0 < log_filename.size()) {
			log_file = log_filename;
		} else {
			time_t now;
			std::time(&now);
			struct tm* d = localtime(&now);
			int year = d->tm_year + 1900;
			int month = d->tm_mon + 1;
			int day = d->tm_mday;
			int hour = d->tm_hour;
			int minute = d->tm_min;
			int second = d->tm_sec;
			log_file = "lamp_log_";// = "lamp_log_" + d.strftime("%Y%m%d") + "_" + d.strftime("%H%M%S") + ".txt"
			log_file += std::to_string(year) + std::to_string(month) + std::to_string(day) + "_";
			log_file += std::to_string(hour) + std::to_string(minute) + std::to_string(second) + ".txt";
		}
		LampCore lampCore;
		lampCore.run(transaction_file, flag_file, threshold, pvalue_procedure,
						  max_comb, log_file, alternative);
		lampCore.print(transaction_file, flag_file, threshold, pvalue_procedure,
						  alternative);
		
		term_code = 0;
	} catch (boost::program_options::error &e) {
		std::cerr << "option parsing error: " << e.what() << std::endl << opt << std::endl;
	} catch (std::string &msg) {
		std::cerr << msg << std::endl;
	} catch (...) {
		std::cerr << "Error: An unexpected error occurred." << std::endl;
	}
	
	return term_code;
}
