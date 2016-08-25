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
#include <sstream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <set>
#include <getopt.h>

using namespace std;

/** @file */
/**
 * output usage
 */
void print_usage(string program_name){
	cerr << "Usage: " << program_name << " filename" << endl;
	cerr << "Options:" << endl;
	cerr << "  -h                  show this help message and exit." << endl;
	cerr << "  -o OUTPUT_FILENAME  The file name to output the eliminated result." << endl;
}

/** Struct to manage LAMP data */
struct detection{
	int rank; /**< The rank ordered by the p-value. */
	double p_val_raw; /**< The p-value calculated using p-value-procedure. */
	double p_val_adjusted; /**< The adjusted p-value. */
	string combination; /**< The significant combination. */
	int arity; /**< The number of elements in the combination. */
	int arity_of_target_rows; /**< The number of target genes of the combination. */
	int arity_positives_in_the_targets; /**< The value used to compute the p-value. */
	
	/**
	 * Comparator
	 */
	bool operator<(const detection &right) const{
		if(p_val_raw < right.p_val_raw){
			return true;
		}else if (p_val_raw > right.p_val_raw) {
			return false;
		}else if(arity < right.arity){
			return true;
		}else{
			return false;
		}
	}
};
/** Struct to manage result */
struct result{
	vector<detection> detections_list; /**< Store the combinations */
	vector<string> meta_line_list;/**< Store the meta lines */
	string time_line; /**< Store the line containing running time */
	vector<string> minp_line_list; /**< Store the minimum P-value lines */
};

/**
 * read the result of LAMP.
 * @param in_filename the result filename of LAMP
 * @param res result data
 */
void readResult(string *in_filename, result *res){
	bool flag_broken = true; // a flag to decide whether the file is broken or not.
	bool flag_comb = false;
	bool flag_minp_dist = false; // a flag for checking that the line is minimum P-value distribution.
//	int detection_size = 0; // number of detection sets.
	
	ifstream ifs(in_filename->c_str(), ios::in);
	if (ifs.fail()) {
		cerr << *in_filename << " cannot be opened." << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	while (getline(ifs, line)) {
		if ( line.find("Time (sec.): ", 0) == 0 ) {
			res->time_line = line;
			flag_broken = false;
			flag_comb = false;
			continue;
		}
		if ( line == "--- minimum P-values ---") {
			res->minp_line_list.push_back(line);
			flag_minp_dist = true;
			continue;
		}
		if (flag_minp_dist){
			res->minp_line_list.push_back(line);
			continue;
		}
		
		if (!flag_comb and !flag_minp_dist){
			res->meta_line_list.push_back(line);
			if ( line.find("Rank", 0) == 0) {
				flag_comb = true;
			}
			continue;
		}
		
		vector<string> detections;
		stringstream ss(line);
		string buffer;
		while(getline(ss, buffer, '\t')){
			detections.push_back(buffer);
		}
		
		detection det;
		det.rank = atoi(detections.at(0).c_str());
		det.p_val_raw = atof(detections.at(1).c_str());
		det.p_val_adjusted = atof(detections.at(2).c_str());
		det.combination = detections.at(3);
		det.arity = atoi(detections.at(4).c_str());
		det.arity_of_target_rows = atoi(detections.at(5).c_str());
		det.arity_positives_in_the_targets = atoi(detections.at(6).c_str());
		res->detections_list.push_back(det);
	}
	ifs.close();
	
	if (flag_broken) {
		cerr << in_filename << " is broken." << endl;
		exit(EXIT_FAILURE);
	}
}

/**
 * Return true, if vec_b is contained in vec_a
 * @param vec_a
 * @param vec_b
 * @return bool
 */
bool isSubset(vector<string> *vec_a, vector<string> *vec_b){
	bool res = true;
	for ( int i = 0; i< vec_b->size(); ++i){
		int hitnum = count( vec_a->begin(), vec_a->end(), vec_b->at(i) );
		if (hitnum == 0){
			res = false;
		}
	}
	return res;
}

/**
 * Eliminate the redundant combinations.
 * @param detections_list a list of combinations made by readResult function.
 * @param merged_list
 */
void mergeResult(vector<detection> *detections_list, vector<detection> *merged_list){
	vector< vector<string> > upper_list;
	
	for ( int i = 0; i< detections_list->size(); ++i){
		detection det = detections_list->at(i);
		set<string> detect_set;
		stringstream ss(det.combination);
		string buffer;
		while(getline(ss, buffer, ',')){
			detect_set.insert(buffer);
		}
		vector<string> detect_vec(detect_set.begin(), detect_set.end());
		
		bool flag = true;
		
		for ( int j = 0; j< upper_list.size(); ++j){
			vector<string> upper_detect_vec = upper_list.at(j);
			if( isSubset(&detect_vec, &upper_detect_vec) || isSubset(&upper_detect_vec, &detect_vec) ){
				flag = false;
				break;
			}
		}
		if (flag) {
			merged_list->push_back(det);
		}
		upper_list.push_back(detect_vec);
	}
}

/**
 * output result to file
 * @param out_filename output filename
 * @param merged_list
 * @param res resulr data
 */
void output(string *out_filename, vector<detection> *merged_list, result *res){
	FILE* outstream;
	// If set output filename, standard output redirect to that file.
	if (out_filename->length() > 0) {
		if((outstream = freopen(out_filename->c_str(), "w", stdout)) == NULL){
			cerr << "Error: Cannot output to " << *out_filename << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	cout << "# Non-redundant combinations" << endl;
	for ( int i = 0; i< res->meta_line_list.size(); ++i){
		string line = res->meta_line_list.at(i);
		cout << line;
		if ( line.find("# # of significant combinations: ", 0) == 0 ) {
			cout << " -> " << merged_list->size();
		}
		cout << endl;
	}
	int rank = 0;
	for ( int i = 0; i< merged_list->size(); ++i){
		rank = ++rank;
		cout << rank << "\t";
		detection det = merged_list->at(i);
		cout << det.p_val_raw << "\t";
		cout << det.p_val_adjusted << "\t";
		cout << det.combination << "\t";
		cout << det.arity << "\t";
		cout << det.arity_of_target_rows << "\t";
		cout << det.arity_positives_in_the_targets << endl;
	}
	cout << res->time_line << endl;
	for ( int i = 0; i< res->minp_line_list.size(); ++i){
		cout << res->minp_line_list.at(i) << endl;
	}
	if (out_filename->length() > 0) {
		fclose(outstream);
	}
}

/**
 * Run function
 * @param in_filename the filename of LAMP
 * @param out_filename the filename to output the eliminated result.
 */
void run(string *in_filename, string *out_filename){
	// Read result
	result res;
	readResult(in_filename, &res);
	
	// Sort result
	sort(res.detections_list.begin(), res.detections_list.end());	
	// Merge result
	vector<detection> merged_list;
	mergeResult(&(res.detections_list), &merged_list);
	
	// Output to file
	output(out_filename, &merged_list, &res);
}

/** 
 * Remove the redundant combinations in the result of LAMP.
 */
int main(int argc, char** argv) {
	string in_filename;
	string out_filename = "";
	
	// parse options
	int c;
	while(( c = getopt(argc, argv, "o:h")) != -1){
		switch (c){
			case 'o':
				out_filename = optarg;
				break;
			case 'h':
				print_usage(argv[0]);
				exit(EXIT_FAILURE);
			default:
				print_usage(argv[0]);
				exit(EXIT_FAILURE);
		}
	}
	
	// check arguments
	if (argc - optind != 1){
		print_usage(argv[0]);
		exit(EXIT_FAILURE);
	}else{
		in_filename = argv[optind++];
	}
	
	run(&in_filename, &out_filename);
	return 0;
}

