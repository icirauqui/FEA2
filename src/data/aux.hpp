#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;


// Load data from csv file with custom delimiter.
vector<vector<float> > get_from_file_vvf (string inputPath, string delimiter);

vector<vector<int> > get_from_file_vvn (string inputPath, string delimiter);


void put_to_file_vvf (string outputPath, string delimiter, vector<vector<float> > vvoutput, bool append);


vector<vector<float> > vector_resize_cols(vector<vector<float> > v1, int n);

#endif