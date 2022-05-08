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
vector<vector<float> > get_from_file_vvf (string inputPath, string delimiter) {
    ifstream file(inputPath);
    vector<vector<string> > datalist;
    string line = "";

    while (getline(file,line)) {
        vector<string> vec;
        boost::algorithm::split(vec,line,boost::is_any_of(delimiter));
        datalist.push_back(vec);
    }
    file.close();

    vector<vector<float> > vPoints;
    for (unsigned int i=0; i<datalist.size(); i++) {
        vector<float> vPoint;
        for (unsigned int j=0; j<datalist[i].size(); j++) {
            string sPoint = datalist[i][j];
            float fPoint = strtof((sPoint).c_str(),0);
            vPoint.push_back(fPoint);
        }
        vPoints.push_back(vPoint);
    }

    return vPoints;
}


vector<vector<int> > get_from_file_vvn (string inputPath, string delimiter) {
    ifstream file(inputPath);
    vector<vector<string> > datalist;
    string line = "";

    while (getline(file,line)) {
        vector<string> vec;
        boost::algorithm::split(vec,line,boost::is_any_of(delimiter));
        datalist.push_back(vec);
    }
    file.close();

    vector<vector<int> > vPoints;
    for (unsigned int i=0; i<datalist.size(); i++) {
        vector<int> vPoint;
        for (unsigned int j=0; j<datalist[i].size(); j++) {
            string sPoint = datalist[i][j];
            int nPoint = stoi((sPoint).c_str(),0);
            vPoint.push_back(nPoint);
        }
        vPoints.push_back(vPoint);
    }

    return vPoints;
}


void put_to_file_vvf (string outputPath, string delimiter, vector<vector<float> > vvoutput, bool append){
    
    ofstream os;
    if (append)
        os.open(outputPath.c_str(), ios::out | ios::app );
    else
        os.open(outputPath.c_str(), ios::out);
    
    for (unsigned int i=0; i<vvoutput.size(); i++) {
        for (unsigned int j=0; j<vvoutput[i].size(); j++) {
            os << vvoutput[i][j] << " ";
        }
        os << endl;
    }
    os << endl;
    
    os.close();
}


vector<vector<float> > vector_resize_cols(vector<vector<float> > v1, int n){
    vector<vector<float> > v2;
    vector<float> v2i;
    for (unsigned int i=0; i<v1.size(); i++){
        for (unsigned int j=0; j<v1[i].size(); j++){
            v2i.push_back(v1[i][j]);
            if (v2i.size() == n){
                v2.push_back(v2i);
                v2i.clear();
            }
        }
    }
    return v2;
}

#endif