#ifndef DATA_HPP
#define DATA_HPP

#include <iostream>
#include <vector>
#include <string>

#include "aux.hpp"


class Dataset {

public:
  Dataset(std::string data_path, std::string element_type);

  std::vector<std::vector<float>> points();
  std::vector<std::vector<int>> elements();
  std::vector<std::vector<float>> forces();
  std::vector<std::vector<int>> dirichlet();



private:
  std::vector<std::vector<float>> vpts_;
  std::vector<std::vector<int>> vElts_;
  std::vector<std::vector<float>> vvF_;
  std::vector<std::vector<int>> vDir_;

  std::string path_;
  std::string element_type_;

};


#endif