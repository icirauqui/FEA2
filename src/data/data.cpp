#include "data.hpp"

Dataset::Dataset(std::string data_path, std::string element_type) {
  path_ = data_path;
  element_type_ = element_type;

  vpts_ = get_from_file_vvf(path_ + "/input_" + element_type_ + "/input_points.csv",",");
  vElts_ = get_from_file_vvn(path_ + "/input_" + element_type_ + "/input_elements.csv",",");
  vvF_ = get_from_file_vvf(path_ + "/input_" + element_type_ + "/input_forces.csv",",");
  vDir_ = get_from_file_vvn(path_ + "/input_" + element_type_ + "/input_Dirichlet.csv",",");

  // Node numbers in Abaqus start at 1, change to 0 for cpp
  for (unsigned int i=0; i<vElts_.size(); i++){
      for (unsigned int j=0; j<vElts_[i].size(); j++){
          vElts_[i][j] -= 1;
      }
  }
}

std::vector<std::vector<float>> Dataset::points() {
  return vpts_;
}

std::vector<std::vector<int>> Dataset::elements() {
  return vElts_;
}

std::vector<std::vector<float>> Dataset::forces() {
  return vvF_;
}

std::vector<std::vector<int>> Dataset::dirichlet() {
  return vDir_;
}