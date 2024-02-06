#ifndef LOADS_2D_HPP
#define LOADS_2D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "loads.hpp"

class Loads2d : public Loads {

public:

  Loads2d(int num_dof, std::vector<Eigen::Vector2d>* nodes);

  void AddNodal(std::vector<unsigned int> &node_ids, std::vector<double> &values);
  void AddNodal(std::vector<double> coords, std::vector<double> values);




private:

  std::vector<Eigen::Vector2d>* _nodes;
  
};



#endif