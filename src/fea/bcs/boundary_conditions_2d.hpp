#ifndef BOUNDARY_CONDITIONS_2D_HPP
#define BOUNDARY_CONDITIONS_2D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "boundary_conditions.hpp"

class BoundaryConditions2d : public BoundaryConditions {

public:

  BoundaryConditions2d(int num_dof, std::vector<Eigen::Vector2d>* nodes);

  void AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values);
  void AddNodal(std::vector<double> coords, std::vector<unsigned int> dof, std::vector<double> values);

  void Encastre(std::vector<double> coords);



private:

  std::vector<Eigen::Vector2d>* _nodes;
  
};



#endif