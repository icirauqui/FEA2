#ifndef BOUNDARY_CONDITIONS_3D_HPP
#define BOUNDARY_CONDITIONS_3D_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "boundary_conditions.hpp"

class BoundaryConditions3d : public BoundaryConditions {

public:

  BoundaryConditions3d(int num_dof, std::vector<Eigen::Vector3d>* nodes);

  void AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values);
  void AddNodal(Eigen::Vector3d coords, std::vector<unsigned int> dof, std::vector<double> values);



private:

  std::vector<Eigen::Vector3d>* _nodes;
  
};



#endif