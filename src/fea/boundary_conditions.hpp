#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

class BoundaryConditions {

public:

  BoundaryConditions(int num_dof, std::vector<Eigen::Vector3d>* nodes);

  void AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values);
  void AddNodal(Eigen::Vector3d coords, std::vector<unsigned int> dof, std::vector<double> values);

  // Accessors
  bool NodeIds(unsigned int idx) { return _node_ids[idx]; }
  std::vector<bool>& NodeIds() { return _node_ids; }
  std::vector<unsigned int>* Dof(unsigned int idx) { return _dof[idx]; }
  std::vector<std::vector<unsigned int>*> Dof() { return _dof; }
  std::vector<double>& Values(unsigned int idx) { return _values[idx]; }
  std::vector<std::vector<double>>& Values() { return _values; }

  void Report();


private:

  int _num_nodes;
  int _num_dof;
  int _num_dof_constrained;
  int _num_dof_free;
  double _tolerance = 1e-6;

  std::vector<Eigen::Vector3d>* _nodes;

  std::vector<bool> _node_ids;
  std::vector<std::vector<unsigned int>*> _dof;
  std::vector<std::vector<double>> _values;
};



#endif