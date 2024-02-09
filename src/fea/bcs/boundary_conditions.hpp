#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

class BoundaryConditions {

public:

  BoundaryConditions(int num_dof) {
    _num_dof = num_dof;
    _num_dof_constrained = 0;
  }

  //virtual void AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values) = 0;
  //virtual void AddNodal(Eigen::Vector3d coords, std::vector<unsigned int> dof, std::vector<double> values) = 0;

  // Accessors
  bool NodeIds(unsigned int idx) { return _node_ids[idx]; }
  std::vector<bool>& NodeIds() { return _node_ids; }
  int NumDof() { return _num_dof; }
  std::vector<unsigned int>* Dof(unsigned int idx) { return _dof[idx]; }
  std::vector<std::vector<unsigned int>*> Dof() { return _dof; }
  std::vector<double>& Values(unsigned int idx) { return _values[idx]; }
  std::vector<std::vector<double>>& Values() { return _values; }

  void Report() {
    std::cout << std::endl;
    std::cout << "BoundaryConditions3d Report" << std::endl;
    std::cout << "  _num_dof: " << _num_dof << std::endl;
    std::cout << "  _num_dof_constrained: " << _num_dof_constrained << std::endl;
    std::cout << "  _num_dof_free: " << _num_dof_free << std::endl;
    std::cout << "  _node_ids.size(): " << _node_ids.size() << std::endl;
    std::cout << "  _dof.size(): " << _dof.size() << std::endl;
    std::cout << "  _values.size(): " << _values.size() << std::endl;

    unsigned int num_bcs = 0;
    std::cout << "  Boundary Conditions:" << std::endl;
    for (unsigned int i = 0; i < _node_ids.size(); i++) {
      if (!_node_ids[i]) 
        continue;
      num_bcs++;
      
      std::cout << "    " << i << ":";
      std::cout << "  dof(";
      for (unsigned int j = 0; j < _num_dof; j++) {
        std::cout << " " << (*_dof[i])[j];
      }
      std::cout << " )  values(";
      for (unsigned int j = 0; j < _num_dof; j++) {
        std::cout << " " << _values[i][j];
      }
      std::cout << " )" << std::endl;
    }

    std::cout << "  Number of nodes with BCs: " << num_bcs << std::endl;
  }


protected:

  int _num_nodes;
  int _num_dof;
  int _num_dof_constrained;
  int _num_dof_free;
  double _tolerance = 1e-6;

  //std::vector<Eigen::Vector3d>* _nodes;

  std::vector<bool> _node_ids;
  std::vector<std::vector<unsigned int>*> _dof;
  std::vector<std::vector<double>> _values;
};



#endif