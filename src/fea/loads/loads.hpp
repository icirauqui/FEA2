#ifndef LOADS_HPP
#define LOADS_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

class Loads {

public:

  Loads(int num_dof) {
    _num_dof = num_dof;
  }

  // Accessors
  bool NodeIds(unsigned int idx) { return _node_ids[idx]; }
  std::vector<bool>& NodeIds() { return _node_ids; }
  std::vector<double>& Values(unsigned int idx) { return _values[idx]; }
  std::vector<std::vector<double>>& Values() { return _values; }

  void Report() {
    std::cout << std::endl;
    std::cout << "Loads3d Report" << std::endl;
    std::cout << "  _num_dof: " << _num_dof << std::endl;
    std::cout << "  _node_ids.size(): " << _node_ids.size() << std::endl;
    std::cout << "  _values.size(): " << _values.size() << std::endl;

    unsigned int num_bcs = 0;
    std::cout << "  Boundary Conditions:" << std::endl;
    for (unsigned int i = 0; i < _node_ids.size(); i++) {
      if (!_node_ids[i]) 
        continue;
      num_bcs++;
      
      std::cout << "    " << i << ":";
      std::cout << " values(";
      for (unsigned int j = 0; j < _num_dof; j++) {
        std::cout << " " << _values[i][j];
      }
      std::cout << " )" << std::endl;
    }

    std::cout << "  Number of nodes with BCs: " << num_bcs << std::endl;
  }


protected:
  int _num_dof;
  double _tolerance = 1e-6;

  std::vector<bool> _node_ids;
  std::vector<std::vector<double>> _values;
};



#endif