#include "boundary_conditions.hpp"


BoundaryConditions::BoundaryConditions(int num_dof, std::vector<Eigen::Vector3d>* nodes) 
  : _num_dof(num_dof), _num_dof_constrained(0), _nodes(nodes) {
    _num_dof_free = nodes->size();
  }

void BoundaryConditions::AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values) {
  for (auto node : node_ids) {
    _node_ids.push_back(node);
    _dof.push_back(dof);
    _values.push_back(values);
  }
}

void BoundaryConditions::AddNodal(Eigen::Vector3d coords, std::vector<double> &values) {
  std::vector<unsigned int> dof_vec(_num_dof, 0);
  for (unsigned int i = 0; i < _num_dof; i++) {
    if (coords(i) >= 0.0) {
      dof_vec[i] = 1;
    }
  }

  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int dof = 0; dof < _num_dof; dof++) {
      if (coords(dof) < 0.0)
        continue;

      if ((*_nodes)[node](dof) - coords(dof) < _tolerance) {
        _node_ids.push_back(node);
        _dof.push_back(dof_vec);
        _values.push_back(values);
        break;
      }
    }
  }
}