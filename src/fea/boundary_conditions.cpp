#include "boundary_conditions.hpp"


BoundaryConditions::BoundaryConditions(int num_dof, std::vector<Eigen::Vector3d>* nodes) 
  : _num_dof(num_dof), _num_dof_constrained(0), _nodes(nodes) {
    _num_dof_free = nodes->size();

    _node_ids = std::vector<unsigned int>(_num_dof_free, -1);
    _dof = std::vector<std::vector<unsigned int>*>(_num_dof_free, nullptr);
    _values = std::vector<std::vector<double>>(_num_dof_free, std::vector<double>(_num_dof, 0.0));
  }

void BoundaryConditions::AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values) {
  std::vector<unsigned int>* dof_vec = new std::vector<unsigned int>(_num_dof, 0);
  for (unsigned int i = 0; i < _num_dof; i++) {
    if (dof[i] == 1) {
      (*dof_vec)[i] = 1;
      _num_dof_constrained++;
      _num_dof_free--;
    }
  }
  for (auto node : node_ids) {
    _node_ids.push_back(node);
    _dof.push_back(dof_vec);
    _values.push_back(values);
  }
}

void BoundaryConditions::AddNodal(Eigen::Vector3d coords, std::vector<double> &values) {
  std::vector<unsigned int>* dof_vec = new std::vector<unsigned int>(_num_dof, 0);
  for (unsigned int i = 0; i < _num_dof; i++) {
    if (coords(i) >= 0.0) {
      (*dof_vec)[i] = 1;
      _num_dof_constrained++;
      _num_dof_free--;
    }
  }

  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int dof = 0; dof < _num_dof; dof++) {
      if (coords(dof) < 0.0)
        continue;

      if ((*_nodes)[node](dof) - coords(dof) < _tolerance) {

        // check if node is in _node_ids
        bool found = false;
        for (auto n : _node_ids) {
          if (n == node) {
            found = true;
            break;
          }
        }


        _node_ids.push_back(node);
        _dof[node] = dof_vec; //.push_back(dof_vec);
        _values[node] = values; //.push_back(values);
        break;
      }
    }
  }
}