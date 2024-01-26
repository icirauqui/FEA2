#include "boundary_conditions.hpp"


BoundaryConditions::BoundaryConditions(int num_dof, std::vector<Eigen::Vector3d>* nodes) 
  : _num_dof(num_dof), _num_dof_constrained(0), _nodes(nodes) {
    _num_dof_free = nodes->size();

    _node_ids = std::vector<bool>(_num_dof_free, false);
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
    _node_ids[node] = true; //.push_back(node);
    _dof[node] = dof_vec; //.push_back(dof_vec);
    _values[node] = values; //.push_back(values);
  }

  std::cout << "Added " << node_ids.size() << " new boundary conditions" << std::endl;
}

void BoundaryConditions::AddNodal(Eigen::Vector3d coords, std::vector<unsigned int> dof, std::vector<double> values) {
  // Transform dof into a pointer to a vector
  std::vector<unsigned int>* pdof = new std::vector<unsigned int>(_num_dof, 0);
  for (unsigned int i = 0; i < _num_dof; i++) {
    (*pdof)[i] = dof[i];
  }

  int num_new_bcs = 0;
  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int dof = 0; dof < _num_dof; dof++) {
      if (coords(dof) < 0.0)
        continue;

      if (abs((*_nodes)[node](dof) - coords(dof)) < _tolerance) {
        //std::cout << "node - coords = " << (*_nodes)[node](dof) << " - " << coords(dof) << " = " << (*_nodes)[node](dof) - coords(dof) << std::endl;
        //std::cout << " Press enter to continue " << std::endl;
        //std::cin.get();
        _node_ids[node] = true; //.push_back(node);
        _dof[node] = pdof; //.push_back(dof_vec);
        _values[node] = values; //.push_back(values);
        num_new_bcs++;
        break;
      }
    }
  }

  std::cout << "   Added " << num_new_bcs << " new boundary conditions" << std::endl;
}




void BoundaryConditions::Report() {
  std::cout << std::endl;
  std::cout << "BoundaryConditions Report" << std::endl;
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