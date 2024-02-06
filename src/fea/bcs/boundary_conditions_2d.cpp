#include "boundary_conditions_2d.hpp"


BoundaryConditions2d::BoundaryConditions2d(int num_dof, std::vector<Eigen::Vector2d>* nodes) 
  : BoundaryConditions(num_dof) {

    _num_dof_free = nodes->size();
    _nodes = nodes;

    _node_ids = std::vector<bool>(_num_dof_free, false);
    _dof = std::vector<std::vector<unsigned int>*>(_num_dof_free, nullptr);
    _values = std::vector<std::vector<double>>(_num_dof_free, std::vector<double>(_num_dof, 0.0));
  }

void BoundaryConditions2d::AddNodal(std::vector<unsigned int> &node_ids, std::vector<unsigned int> &dof, std::vector<double> &values) {
  std::vector<unsigned int>* dof_vec = new std::vector<unsigned int>(_num_dof, 0);
  for (unsigned int i = 0; i < _num_dof; i++) {
    if (dof[i] == 1) {
      (*dof_vec)[i] = 1;
      _num_dof_constrained++;
      _num_dof_free--;
    }
  }

  for (auto node : node_ids) {
    _node_ids[node] = true;
    _dof[node] = dof_vec;
    _values[node] = values;
  }

  std::cout << "   Added " << node_ids.size() << " new boundary conditions" << std::endl;
}



void BoundaryConditions2d::AddNodal(std::vector<double> coords, std::vector<unsigned int> dof, std::vector<double> values) {
  // Transform dof into a pointer to a vector
  std::vector<unsigned int>* pdof = new std::vector<unsigned int>(_num_dof, 0);
  for (unsigned int i = 0; i < _num_dof; i++) {
    (*pdof)[i] = dof[i];
  }

  int num_new_bcs = 0;
  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int dof = 0; dof < _num_dof; dof++) {
      if (coords[dof] < 0.0)
        continue;

      if (abs((*_nodes)[node](dof) - coords[dof]) < _tolerance) {
        _node_ids[node] = true;
        _dof[node] = pdof;
        _values[node] = values;
        num_new_bcs++;
        break;
      }
    }
  }

  std::cout << "   Added " << num_new_bcs << " new boundary conditions" << std::endl;
}





void BoundaryConditions2d::Encastre(std::vector<double> coords) {
  // Transform dof into a pointer to a vector
  std::vector<unsigned int>* pdof = new std::vector<unsigned int>(_num_dof, 0);
  std::vector<double> values = {0.0, 0.0};

  int num_new_bcs = 0;
  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int dof = 0; dof < _num_dof; dof++) {
      if (coords[dof] < 0.0)
        continue;

      if (abs((*_nodes)[node](dof) - coords[dof]) < _tolerance) {
        _node_ids[node] = true;
        _dof[node] = pdof;
        _values[node] = values;
        num_new_bcs++;
        break;
      }
    }
  }

  std::cout << "   Added encastre in " << num_new_bcs << " nodes" << std::endl;
}

