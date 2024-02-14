#include "boundary_conditions.hpp"



BoundaryConditions::BoundaryConditions(int num_dof, bool bypass) {
  _num_dof = num_dof;
  _num_nodes_constrained = 0;
  _bypass = bypass;
}

void BoundaryConditions::AddNodalByNodeIds(std::vector<unsigned int> &node_ids, std::vector<double> &values) {
  for (auto node : node_ids) {
    if (_node_ids[node]) {
      std::cout << "   Warning: Node " << node << " already has a load: [";
      for (unsigned int i = 0; i < _num_dof; i++) {
        std::cout << " " << _values[node][i];
      }
      std::cout << " ]" << std::endl;
      if (!_bypass)
        continue;
    }
    _node_ids[node] = true;
    _values[node] = values;
    _num_nodes_constrained++;
  }
  std::cout << "   Added loads in " << node_ids.size() << " nodes" << std::endl;
}

void BoundaryConditions::AddNodalX(float x, float f) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  std::vector<double> values = std::vector<double>(_num_dof, 0.0);
  values[0] = f;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalY(float y, float f) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[1] = true;
  std::vector<double> values = std::vector<double>(_num_dof, 0.0);
  values[1] = f;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalZ(float z, float f) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[2] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[2] = true;
  std::vector<double> values = std::vector<double>(_num_dof, 0.0);
  values[2] = f;
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalXY(float x, float y, float fx, float fy) {
  std::vector<double> coords = {x, y};
  std::vector<bool> dof = {true, true};
  std::vector<double> values = {fx, fy};
  AddNodalByCoords(coords, dof, values);
}

void BoundaryConditions::AddNodalXYZ(float x, float y, float z, float fx, float fy, float fz) {
  std::vector<double> coords = {x, y, z};
  std::vector<bool> dof = {true, true, true};
  std::vector<double> values = {fx, fy, fz};
  AddNodalByCoords(coords, dof, values);
}


void BoundaryConditions::EncastreX(float x) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  Encastre(coords, dof);
}

void BoundaryConditions::EncastreY(float y) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[1] = true;
  Encastre(coords, dof);
}

void BoundaryConditions::EncastreZ(float z) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[1] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[1] = true;
  Encastre(coords, dof);
}

void BoundaryConditions::EncastreXY(float x, float y) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  coords[1] = y;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  dof[1] = true;
  Encastre(coords, dof);
}

void BoundaryConditions::EncastreXYZ(float x, float y, float z) {
  std::vector<double> coords = std::vector<double>(_num_dof, 0.0);
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
  std::vector<bool> dof = std::vector<bool>(_num_dof, false);
  dof[0] = true;
  dof[1] = true;
  dof[2] = true;
  Encastre(coords, dof);
}



void BoundaryConditions::Report() {
  std::cout << std::endl;
  std::cout << "BoundaryConditions3d Report" << std::endl;
  std::cout << "  _num_dof: " << _num_dof << std::endl;
  std::cout << "  _node_ids.size(): " << _node_ids.size() << std::endl;
  std::cout << "  _num_nodes_constrained: " << _num_nodes_constrained << std::endl;
  std::cout << "  _num_nodes_free: " << _node_ids.size() - _num_nodes_constrained << std::endl;

  unsigned int num_bcs = 0;
  std::cout << "  Boundary Conditions:" << std::endl;
  for (unsigned int i = 0; i < _node_ids.size(); i++) {
    if (!_node_ids[i]) 
      continue;
    num_bcs++;
    
    std::cout << "    " << i << ":  values(";
    for (unsigned int j = 0; j < _num_dof; j++) {
      std::cout << " " << _values[i][j];
    }
    std::cout << " )" << std::endl;
  }
  std::cout << "  Number of nodes with BCs: " << num_bcs << std::endl;
}



