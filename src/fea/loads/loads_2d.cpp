#include "loads_2d.hpp"


Loads2d::Loads2d(int num_dof, std::vector<Eigen::Vector2d>* nodes) 
  : Loads(num_dof) {

    _nodes = nodes;

    _node_ids = std::vector<bool>(nodes->size(), false);
    _values = std::vector<std::vector<double>>(nodes->size(), std::vector<double>(_num_dof, 0.0));
  }


void Loads2d::AddNodal(std::vector<unsigned int> &node_ids, std::vector<double> &values) {
  for (auto node : node_ids) {
    _node_ids[node] = true;
    _values[node] = values;
  }

  std::cout << "   Added loads in " << node_ids.size() << " nodes" << std::endl;
}



void Loads2d::AddNodal(std::vector<double> coords, std::vector<double> values) {
  std::vector<unsigned int> nodes;

  for (unsigned int node = 0; node < _nodes->size(); node++) {
    for (unsigned int dof = 0; dof < _num_dof; dof++) {
      if (coords[dof] < 0.0)
        continue;

      if (abs((*_nodes)[node](dof) - coords[dof]) < _tolerance) {
        nodes.push_back(node);
        break;
      }
    }
  }

  AddNodal(nodes, values);
}


