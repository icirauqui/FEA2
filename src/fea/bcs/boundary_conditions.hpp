#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>

class BoundaryConditions {

public:

  BoundaryConditions(int num_dof, bool bypass = false);



  void AddNodalByNodeIds(std::vector<unsigned int> &node_ids, std::vector<double> &values);

  virtual void AddNodalByCoords(std::vector<double> coords, std::vector<bool> dof, std::vector<double> values) = 0;

  void AddNodalX(float x, float f);

  void AddNodalY(float y, float f);

  void AddNodalZ(float z, float f);

  void AddNodalXY(float x, float y, float fx, float fy);

  void AddNodalXYZ(float x, float y, float z, float fx, float fy, float fz);



  virtual void Encastre(std::vector<double> coords, std::vector<bool> dof) = 0;

  void EncastreX(float x);

  void EncastreY(float y);

  void EncastreZ(float z);

  void EncastreXY(float x, float y);

  void EncastreXYZ(float x, float y, float z);


  void Report();



  // Accessors
  bool NodeIds(unsigned int idx) { return _node_ids[idx]; }
  std::vector<bool>& NodeIds() { return _node_ids; }
  int NumDof() { return _num_dof; }
  std::vector<double>& Values(unsigned int idx) { return _values[idx]; }
  std::vector<std::vector<double>>& Values() { return _values; }


protected:

  int _num_nodes;
  int _num_dof;
  int _num_nodes_constrained;
  int _num_dof_free;
  double _tolerance = 1e-6;
  double _bypass = false;

  std::vector<bool> _node_ids;
  std::vector<std::vector<double>> _values;
};



#endif