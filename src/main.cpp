#include <iostream>
#include <vector>
#include <math.h>

#include "dataset/test_models.cpp"
#include "fea/fea.hpp"
#include "vis/vis.hpp"

// Parameters
double E = 10000.0;
double nu = 0.495;
float Klarge = 100000000.0;



int main(int argc, char** argv) {
  AbaqusC3D8_2 model;

  std::cout << "\nBuild FEA" << std::endl;
  FEA fea("C3D8", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in z = 0" << std::endl;
  bc.AddNodal(Eigen::Vector3d(-1.0, -1.0, 0.0), {0, 0, 0}, {0.0, 0.0, 0.0});
  std::cout << " - Force of magnitude 1 in direction of z on nodes in z = 25" << std::endl;
  bc.AddNodal(Eigen::Vector3d(-1.0, -1.0, 5.0), {1, 1, 1}, {0.0, 0.0, 1.0});
  //bc.Report();


  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);

  std::cout << "\nComputeDisplacements" << std::endl;
  fea.ComputeDisplacements();

  std::cout << "\nCompute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  double scale = 100.0;
  model.ApplyDisplacements(u, scale);

  //std::cout << "\nLen nodes: " << model._nodes.size() << "; Len _nodes_deformed: " << model._nodes_deformed.size() << std::endl;
  ////for (unsigned int i=0; i<u.size(); i++) {
  //for (unsigned int i=0; i<5; i++) {
  //  std::cout << "u[" << i << "] = " 
  //            << model._nodes[i].transpose() << " ->\t"
  //            << scale * u[i].transpose() << " ->\t"
  //            << model._nodes_deformed[i].transpose() << std::endl;
  //}


  std::cout << "\nReport" << std::endl;
  fea.ReportNodes("../data/abaqus_c3d8_2/nodes.csv");
  fea.ExportK("../data/abaqus_c3d8_2/K.csv");


  //std::cout << "\nVisualization" << std::endl;
  //PCLViewer viewer(true);
  //viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.0, 0.0, 1.0));
  //viewer.AddEdges(model._elements, "original", Eigen::Vector3d(0.0, 0.0, 1.0));
  //viewer.AddLoads(bc.NodeIds(), bc.Values(), 1.0);
  //viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(1.0, 0.0, 0.0));
  //viewer.AddEdges(model._elements, "deformed", Eigen::Vector3d(1.0, 0.0, 0.0));
  //viewer.Render();

  
  return 0;
}
