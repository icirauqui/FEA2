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
  AbaqusC3D8_3 model;

  std::cout << "\nBuild FEA" << std::endl;
  FEA fea("C3D8", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in z = 0" << std::endl;
  bc.AddNodal(Eigen::Vector3d(-1.0, -1.0, 0.0), {0, 0, 0}, {0.0, 0.0, 0.0});
  std::cout << " - Force of magnitude 1 in direction of z on nodes in z = 3" << std::endl;
  bc.AddNodal(Eigen::Vector3d(-1.0, -1.0, 3.0), {1, 1, 1}, {0.0, 0.0, 1.0});
  //bc.Report();


  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/abaqus_c3d8_3/K_a.csv");
  
  return 0;

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ExportK("../data/abaqus_c3d8_3/K_b.csv");



  //std::cout << "\nComputeDisplacements" << std::endl;
  //fea.ComputeDisplacements();



  // FEA Solver
  std::cout << "\nSolve" << std::endl;
  fea.Solve();



  std::cout << "\nSolver Done" << std::endl;

  std::cout << "\nCompute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  std::cout << " - u.size(): " << u.size() << std::endl;
  std::cout << " - model._nodes.size(): " << model._nodes.size() << std::endl;
  double scale = 10.0;
  model.ApplyDisplacements(u, scale);

  std::cout << "\nExport All" << std::endl;
  fea.ExportAll("../data/abaqus_c3d8_2/abaqus_c3d8_2.csv");




  std::cout << "\nLen nodes/deformed: " << model._nodes.size() << " / " << model._nodes_deformed.size() << std::endl;
  for (unsigned int i=0; i<5; i++) {
    std::cout << "u[" << i << "] = " 
              << model._nodes[i].transpose() << " ->\t"
              << scale * u[i].transpose() << " ->\t"
              << model._nodes_deformed[i].transpose() << std::endl;
  }


  //std::cout << "\nReport" << std::endl;
  //fea.ReportNodes("../data/abaqus_c3d8_2/nodes.csv");
  //fea.ExportK("../data/abaqus_c3d8_2/K.csv");


  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.0, 0.0, 1.0));
  viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(1.0, 0.0, 0.0));
  viewer.AddLoads(bc.NodeIds(), bc.Values());
  viewer.AddEdges(model._elements);
  viewer.Render();

  
  return 0;
}
