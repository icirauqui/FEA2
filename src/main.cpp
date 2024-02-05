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


void test_c3d8() {
  AbaqusC3D8_2 model;

  std::cout << "\nBuild FEA" << std::endl;
  FEA fea("C3D8", E, nu, true);

  std::cout << "\nBuild BoundaryConditions" << std::endl;
  BoundaryConditions3d bc(fea.NumDof(), &model._nodes);
  std::cout << " - Encastre in z = 0" << std::endl;
  bc.AddNodal(Eigen::Vector3d(-1.0, -1.0, 0.0), {0, 0, 0}, {0.0, 0.0, 0.0});
  std::cout << " - Force of magnitude 1 in direction of z on nodes in z = " << model.LoadLocation() << std::endl;
  bc.AddNodal(Eigen::Vector3d(-1.0, -1.0, model.LoadLocation()), {1, 1, 1}, {0.0, 0.0, 1.0});
  //bc.Report();


  std::cout << "\nMatAssembly" << std::endl;
  fea.MatAssembly(model._nodes, model._elements);
  fea.ExportK("../data/" + model.Name() + "/K_a.csv");

  //return 0;

  std::cout << "\nApplyBoundaryConditions" << std::endl;
  fea.ApplyBoundaryConditions(bc);
  fea.ExportK("../data/" + model.Name() + "/K_b.csv");

  //std::cout << "\nComputeDisplacements" << std::endl;
  //fea.ComputeDisplacements();

  // FEA Solver
  std::cout << "\nSolve" << std::endl;
  fea.Solve("BiCGSTAB");

  std::cout << "\nCompute deformed positions" << std::endl;
  Eigen::MatrixXd U = fea.U();
  std::vector<Eigen::Vector3d> u;
  for (unsigned int n=0; n<model._nodes.size(); n++) {
    u.push_back(U.block(n*3, 0, 3, 1));
  }
  std::cout << " - u.size(): " << u.size() << std::endl;
  std::cout << " - model._nodes.size(): " << model._nodes.size() << std::endl;
  double scale = 1000.0;
  model.ApplyDisplacements(u, scale);

  fea.ExportAll("../data/" + model.Name() + "/KFU.csv");
  fea.ReportNodes("../data/" + model.Name() + "/nodes.csv");

  std::cout << "\nVisualization" << std::endl;
  PCLViewer viewer;
  viewer.AddNodes(model._nodes, "original", Eigen::Vector3d(0.2, 0.2, 1.0));
  viewer.AddNodes(model._nodes_deformed, "deformed", Eigen::Vector3d(0.0, 0.7, 0.0));
  viewer.AddLoads(bc.NodeIds(), bc.Values());
  viewer.AddEdges(model._elements);
  viewer.Render();
}


//void test_c2d4() {
//  AbaqusC2D4_1 model("../py/Dogbone_Tension.input");
//
//  std::cout << "\nBuild FEA" << std::endl;
//  FEA fea("C2D4", E, nu, true);
//
//  std::cout << "\nBuild BoundaryConditions" << std::endl;
//  BoundaryConditions3d bc(fea.NumDof(), &model._nodes);
//}


int main(int argc, char** argv) {

  test_c3d8();
  
  return 0;
}
