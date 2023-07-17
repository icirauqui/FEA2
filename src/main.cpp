#include <iostream>
#include <vector>
#include <math.h>

#include "fea/fem.hpp"
#include "fea/fea.hpp"
#include "fea/pos.hpp"
#include "dataset/dataset.hpp"

#include <chrono>


// Parameters
int image_id = 0;
std::string element = "C3D6";
float E = 3500.0;
float nu = 0.495;
float depth = 1.0;
float fg = 0.577350269;
float Klarge = 100000000.0;


void test_fea() {
  // Load data
  Dataset ds("../data", element);
  std::vector<std::vector<float>> vpts = ds.points();
  std::vector<std::vector<int>>   vElts = ds.elements();
  std::vector<std::vector<float>> vvF = ds.forces();
  std::vector<std::vector<int>>   vDir = ds.dirichlet();

  // Test new class
  FEA fea(0, element, E, nu, depth, fg, false);
  
  auto start1 = std::chrono::high_resolution_clock::now();
  fea.MatAssembly(vpts, vElts);
  fea.SetForces(vvF);
  fea.ImposeDirichletEncastre(vDir, Klarge);
  fea.ComputeDisplacements();
  fea.ComputeStrainEnergy();
  auto stop1 = std::chrono::high_resolution_clock::now();

  std::cout << " - Strain energy = " << fea.StrainEnergy() << std::endl;
  std::cout << "   Time fea1 = " << std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1).count() << std::endl;
}

void test_fe() {

  // Load data
  Dataset ds("../data", element);
  std::vector<std::vector<float>> vpts = ds.points();

  FEM fem(element);

  // Add points
  for (int i = 0; i < vpts.size(); i++) {
    Eigen::Vector3d pt(vpts[i][0], vpts[i][1], vpts[i][2]);
    fem.AddPoint(pt);
  }

  fem.Compute(true);

  fem.ViewMesh();

  std::vector<std::vector<float>> nodes = fem.GetNodes();
  std::vector<std::vector<int>> elements = fem.GetElements();
  
  FEA fea(0, element, E, nu, depth, fg, false);





}


int main(int argc, char** argv) {
  test_fe();
  
  return 0;
}
