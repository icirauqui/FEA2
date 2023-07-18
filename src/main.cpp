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

float noise_multiplier = 0.3;
Eigen::Vector3d model_offset(1.0, 2.0, 3.0); 

Eigen::Vector3d ApproximatePose(std::vector<Eigen::Vector3d> pts) {
  Eigen::Vector3d centroid(0.0, 0.0, 0.0);
  for (int i = 0; i < pts.size(); i++) {
    centroid += pts[i];
  }
  centroid /= pts.size();

  double dist = sqrt(pow(centroid(0), 2) + pow(centroid(1), 2) + pow(centroid(2), 2));
  centroid(2) += dist/2;

  return centroid;
}

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
  FEM fem2(element);

  // Add points
  for (int i = 0; i < vpts.size(); i++) {
    Eigen::Vector3d pt(vpts[i][0], vpts[i][1], vpts[i][2]);
    Eigen::Vector3d pt2(vpts[i][0], vpts[i][1], vpts[i][2]);
    for (unsigned int j=0; j<3; j++) {
      pt2(j) += noise_multiplier * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    }
    pt2 += model_offset;
    fem.AddPoint(pt);
    fem2.AddPoint(pt2);
  }

  fem.Compute(true);
  fem2.InitCloud();
  fem.ComputeExtrusion();

  
  FEA fea(0, element, E, nu, depth, fg, false);

  Eigen::Vector3d pose1 = ApproximatePose(fem.GetEigenNodes());

  fem.ViewMesh(true, fem2.GetCloud(), pose1);
}


int main(int argc, char** argv) {
  test_fe();
  
  return 0;
}
