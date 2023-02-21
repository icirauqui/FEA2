#include <iostream>
#include <vector>
#include <math.h>

#include "data/aux.hpp"
#include "fea/fea.hpp"
#include "data/data.hpp"
#include "fea/fea2.hpp"

#include <chrono>


// Parameters
int image_id = 0;
std::string element = "C3D8";
float E = 3500.0;
float nu = 0.495;
float depth = 1.0;
float fg = 0.577350269;
float Klarge = 100000000.0;


int main(int argc, char** argv) {
  // Load data
  Dataset ds("../data", element);
  std::vector<std::vector<float>> vpts = ds.points();
  std::vector<std::vector<int>>   vElts = ds.elements();
  std::vector<std::vector<float>> vvF = ds.forces();
  std::vector<std::vector<int>>   vDir = ds.dirichlet();

  // Test new class
  FEA fea(0, element, E, nu, depth, fg, false);
  FEA2 fea2(1,E,nu,1,fg,false);

  auto start1 = std::chrono::high_resolution_clock::now();
  fea.MatAssembly(vpts, vElts);
  fea.SetForces(vvF);
  fea.ImposeDirichletEncastre(vDir, Klarge);
  fea.ComputeDisplacements();
  fea.ComputeStrainEnergy();
  auto stop1 = std::chrono::high_resolution_clock::now();

  std::cout << " - Strain energy = " << fea.StrainEnergy() << std::endl;
  std::cout << "   Time fea1 = " << std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1).count() << std::endl;


  // Test old class
  auto start2 = std::chrono::high_resolution_clock::now();
  fea2.K = fea2.MatrixAssemblyC3D8(vpts,vElts);
  fea2.vF = vector_resize_cols(vvF,1);
  fea2.ImposeDirichletEncastre(vDir,Klarge);
  fea2.K1 = fea2.InvertMatrixEigen(fea2.K);
  fea2.vU = fea2.MultiplyMatricesEigen(fea2.K1,fea2.vF);
  std::vector<std::vector<float> > vvU2 = vector_resize_cols(fea2.vU,3);
  std::vector<std::vector<float> > vUt2;
  std::vector<float> vUti2;
  for (unsigned int i=0; i<fea2.vU.size(); i++){
      vUti2.push_back(fea2.vU[i][0]);
  }
  vUt2.push_back(vUti2);
  std::vector<std::vector<float> > sE2 = fea2.MultiplyMatricesEigen(vUt2,fea2.vF);
  auto stop2 = std::chrono::high_resolution_clock::now();

  std::cout << " - Strain energy = " << sE2[0][0] << std::endl;
  std::cout << "   Time fea2 = " << std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2).count() << std::endl;



  return 0;
}
