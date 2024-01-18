#include <iostream>
#include <vector>
#include <math.h>

//#include "fea/fem.hpp"
//#include "fea/fea.hpp"
//#include "fea/pos.hpp"
//#include "dataset/dataset.hpp"
//#include "nlo/levenberg_marquardt.hpp"


#include "test_models.cpp"

#include "fea/fea.hpp"
#include "fea/elts/c3d6.hpp"
#include "fea/elts/c3d8.hpp"

// include for std::ofstream
#include <fstream>

#include <chrono>

// Parameters
double E = 10000.0;
double nu = 0.495;
float Klarge = 100000000.0;

Eigen::Vector3d model_offset(0.2, 0.2, 0.2); 

std::pair<Eigen::Vector4d, Eigen::Vector3d> ApproximatePose(std::vector<Eigen::Vector3d> pts) {
  Eigen::Vector3d centroid(0.0, 0.0, 0.0);
  for (int i = 0; i < pts.size(); i++) {
    centroid += pts[i];
  }
  centroid /= pts.size();

  double dist = sqrt(pow(centroid(0), 2) + pow(centroid(1), 2) + pow(centroid(2), 2));
  Eigen::Vector3d tvec = centroid + Eigen::Vector3d(0, 0, dist/2);

  // Compute rotation, from t to centroid, as a quaternion
  Eigen::Vector3d direction = centroid - tvec;
  direction.normalize();
  Eigen::Quaterniond quaternion;
  quaternion.setFromTwoVectors(-Eigen::Vector3d::UnitX(), direction);
  Eigen::Vector4d qvec = quaternion.coeffs();

  return std::make_pair(qvec, tvec);
}



int main(int argc, char** argv) {
  AbaqusC3D8_1 model;
  C3D8 fea_obj = C3D8(E, nu);
  Eigen::MatrixXd D = fea_obj.getD();
  std::cout << "Elasticity Matrix: \n" << D << std::endl;
  std::cout << "E: " << fea_obj.getE() << std::endl;
  std::cout << "nu: " << fea_obj.getNu() << std::endl;

  //Eigen::MatrixXd K = fea_obj.matAssembly(model._nodes, model._elements);
  //std::cout << "Stiffness Matrix: \n" << K << std::endl;

  FEA fea("C3D8", E, nu, false);

  // Save to csv
  //std::ofstream file;
  //file.open("stiffness_matrix.csv");
  //file << K << std::endl;
  //file.close();


  //Eigen::MatrixXd K = c3d8::matAssembly(model._nodes, model._elements, E, nu);
  //std::cout << "Stiffness Matrix: \n" << K << std::endl;
  
  return 0;
}
