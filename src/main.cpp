#include <iostream>
#include <vector>
#include <math.h>

#include "fea/fem.hpp"
#include "fea/fea.hpp"
#include "fea/pos.hpp"
#include "dataset/dataset.hpp"
#include "nlo/levenberg_marquardt.hpp"

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

  // 1. Load data
  Dataset ds("../data", element);
  std::vector<std::vector<float>> vpts = ds.points();



  // 2. Create FEM objects and add points
  FEM fem(element);
  FEM fem2(element);

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



  // 3. Triangulate andcompute poses.
  fem.Compute(true);
  fem2.InitCloud();
  fem.ComputeExtrusion();
  fem2.SetExtrusion(fem.GetExtrusionDelta(), fem.GetElementHeight());

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1 = ApproximatePose(fem.GetEigenNodes());
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = ApproximatePose(fem2.GetEigenNodes());

  fem.ViewMesh(true, fem2.GetCloud(), fem2.GetExtrusion(), pose1, pose2);



  // 4. Transform Pose 2 for simulation, impose a rotation of x degrees around each axis
  POS pos(fem2.GetNodes(), pose2);
  double ang = 30*M_PI/180;
  Eigen::Vector3d axis(1,1,1);
  Eigen::Vector4d imposed_angle_q = pos.QuaternionFromAngleAxis(axis, ang);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> latest_pose = pos.GetPose();
  pos.Transform(imposed_angle_q, Eigen::Vector3d(0, 0, 0), 1.0);

  std::vector<Eigen::Vector3d> new_nodes = pos.GetPoints();
  std::vector<Eigen::Vector3d> new_nodes_front, new_nodes_back;
  for (unsigned int i=0; i<new_nodes.size(); i++) {
    if (i < new_nodes.size()/2) {
      new_nodes_front.push_back(new_nodes[i]);
    } else {
      new_nodes_back.push_back(new_nodes[i]);
    }
  }
  std::pair<Eigen::Vector4d, Eigen::Vector3d> new_pose = pos.GetPose();
  fem.ViewMesh(true, new_nodes_front, new_nodes_back, pose1, new_pose);

  Eigen::Vector4d applied_rotation = pos.ComputeQuaternionRotation(latest_pose.first, new_pose.first);



  // 5. Initialize FEA object, use conectivity to compute stiffness matrix
  FEA fea(0, element, E, nu, depth, fg, false);
  std::vector<std::vector<float>> nodes = fem.GetNodes();
  std::vector<std::vector<int>> elements = fem.GetElements();
  fea.MatAssembly(nodes, elements);



  // 6. Simulate 5 steps of translation and rotation, model 2 will move towards model 1.
  int num_steps = 5;
  Eigen::Vector3d step = (pose2.second-pose1.second) / num_steps;
  Eigen::Vector4d step_rot = pos.QuaternionFromAngleAxis(axis, -ang/num_steps);
  
  for (unsigned int s=0; s<num_steps; s++) {
    std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = pos.GetPose();
    std::vector<Eigen::Vector3d> nodes_k1 = pos.GetPoints();

    pos.Transform(step_rot, -step, 1.0);

    std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k = pos.GetPose();

    std::vector<Eigen::Vector3d> nodes_k = pos.GetPoints();
    std::vector<Eigen::Vector3d> nodes_k_front, nodes_k_back;
    for (unsigned int i=0; i<nodes_k.size(); i++) {
      if (i < nodes_k.size()/2) {
        nodes_k_front.push_back(nodes_k[i]);
      } else {
        nodes_k_back.push_back(nodes_k[i]);
      }
    }

    double sE = fea.ComputeStrainEnergy(nodes_k1, nodes_k);

    std::cout << std::endl;
    std::cout << "Step " << s << std::endl;
    std::cout << "  Translation = " << pose_k1.first.transpose() << " -> " << pose_k.first.transpose() << std::endl;
    std::cout << "  Orientation = " << pose_k1.second.transpose() << " -> " << pose_k.second.transpose() << std::endl;
    std::cout << "  Strain energy = " << sE << std::endl;

    fem.ViewMesh(true, nodes_k_front, nodes_k_back, pose1, pose_k);
  }
}


int main(int argc, char** argv) {
  //test_fea();
  test_fe();
  
  return 0;
}
