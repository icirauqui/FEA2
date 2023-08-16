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
  quaternion.setFromTwoVectors(Eigen::Vector3d::UnitX(), direction);
  Eigen::Vector4d qvec = quaternion.coeffs();

  return std::make_pair(qvec, tvec);
}

std::vector<Eigen::Vector3d> SimulateSteps(Eigen::Vector3d offset, int steps) {
  std::vector<Eigen::Vector3d> poses;
  Eigen::Vector3d delta = offset / steps;
  for (int i = 0; i < steps; i++) {
    poses.push_back(offset - (delta * i));
  }
  return poses;
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
  fem2.SetExtrusion(fem.GetExtrusionDelta(), fem.GetElementHeight());

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose1 = ApproximatePose(fem.GetEigenNodes());
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose2 = ApproximatePose(fem2.GetEigenNodes());
  std::cout << "Pose1 = " << pose1.second.transpose() << std::endl;
  std::cout << "Pose2 = " << pose2.second.transpose() << std::endl;

  fem.ViewMesh(true, fem2.GetCloud(), fem2.GetExtrusion(), pose1, pose2);

  FEA fea(0, element, E, nu, depth, fg, false);
  std::vector<std::vector<float>> nodes = fem.GetNodes();
  std::vector<std::vector<int>> elements = fem.GetElements();
  fea.MatAssembly(nodes, elements);


  // Create POS object
  POS pos(fem2.GetNodes(), pose2);






  // Impose a rotation of x degrees around each axis
  double imposed_angle = pos.DegToRad(30.0);
  //Eigen::Vector3d imposed_angle_v(imposed_angle, imposed_angle, imposed_angle);
  //Eigen::Vector3d imposed_angle_v(imposed_angle, 0.0, 0.0);
  Eigen::Vector3d imposed_angle_v(0.0, imposed_angle, 0.0);
  //Eigen::Vector3d imposed_angle_v(0.0, 0.0, imposed_angle);
  //Eigen::Vector3d imposed_angle_v(imposed_angle, imposed_angle, imposed_angle);
  Eigen::Vector4d imposed_angle_q = pos.EulerToQuaternion(imposed_angle_v);

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

  std::cout << "Requested rotation = " << imposed_angle_q.transpose() << std::endl;
  std::cout << "Applied rotation = " << applied_rotation.transpose() << std::endl;

  int num_steps = 5;
  std::vector<Eigen::Vector3d> steps = SimulateSteps(pose2.second-pose1.second, num_steps);
  Eigen::Vector3d delta = (pose2.second-pose1.second) / num_steps;
  double step_size = 1.0 / static_cast<double>(num_steps);

  std::vector<Eigen::Vector4d> steps_rot;
  
  
//  for (auto step: steps) {
//    //Eigen::Vector4d identity_quaternion = Eigen::Quaterniond::Identity().coeffs();
//    Eigen::Vector4d identity_quaternion(1.0, 0.0, 0.0, 0.0);
//    std::pair<Eigen::Vector4d, Eigen::Vector3d> latest_pose = pos.GetPose();
//    std::cout << "  Translation = " << delta.transpose() << std::endl;
//
//    pos.Transform(identity_quaternion, identity_quaternion, -delta, 1.0);
//
//    std::pair<Eigen::Vector4d, Eigen::Vector3d> new_pose = pos.GetPose();
//    std::cout << "  Pose = " << latest_pose.second.transpose() << " -> " << new_pose.second.transpose() << std::endl;
//
//
//    std::vector<Eigen::Vector3d> new_nodes = pos.GetPoints();
//    std::vector<Eigen::Vector3d> new_nodes_front, new_nodes_back;
//    for (unsigned int i=0; i<new_nodes.size(); i++) {
//      if (i < new_nodes.size()/2) {
//        new_nodes_front.push_back(new_nodes[i]);
//      } else {
//        new_nodes_back.push_back(new_nodes[i]);
//      }
//    }
//    fem.ViewMesh(true, new_nodes_front, new_nodes_back, pose1, new_pose);
//    std::cout << std::endl;
//  }
}


int main(int argc, char** argv) {
  //test_fea();
  test_fe();
  
  return 0;
}
