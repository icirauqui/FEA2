#ifndef FEA_HPP
#define FEA_HPP

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <chrono>
#include <bits/stdc++.h>
#include <eigen3/Eigen/Dense>

#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/point_types.h>

#include <pcl/common/common.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/surface/mls.h>
#include <pcl/surface/impl/mls.hpp>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/visualization/cloud_viewer.h>

#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>

#include <pcl/surface/gp3.h>

#include <pcl/visualization/vtk.h>

#include <pcl/console/parse.h>
#include <pcl/io/vtk_lib_io.h>

#include <elements/c3d6.hpp>
#include <elements/c3d8.hpp>


class FEA {

public:

  FEA(int frame_id, 
      std::string element_type,
      float young_modulus, float poisson_coefficient, float element_depth, 
      float gauss_point, 
      bool debug_mode);
    
  // Finite Element Analysis

  void MatAssembly(std::vector<std::vector<float> > &vpts, 
                   std::vector<std::vector<int> > &velts);

  void Print_K();
  void Eigenvalues_K();

  void ComputeForces();

  void SetForces(std::vector<std::vector<float>> &vF);


  void EncastreBackLayer(float k_large = 1e8);
  void ImposeDirichletEncastre(std::vector<int> &dir, float k_large = 1e8);
  void ImposeDirichletEncastre(std::vector<std::vector<int>> &dir, float k_large = 1e8);

  void ComputeDisplacements();

  void ComputeStrainEnergy();
  double ComputeStrainEnergy(std::vector<Eigen::Vector3d> &u0,
                             std::vector<Eigen::Vector3d> &u1);

  // Accessors
  Eigen::MatrixXd K();
  Eigen::MatrixXd F();
  Eigen::MatrixXd U();
  float StrainEnergy();


private:


  // Finite Element Analysis

  void InitC3D6();

  void InitC3D8();

  void InitGaussPoints(float fg);
  
  //void ComputeKei(std::vector<std::vector<float>> &vfPts);
  
  void dNdgs(float xi, float eta, float zeta, int dim);








  int frame_id_ = 0;
  std::string element_;

  Eigen::MatrixXd D_ = Eigen::MatrixXd::Zero(6, 6);
  Eigen::Matrix<float, 8, 3> gs_;

  Eigen::MatrixXd K_;
  Eigen::MatrixXd K1_;
  Eigen::MatrixXd Kei_;

  Eigen::MatrixXd F_;
  Eigen::MatrixXd U_;

  float sE_ = 0.0;

  float E_ = 1.0;
  float nu_ = 0.499;
  float h_ = 1.0;

  float lambda_ = 0.0;
  float G_ = 0.0;
  int base_size_ = 0;

  Eigen::MatrixXd dndgs_;
  Eigen::Matrix3d J_;
  Eigen::Matrix3d J1_;
  Eigen::MatrixXd B_;

  bool debug_mode_ = false;
};

#endif