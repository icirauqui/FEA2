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

#include <elts/element.hpp>
#include <elts/element3d.hpp>
#include <elts/c3d6.hpp>
#include <elts/c3d8.hpp>

#include <boundary_conditions.hpp>



class FEA {

public:


  FEA(std::string element_type,
      float young_modulus, float poisson_coefficient,
      bool debug_mode);
    
  void MatAssembly(std::vector<Eigen::Vector3d> &vpts, 
                   std::vector<std::vector<unsigned int>> &velts);

  void ApplyBoundaryConditions(BoundaryConditions &bc);




  // Legacy fea

  void EncastreBackLayer();
  void ImposeDirichletEncastre(std::vector<int> &dir);
  void ImposeDirichletEncastre(std::vector<std::vector<int>> &dir);

  void ComputeForces();
  void SetForces(std::vector<std::vector<float>> &vF);

  void ComputeDisplacements();

  void ComputeStrainEnergy();
  double ComputeStrainEnergy(std::vector<Eigen::Vector3d> &u0,
                             std::vector<Eigen::Vector3d> &u1);

  // Reports

  void ReportNodes(std::string filename);
  void ExportK(std::string filename);
  void PrintK();
  void PrintEigenvaluesK();

  // Accessors

  Eigen::MatrixXd K() { return K_; }
  Eigen::MatrixXd F() { return F_; }
  Eigen::MatrixXd U() { return U_; }
  float StrainEnergy() { return sE_; }
  int NumNodes() { return base_size_; }
  int NumDof() { return element_->getDofPerNode(); }





private:

  Element* element_;

  Eigen::MatrixXd K_;
  Eigen::MatrixXd K1_;
  Eigen::MatrixXd Kei_;

  Eigen::MatrixXd F_;
  Eigen::MatrixXd U_;

  float k_large_ = 1e8;

  float sE_ = 0.0;

  float E_ = 1.0;
  float nu_ = 0.499;

  bool debug_mode_ = false;
};

#endif