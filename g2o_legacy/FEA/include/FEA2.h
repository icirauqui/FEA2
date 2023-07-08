/**
 * This file is an addon to ORB-SLAM2.
 *
 * Its goal is to create a mesh from the 3D points within the computed map.
 * Also, it will assembly a rigidity matrix for such mesh, and compute the
 * strain deformation energy.
 *
 * Copyright (C) 2017 Íñigo Cirauqui Viloria <icirauquiviloria at gmail dot com>
 * (University of Zaragoza)
 *
 * This addon, as ORB-SLAM2, is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This addon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * This file is part of ORB-SLAM2.
 *
 * Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University
 * of Zaragoza) For more information see <https://github.com/raulmur/ORB_SLAM2>
 *
 * ORB-SLAM2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ORB-SLAM2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FEA2_H
#define FEA2_H

#include <math.h>

#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#include "../../types/types_seven_dof_expmap.h"

// ORB_SLAM2 Libraries
#include "../../../../../src/core/point3d.cpp"

// OpenCV Libraries
#include <opencv4/opencv2/core/core.hpp>
#include <opencv4/opencv2/highgui/highgui.hpp>
#include <opencv4/opencv2/opencv.hpp>

// PCL Libraries
#include <pcl/common/common.h>
#include <pcl/common/common_headers.h>
#include <pcl/console/parse.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_types.h>
#include <pcl/surface/gp3.h>
#include <pcl/surface/mls.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <boost/thread/thread.hpp>
#include <pcl/surface/impl/mls.hpp>

// Other Libraries
#include <pcl/visualization/vtk.h>

#include <Eigen/StdVector>
#include <algorithm>

using namespace std;
using namespace pcl;

class FEA2 {
 public:  // FUNCTIONS
          // Constructor & Destructor
  FEA2(unsigned int input_nFrameId, unsigned int input_E, float input_nu,
       float input_h, float input_fg1, int input_nElType, bool bSetDebug);
  ~FEA2();

  // Build mesh and compute K
  bool Compute();
  void ComputeTestMesh();
  bool Compute(int nMode = 1);
  bool ComputePlanar();

  // Get coordinates and normals from Points
  void GetMapPointCoordinates();
  void ClearMapPointCoordinates();
  void AddMapPointCoordinates(std::vector<float> vPoint);

  // Load Points into cloud
  bool LoadMPsIntoCloud(int nMode = 1);

  // Moving Least Squares
  bool MLS(int nMode = 1);

  // Triangulation, greedy projection (meshing)
  bool ComputeMesh(int nMode = 1);

  // Triangulate over a planar mesh
  bool Triangulate();

  std::vector<std::vector<int>> GetMeshTriangles(int nMode = 1);

  void CalculateGP3Parameters(pcl::PointCloud<pcl::PointNormal>::Ptr ppc,
                              int *pnInMu, int *pnSearchRad, int *pnMaxNeig,
                              int *pnSurAng, int *pnMinAng, int *pnMaxAng);

  // Turns a triangle mesh into a quadrilateral mesh, 1 tri = 3 quads
  // The function adds a node in the middle of each side of the original
  // triangle, and an extra node in the ortocenter.
  void tri2quad();

  // Duplicates the previously generated mesh at a distance, builds a 3D
  // structure.
  void SetSecondLayer(int nMode = 1);

  vector<vector<float>> ComputeKeiC3D8(vector<vector<float>> vfPts);
  vector<vector<float>> ComputeKeiC3D6(vector<vector<float>> vfPts);

  bool MatrixAssemblyC3D8(int nMode = 1);
  bool MatrixAssemblyC3D6(int nMode = 1);

  void ImposeDirichletEncastre_K(int nMode, vector<vector<int>> vD,
                                 float Klarge);
  void ImposeDirichletEncastre_a(vector<vector<int>> vD, float Klarge);

  vector<vector<float>> InvertMatrixEigen(vector<vector<float>> m1);
  vector<vector<float>> MultiplyMatricesEigen(vector<vector<float>> m1,
                                              vector<vector<float>> m2);

  void Set_u0(vector<Point3D *> vpMPs, int nMode = 1);
  void Set_uf_self();
  void Set_uf(vector<vector<float>> vPoints);
  void SetDisplacement(std::vector<std::vector<float>> vDisplacement);
  void ComputeDisplacement();
  void ComputeNewDisplacement();

  void ComputeForces();

  float ComputeStrainEnergy();
  float GetStrainEnergy();
  float NormalizeStrainEnergy();
  float GetNormalizedStrainEnergy();

  void ComputeInitialStrainEnergy(vector<vector<float>> vPoints);

  void UpdateForces();

  vector<vector<float>> vector_resize_cols(vector<vector<float>> v1,
                                           unsigned int n);
  vector<vector<int>> vector_resize_cols_int(vector<vector<int>> v1,
                                             unsigned int n);
  vector<float> ComputeScaledDef(vector<vector<float>> vf, float scalerange);

  void set_camera_pose(int model, Eigen::Vector4d qvec, Eigen::Vector3d tvec);
  void add_point3d(int model, Point3D *pMP);

  void setbfea(bool bSet);
  void setbfea2(bool bSet);

  // Visualization
  void DisplayPolygons(std::string &img_path,
                       std::vector<std::vector<int>> &tri,
                       std::vector<std::vector<float>> &points_2d,
                       std::vector<int> color = {0, 150, 0});

  cv::Mat SetTransparentColor(cv::Mat &img,
                              std::vector<std::vector<cv::Point>> &roi,
                              double alpha1, double alpha2);

  cv::Scalar SetColor(float floatValue);

  void set_mesh_settings(float mls_search_radius, int mls_polynomial_order,
                         float mesh_mu, float mesh_search_radius,
                         int mesh_max_neighbours, int mesh_surf_angle,
                         int mesh_min_angle,
                         int mesh_max_angle);

public :  // VARIABLES
                vector<Point3D *> vpMPs_t;
  vector<Point3D *> vpMPs_ut;
  vector<cv::KeyPoint *> vpKPs_t;

  std::vector<Eigen::Vector4d> qvecs_;
  std::vector<Eigen::Vector3d> tvecs_;

  std::vector<std::vector<int>> vpMPs_pairs_;

  vector<vector<float>> vMPsXYZN_t;
  vector<vector<float>> vMPsXYZN_t2;
  vector<vector<float>> vMPsXYZN_ut;
  vector<vector<float>> vMPsXYZN_ut2;
  vector<vector<float>> vMPsXYZN_u;
  vector<vector<float>> vMPsXYZN_u2;

  // IndexTrackedMpInFEM - IndexTrackedMpInFrame
  vector<int> idxMpF;

  vector<g2o::VertexSBAPointXYZ *> vVertices;

  vector<vector<Point3D *>> vpMPs2Draw;
  vector<float> vpMPs2DrawWgt;
  vector<vector<Point3D *>> vpMPs2Drawu;
  vector<int> vpMPs2DrawuWgt;
  vector<vector<cv::KeyPoint *>> vpKPs2Draw;
  vector<vector<cv::KeyPoint *>> vpKPs2Drawu;

  vector<vector<int>> vvDir_t;
  vector<vector<int>> vvDir_u;

  bool bEInverse = false;

  // Set Points as active or not
  vector<bool> vbtMPsActive;
  vector<bool> vbuMPsActive;

  // Staring PointCloud
  pcl::PointCloud<pcl::PointXYZ> pc_t_0;
  pcl::PointCloud<pcl::PointXYZ>::Ptr pc_0;

  // Smoothed PointCloud (MLS)
  pcl::PointCloud<pcl::PointNormal> pc_t_1;
  pcl::PointCloud<pcl::PointNormal> pc_u_1;
  vector<int> vtindices;
  vector<int> vuindices;

  // Reconstructed mesh
  PolygonMesh mesh_t;
  PolygonMesh mesh_u;

  vector<vector<int>> triangles_t;
  vector<vector<int>> triangles_u;
  vector<vector<vector<int>>> edgespertriangle_t;
  vector<vector<vector<int>>> edgespertriangle_u;

  vector<vector<int>> vIdxNewVertices;
  vector<vector<int>> vIdxNewVertices_u;
  vector<vector<int>> vNewPointsBase;
  vector<vector<int>> vNewPointsBase_u;

  vector<vector<int>> quads_t;
  vector<vector<int>> quads_u;

  // Frame
  unsigned int nFrameId;

  // Young Modulus [Pa] & Poisson Coefficient
  unsigned int E;
  float nu;

  // Lamé parameters & Behaviour matrix
  float lambda = 0.0;
  float G = 0.0;
  vector<vector<float>> D;

  // Element depth
  float h;

  // Gauss Points
  float fg;
  vector<vector<float>> gs;

  // Elemental matrix
  vector<vector<float>> Ke;

  // Stiffness matrix
  unsigned int Ksize = 0;
  vector<vector<float>> K;
  float DetK = 0.0;

  unsigned int Kusize = 0;
  vector<vector<float>> Ku;
  float DetKu = 0.0;

  vector<vector<float>> K1;
  float DetK1 = 0.0;

  // Displacement
  vector<float> u0;   // Starting position of mappoints
  vector<float> u0u;  // Starting position of mappoints
  vector<float> uf;   // New position after being modified by the optimization
  vector<vector<float>> vva;
  vector<vector<float>> vva2;

  // Forces
  vector<vector<float>> vvf;
  vector<int> svvf;

  // Strain energy
  float sE, sE0;
  float nsE, nsE0;
  float StartingSE;
  float CurrentSE;
  float TempSE;
  float LastSE = 0.0;

  // Normalization vector
  float fNormFactor;

  bool bInFEA = false;
  bool bInFEA2 = false;

  bool it0 = false;

  // Interface
  float mls_search_radius_ = 1.0;
  int mls_polynomial_order_ = 3;
  float mesh_mu_ = 2.5;
  float mesh_search_radius_ = 1.0;
  int mesh_max_neighbours_ = 25;
  int mesh_surf_angle_ = 150;
  int mesh_min_angle_ = 5;
  int mesh_max_angle_ = 85;

  // Element Type
  // nElType = 1    -> C3D8   Hexahedron
  // nElType = 2    -> C3D6   Triangular prism
  int nElType;

  // Debug mode
  bool bDebugMode = false;

  int ifmin_ = 0;
  int ifmax_ = 0;
  float vfmin_ = 0.0;
  float vfmax_ = 0.0;
  float vamin_ = 0.0;
  float vamax_ = 0.0;
  float emin_ = 0.0;
  float emax_ = 0.0;

  std::vector<float> ve_;

};  // class FEA2

#endif  // FEA2_H
