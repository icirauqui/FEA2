#ifndef FE_HPP
#define FE_HPP


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include <eigen3/Eigen/Dense>


// PCL Libraries
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



class FE {

public:

  FE(std::string element);

  void AddPoint(Eigen::Vector3d point);

  bool InitCloud();

  bool MovingLeastSquares();
  
  bool Triangulate();

  bool Compute(bool moving_least_squares = true);



private:

  std::string element_;

  std::vector<Eigen::Vector3d> points_;

  pcl::PointCloud<pcl::PointXYZ> pc0_;

  std::vector<std::vector<int>> triangles_;


  // Interface 
  float mls_search_radius_ = 1.0;
  int mls_polynomial_order_ = 3;
  float mesh_mu_ = 2.5;
  float mesh_search_radius_ = 1.0;
  int mesh_max_neighbours_ = 25;
  int mesh_surf_angle_ = 150;
  int mesh_min_angle_ = 5;
  int mesh_max_angle_ = 85;


};
















#endif