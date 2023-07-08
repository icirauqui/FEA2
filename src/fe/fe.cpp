#include "fe.hpp"






FE::FE(std::string element): element_(element) {
  std::cout << "FE constructor" << std::endl;
}

void FE::AddPoint(Eigen::Vector3d point) {
  points_.push_back(point);
}


bool FE::InitCloud() {
  pc0_.width = points_.size();
  pc0_.height = 1;
  pc0_.points.resize(pc0_.width * pc0_.height);

  for (size_t i = 0; i < pc0_.points.size(); ++i) {
    pc0_.points[i].x = points_[i](0);
    pc0_.points[i].y = points_[i](1);
    pc0_.points[i].z = points_[i](2);
  }

  return pc0_.width > 0;
}


bool FE::MovingLeastSquares() {
  std::cout << "MovingLeastSquares" << std::endl;
  return true;
}


bool FE::Triangulate() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc0_));

  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud(cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);

  // Concatenate the XYZ and normal fields*
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);

    // Create search tree*
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud (cloud_with_normals);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
  pcl::PolygonMesh mesh_out;

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius (500.0);

  // Set typical values for the parameters
  gp3.setMu (2.5);
  gp3.setMaximumNearestNeighbors (100);
  gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
  gp3.setMinimumAngle(M_PI/18); // 10 degrees
  gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
  gp3.setNormalConsistency(false);

  // Get result
  gp3.setInputCloud (cloud_with_normals);
  gp3.setSearchMethod (tree2);
  gp3.reconstruct (mesh_out);

  // Additional vertex information
  // To which part (cloud) does each vertex belong
  std::vector<int> parts = gp3.getPartIDs();
  // Whether the vertex status is [-1,0,1,2,3] = [NONE,FREE,FRINGE,BOUNDARY,COMPLETED]
  std::vector<int> states = gp3.getPointStates();

  // Get largest part
  std::unordered_map<int, int> map;
  int max_val = -1;
  int max_key = -1;
  for (unsigned int i=0; i<parts.size(); i++) {
    if (map.find(parts[i]) == map.end()) {
      map[parts[i]] = 1;
      if (max_val < 1) {
        max_val = 1;
        max_key = parts[i];
      }
    }
    else {
      map[parts[i]]++;
      if (map[parts[i]] > max_val) {
        max_val = map[parts[i]];
        max_key = parts[i];
      }
    }
  }

  // Get triangles of largest part only
  triangles_.clear();
  for (unsigned int i=0; i<mesh_out.polygons.size(); i++) {
    unsigned int nver0 = mesh_out.polygons[i].vertices[0];
    unsigned int nver1 = mesh_out.polygons[i].vertices[1];
    unsigned int nver2 = mesh_out.polygons[i].vertices[2];

    if (parts[nver0] == max_key && parts[nver1] == max_key && parts[nver2] == max_key) {
      std::vector<int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_.push_back(triangle);
    }
  }
  
  return triangles_.size() > 0 ? true : false;
}


bool FE::Compute(bool moving_least_squares) {
  std::cout << "FE Compute" << std::endl;

  bool ok = true;

  if (ok) 
    ok = InitCloud();

  if (ok && moving_least_squares) {
    ok = MovingLeastSquares();
  }

  if (ok)
    ok = Triangulate();

  return ok;
}





