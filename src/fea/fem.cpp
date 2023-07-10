#include "fem.hpp"






FE::FE(std::string element): element_(element) {}

void FE::AddPoint(Eigen::Vector3d point) {
  points_.push_back(point);
  points_alive_.push_back(true);
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
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc0_));

  // Create a KD-Tree
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);

  // Output has the PointNormal type in order to store the normals calculated by MLS
  pcl::PointCloud<pcl::PointNormal> mls_points;

  // Init object (second point type is for the normals, even if unused)
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
 
  mls.setComputeNormals (true);

  // Set parameters
  mls.setInputCloud (cloud);
  mls.setPolynomialOrder (2);
  mls.setSearchMethod (tree);
  mls.setSearchRadius (2);

  // Reconstruct
  mls.process (mls_points);

  if (mls_points.size() == 0) {
    std::cout << "MLS failed" << std::endl;
    return false;
  }

  // Get corresponding indexes: for each output point, returns the index of the
  // input one.
  mls_indices_.clear();

  pcl::PointIndicesPtr pIdx1 = mls.getCorrespondingIndices();
  for (unsigned int i = 0; i < pIdx1->indices.size(); i++)
    mls_indices_.push_back(pIdx1->indices[i]);

  int idxit = 0;
  for (unsigned int i = 0; i < points_.size(); i++) {
    int currentpos = i;
    if (currentpos == mls_indices_[idxit]) {
      idxit++;
    } else {
      points_alive_[i] = false;
    }
  }

  // Replace pc0_ with mls_points
  pc0_.width = mls_points.size();
  pc0_.height = 1;
  pc0_.points.resize(pc0_.width * pc0_.height);

  for (size_t i = 0; i < pc0_.points.size(); ++i) {
    pc0_.points[i].x = mls_points[i].x;
    pc0_.points[i].y = mls_points[i].y;
    pc0_.points[i].z = mls_points[i].z;
  }

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

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius (1.5);

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
  gp3.reconstruct (mesh_);

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
  for (unsigned int i=0; i<mesh_.polygons.size(); i++) {
    unsigned int nver0 = mesh_.polygons[i].vertices[0];
    unsigned int nver1 = mesh_.polygons[i].vertices[1];
    unsigned int nver2 = mesh_.polygons[i].vertices[2];

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
  bool ok = InitCloud();

  if (ok && moving_least_squares) {
    ok = MovingLeastSquares();
  }

  if (ok)
    ok = Triangulate();

  return ok;
}


void FE::ViewMesh() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc0_));

  pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloud);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 0, 0);
  
  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloud->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloud->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloud->points[triangles_[i][2]];
    std::string name = "triangle" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 1, 0, 0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 1, 0, 0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 1, 0, 0, name + "c");
  }

  viewer.spin();
  viewer.close();
}



