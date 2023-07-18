#include "fem.hpp"

FEM::FEM(std::string element): element_(element) {}

void FEM::AddPoint(Eigen::Vector3d point) {
  points_.push_back(point);
  points_alive_.push_back(true);
}


bool FEM::InitCloud() {
  pc_.width = points_.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = points_[i](0);
    pc_.points[i].y = points_[i](1);
    pc_.points[i].z = points_[i](2);
  }

  return pc_.width > 0;
}

bool FEM::MovingLeastSquares() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

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

  // Replace pc_ with mls_points
  pc_.width = mls_points.size();
  pc_.height = 1;
  pc_.points.resize(pc_.width * pc_.height);

  for (size_t i = 0; i < pc_.points.size(); ++i) {
    pc_.points[i].x = mls_points[i].x;
    pc_.points[i].y = mls_points[i].y;
    pc_.points[i].z = mls_points[i].z;
  }

  return true;
}


bool FEM::Triangulate() {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
    new pcl::PointCloud<pcl::PointXYZ> (pc_));

  // Normal estimation
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud (cloud);
  n.setInputCloud(cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);

  // Store normals in normals_ as Eigen::Vector3d
  normals_.clear();
  for (unsigned int i = 0; i < normals->size(); i++) {
    Eigen::Vector3d normal;
    normal << normals->points[i].normal_x,
              normals->points[i].normal_y,
              normals->points[i].normal_z;
    normals_.push_back(normal);
  }

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



bool FEM::Compute(bool moving_least_squares) {
  bool ok = InitCloud();

  if (ok && moving_least_squares) {
    ok = MovingLeastSquares();
  }

  if (ok)
    ok = Triangulate();

  return ok;
}

void FEM::ComputeExtrusion() {
  std::vector<double> distances;
  for (std::vector<int> triangle : triangles_) {
    pcl::PointXYZ p0 = pc_.points[triangle[0]];
    pcl::PointXYZ p1 = pc_.points[triangle[1]];
    pcl::PointXYZ p2 = pc_.points[triangle[2]];

    double d01 = pcl::euclideanDistance(p0, p1);
    double d12 = pcl::euclideanDistance(p1, p2);
    double d20 = pcl::euclideanDistance(p2, p0);

    distances.push_back(d01);
    distances.push_back(d12);
    distances.push_back(d20);
  }

  // Get median
  std::sort(distances.begin(), distances.end());
  element_height_ = distances[distances.size()/2];
  
  // Compute second layer at a distance of 1/2 element height and with the direction of the normal vector
  points2_.clear();
  pc2_.width = points_.size();
  pc2_.height = 1;
  pc2_.points.resize(pc2_.width * pc2_.height);

  // Average normal
  Eigen::Vector3d normal;
  for (Eigen::Vector3d n : normals_) {
    normal -= n;
  }
  normal /= normals_.size();
  
  for (unsigned int i=0; i<points_.size(); i++) {
    Eigen::Vector3d point2 = points_[i] - element_height_/2 * normal;
    points2_.push_back(point2);

    pc2_.points[i].x = point2(0);
    pc2_.points[i].y = point2(1);
    pc2_.points[i].z = point2(2);
  }
}


void FEM::ViewMesh(bool extrusion, 
                   pcl::PointCloud<pcl::PointXYZ> cloud2,
                   Eigen::Vector3d pose1,
                   Eigen::Vector3d pose2) {
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ> (pc_));

  pcl::visualization::PCLVisualizer viewer;
  viewer.setBackgroundColor(1, 1, 1);
  viewer.addPointCloud(cloud);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2);
  viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 0.0, 0.0);
  
  for (unsigned int i=0; i<triangles_.size(); i++) {
    pcl::PointXYZ p0 = cloud->points[triangles_[i][0]];
    pcl::PointXYZ p1 = cloud->points[triangles_[i][1]];
    pcl::PointXYZ p2 = cloud->points[triangles_[i][2]];
    std::string name = "triangle_1_" + std::to_string(i);

    // Add line
    viewer.addLine<pcl::PointXYZ>(p0, p1, 0.0, 0.7, 0.0, name + "a");
    viewer.addLine<pcl::PointXYZ>(p1, p2, 0.0, 0.7, 0.0, name + "b");
    viewer.addLine<pcl::PointXYZ>(p2, p0, 0.0, 0.7, 0.0, name + "c");
  }

  if (extrusion) {
    for (unsigned int i=0; i<triangles_.size(); i++) {
      Eigen::Vector3d p0e = points2_[triangles_[i][0]];
      Eigen::Vector3d p1e = points2_[triangles_[i][1]];
      Eigen::Vector3d p2e = points2_[triangles_[i][2]];

      pcl::PointXYZ p0(p0e(0), p0e(1), p0e(2));
      pcl::PointXYZ p1(p1e(0), p1e(1), p1e(2));
      pcl::PointXYZ p2(p2e(0), p2e(1), p2e(2));
      std::string name = "triangle_1_" + std::to_string(i) + "_extrusion";

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name + "_a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.9, 0.7, name + "_b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.9, 0.7, name + "_c");
    }

    for (unsigned int i=0; i<points2_.size(); i++) {
      pcl::PointXYZ p0 = cloud->points[i];
      pcl::PointXYZ p1(points2_[i](0), points2_[i](1), points2_[i](2));
      std::string name = "line_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.9, 0.7, name);
    }
  }

  // If pose1 not zero, then add as a thick point
  if (pose1.norm() > 0) {
    pcl::PointXYZ p1(pose1(0), pose1(1), pose1(2));
    viewer.addSphere(p1, 0.1, 0.0, 0.7, 0.0, "pose1");
  }

  if (cloud2.points.size() > 0) {
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> cloud2_color(cloud2.makeShared(), 0, 0, 255);
    viewer.addPointCloud(cloud2.makeShared(), cloud2_color, "cloud2");
    viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud2");

    for (unsigned int i=0; i<triangles_.size(); i++) {
      pcl::PointXYZ p0 = cloud2.points[triangles_[i][0]];
      pcl::PointXYZ p1 = cloud2.points[triangles_[i][1]];
      pcl::PointXYZ p2 = cloud2.points[triangles_[i][2]];
      std::string name = "triangle_2_" + std::to_string(i);

      // Add line
      viewer.addLine<pcl::PointXYZ>(p0, p1, 0.7, 0.0, 0.0, name + "a");
      viewer.addLine<pcl::PointXYZ>(p1, p2, 0.7, 0.0, 0.0, name + "b");
      viewer.addLine<pcl::PointXYZ>(p2, p0, 0.7, 0.0, 0.0, name + "c");
    }

    // If pose2 not zero, then add as a thick point
    if (pose2.norm() > 0) {
      pcl::PointXYZ p2(pose2(0), pose2(1), pose2(2));
      viewer.addSphere(p2, 0.1, 0.7, 0.0, 0.0, "pose2");
    }
  }

  viewer.spin();
  viewer.close();
}


std::vector<std::vector<float>> FEM::GetNodes() {
  std::vector<std::vector<float>> points;
  for (Eigen::Vector3d pt: points_) {
    std::vector<float> point;
    point.push_back(pt[0]);
    point.push_back(pt[1]);
    point.push_back(pt[2]);
    points.push_back(point);
  }
  return points;
}


std::vector<Eigen::Vector3d> FEM::GetEigenNodes() {
  return points_;
}


std::vector<std::vector<int>> FEM::GetElements() {
  return triangles_;
}

pcl::PointCloud<pcl::PointXYZ> FEM::GetCloud() {
  return pc_;
}