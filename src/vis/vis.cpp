#include "vis.hpp"

PCLViewer::PCLViewer() {
    initializeViewer();
}

void PCLViewer::initializeViewer() {
  viewer_.reset(new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer_->setBackgroundColor(1.0, 1.0, 1.0);
  viewer_->addCoordinateSystem(1.0);
  viewer_->initCameraParameters();
  viewer_->setSize(1070, 820);
}

void PCLViewer::AddNodes(std::vector<Eigen::Vector3d> &pts, std::string name, Eigen::Vector3d color) {
  cloud_ = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
  cloud_->width = pts.size();
  cloud_->height = 1;
  cloud_->is_dense = false;
  cloud_->points.resize(cloud_->width * cloud_->height);

  for (size_t i = 0; i < cloud_->points.size(); ++i) {
      cloud_->points[i].x = pts[i](0);
      cloud_->points[i].y = pts[i](1);
      cloud_->points[i].z = pts[i](2);
  }

  viewer_->addPointCloud<pcl::PointXYZ>(cloud_, name);
  viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, name);
  viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, color(0), color(1), color(2), name);
}

void PCLViewer::AddEdges(std::vector<std::vector<unsigned int>> &elts, std::string name, Eigen::Vector3d color) {

  // Create a vector of tuples with pairs of point indices
  std::vector<std::tuple<unsigned int, unsigned int>> edges = {
    std::make_tuple(0, 1),
    std::make_tuple(1, 2),
    std::make_tuple(2, 3),
    std::make_tuple(3, 0),
    std::make_tuple(4, 5),
    std::make_tuple(5, 6),
    std::make_tuple(6, 7),
    std::make_tuple(7, 4),
    std::make_tuple(0, 4),
    std::make_tuple(1, 5),
    std::make_tuple(2, 6),
    std::make_tuple(3, 7)
  };

  // Create a storage for edges to be created so they are not duplicated
  std::vector<std::tuple<unsigned int, unsigned int>> edges_to_create;

  for (unsigned int e=0; e<elts.size(); e++) {
    //std::string name = "edge_" + name + "_" + std::to_string(e) + "_";
    for (auto edge : edges) {
      unsigned int i0 = std::get<0>(edge);
      unsigned int i1 = std::get<1>(edge);
      unsigned int node0 = elts[e][i0];
      unsigned int node1 = elts[e][i1];
      std::tuple<unsigned int, unsigned int> edge_to_create = std::make_tuple(node0, node1);
      if (std::find(edges_to_create.begin(), edges_to_create.end(), edge_to_create) == edges_to_create.end()) {
        edges_to_create.push_back(edge_to_create);
        viewer_->addLine(cloud_->points[node0], 
                         cloud_->points[node1], 
                         color(0), color(1), color(2),
                         name + std::to_string(node0)+"_"+std::to_string(node1));
      }
    }
  }
}

void PCLViewer::AddLoads(std::vector<bool> &nodes, std::vector<std::vector<double>> &mag, double scale) {
  for (unsigned int i=0; i<nodes.size(); i++) {
    if (!nodes[i]) continue;

    if (mag[i][0] == 0.0 && mag[i][1] == 0.0 && mag[i][2] == 0.0) {
      // Draw a sphere in the node
      viewer_->addSphere<pcl::PointXYZ>(cloud_->points[i], 0.2, 1.0, 1.0, 0.0, "load_" + std::to_string(i));
    } else {
      std::cout << std::endl;
      std::cout << "pts[" << i << "]" << std::endl;
      std::cout << cloud_->points[i].x << ", " << cloud_->points[i].y << ", " << cloud_->points[i].z << std::endl;
      std::cout << mag[i][0] << ", " << mag[i][1] << ", " << mag[i][2] << std::endl;
      std::cout << cloud_->points[i].x + scale * mag[i][0] << ", " << cloud_->points[i].y + scale * mag[i][1] << ", " << cloud_->points[i].z + scale * mag[i][2] << std::endl;

      viewer_->addArrow<pcl::PointXYZ, pcl::PointXYZ>(pcl::PointXYZ(cloud_->points[i].x + scale * mag[i][0],
                                                                    cloud_->points[i].y + scale * mag[i][1],
                                                                    cloud_->points[i].z + scale * mag[i][2]),
                                                      cloud_->points[i],
                                                      1.0, 0.0, 0.0, false, "load_" + std::to_string(i));
    }

  }

}

void PCLViewer::Render() {
  viewer_->spin();
}