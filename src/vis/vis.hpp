#ifndef VIS_HPP
#define VIS_HPP

#include <pcl/common/common_headers.h>
#include <pcl/console/parse.h>
#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/visualization/pcl_visualizer.h>

#if VTK_MAJOR_VERSION > 8
#include <vtk-9.2/QVTKOpenGLNativeWidget.h>
#include <vtk-9.2/vtkRenderWindow.h>
#else
#include <vtk-7.1/QVTKWidget.h>
#include <vtk-7.1/vtkRenderWindow.h>
#endif




class PCLViewer {

public:

  PCLViewer();

  void initializeViewer();

  void AddNodes(std::vector<Eigen::Vector3d> &pts, std::string name, Eigen::Vector3d color);
  void AddEdges(std::vector<std::vector<unsigned int>> &elts, std::string name, Eigen::Vector3d color);
  void AddLoads(std::vector<bool> &nodes, std::vector<std::vector<double>> &mag, double scale = 1.0);

  void Render();

private:
  
    pcl::visualization::PCLVisualizer::Ptr viewer_;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_;

};





#endif