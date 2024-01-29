#include <iostream>
#include "eigen3/Eigen/Dense"

class AbaqusC3D8_1 {
public:
  AbaqusC3D8_1() {
    unsigned int n = 0;
    
    for (int x=int(_x0); x>=0; x-=int(_w)) {
      for (int z=int(_z0); z>=0; z-=int(_w)) {
        for (int y=int(_y0); y>=0; y-=int(_w), n++) {
          //std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
          _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));

          if (x==0 || y==0 || z==0) continue;

          unsigned int offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
          //std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;

          std::vector<unsigned int> element = {
            n,
            n+1,
            n+1+(int(_y0)+1),
            n+(int(_y0)+1),
            n+offset,
            n+1+offset,
            n+1+(int(_y0)+1)+offset,
            n+(int(_y0)+1)+offset
          };

          _elements.push_back(element);
        }
      }
    } 
  }

  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + displacements[n]);
    }

    
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  double _x0 = 5.0;
  double _y0 = 5.0;
  double _z0 = 25.0;
  double _w = 1.0;
};





class AbaqusC3D8_2 {
public:
  AbaqusC3D8_2() {
    unsigned int n = 0;
    
    for (int x=int(_x0); x>=0; x-=int(_w)) {
      for (int z=int(_z0); z>=0; z-=int(_w)) {
        for (int y=int(_y0); y>=0; y-=int(_w), n++) {
          //std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
          _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));

          if (x==0 || y==0 || z==0) continue;

          unsigned int offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
          //std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;

          std::vector<unsigned int> element = {
            n,
            n+1,
            n+1+(int(_y0)+1),
            n+(int(_y0)+1),
            n+offset,
            n+1+offset,
            n+1+(int(_y0)+1)+offset,
            n+(int(_y0)+1)+offset
          };

          _elements.push_back(element);
        }
      }
    } 
  }

  void ApplyDisplacements(std::vector<Eigen::Vector3d> &displacements, double scale = 1.0) {
    _nodes_deformed.clear();
    for (unsigned int n=0; n<_nodes.size(); n++) {
      _nodes_deformed.push_back(_nodes[n] + displacements[n]);
    }

    
  }
    
  std::vector<Eigen::Vector3d> _nodes, _nodes_deformed;
  std::vector<std::vector<unsigned int>> _elements;

  double _x0 = 1.0;
  double _y0 = 1.0;
  double _z0 = 5.0;
  double _w = 1.0;
};



