#include <iostream>
#include "eigen3/Eigen/Dense"

class AbaqusC3D8_1 {
  public:
    AbaqusC3D8_1() {
      build_nodes();
      //build_elements();
    }
    
    ~AbaqusC3D8_1() {}

    void build_nodes() {
      unsigned int n = 0;
      
      for (int x=int(_x0); x>=0; x-=int(_w)) {
        for (int z=int(_z0); z>=0; z-=int(_w)) {
          for (int y=int(_y0); y>=0; y-=int(_w), n++) {
            std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
            _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));
          }
        }
      }  
    }

    void build_elements() {
      unsigned int n = 0;
      
      for (int x=int(_x0); x>=0; x-=int(_w)) {
        for (int z=int(_z0); z>=0; z-=int(_w)) {
          for (int y=int(_y0); y>=0; y-=int(_w), n++) {
            if (x==0 || y==0 || z==0) continue;

            std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
            std::array<unsigned int,8> element;
            //element[0] = 
            
          }
        }
      }  
    }


    std::vector<Eigen::Vector3d> _nodes;
    std::vector<std::array<unsigned int,8>> _elements;

    double _x0 = 5.0;
    double _y0 = 5.0;
    double _z0 = 25.0;
    double _w = 1.0;
};



