#include <iostream>
#include "eigen3/Eigen/Dense"

class AbaqusC3D8_1 {
  public:
    AbaqusC3D8_1() {
      build_nodes();
    }
    
    ~AbaqusC3D8_1() {}

    void build_nodes() {
      unsigned int n = 0;
      
      for (int x=int(_x0); x>=0; x-=int(_w)) {
        for (int z=int(_z0); z>=0; z-=int(_w)) {
          for (int y=int(_y0); y>=0; y-=int(_w), n++) {
            //std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
            _nodes.push_back(Eigen::Vector3d(double(x),double(y),double(z)));

            if (x==0 || y==0 || z==0) continue;

            unsigned int offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
            //std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;

            std::array<unsigned int,8> element = {
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

    void build_elements() {
      unsigned int n = 0;
      unsigned int elt = 0;
      unsigned int offset = 0;
      
      for (int x=int(_x0); x>=0; x-=int(_w)) {
        for (int z=int(_z0); z>=0; z-=int(_w)) {
          for (int y=int(_y0); y>=0; y-=int(_w), n++) {
            if (x==0 || y==0 || z==0) continue;
            offset = (_y0+1)*(_z0+1);//*(_x0+1-x);
            std::cout << "n (offset) = " << n << " (" << offset << ") : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;

            std::array<unsigned int,8> element;
            element[0] = n;
            element[1] = n+1;
            element[2] = n+1+(_y0+1);
            element[3] = n+(_y0+1);
            element[4] = n+offset;
            element[5] = n+1+offset;
            element[6] = n+1+(_y0+1)+offset;
            element[7] = n+(_y0+1)+offset;
            _elements.push_back(element);
            elt++;
            //std::cout << "elt = " << elt << " : ";
            //for (int i=0; i<8; i++) {
            //  std::cout << element[i] << " ";
            //}
            //std::cout << std::endl << std::endl;

            
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



