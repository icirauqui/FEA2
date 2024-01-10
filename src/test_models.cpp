#include <iostream>
#include "eigen3/Eigen/Dense"

std::vector<Eigen::Vector3d> build_abaqus_c3d8_1() {
    double x0 = 5.0;
    double y0 = 5.0;
    double z0 = 25.0;
    unsigned int n = 0;

    std::vector<Eigen::Vector3d> nodes;
    
    for (int x=int(x0); x>=0; x--) {
        for (int z=int(z0); z>=0; z--) {
            for (int y=int(y0); y>=0; y--, n++) {
                std::cout << "n = " << n << " : (x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;
                nodes.push_back(Eigen::Vector3d(x,y,z));
            }
        }
    }

    return nodes;
}