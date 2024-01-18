#ifndef C3D6_HPP
#define C3D6_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element3d.hpp"

class C3D6 : public Element3D {

public:
    
    C3D6(double E, double nu);

    // Function to compute shape functions for a triangular prism
    Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta);

    // Function to compute the derivative of shape functions for a triangular prism
    Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta);

    // Function to compute the Jacobian matrix for a triangular prism
    Eigen::MatrixXd computeJacobian(const std::vector<Eigen::Vector3d>& nodes, double xi, double eta, double zeta);

    // Function to compute the inverse of the Jacobian matrix and its determinant
    std::pair<Eigen::MatrixXd, double> computeInverseJacobianAndDet(const Eigen::MatrixXd& J);

    // Function to compute the Strain-Displacement Matrix (B)
    Eigen::MatrixXd computeStrainDisplacementMatrix(const std::vector<Eigen::Vector3d>& nodes, double xi, double eta, double zeta);

    // Function to compute the stiffness matrix for a triangular prism
    Eigen::MatrixXd computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes);

    // Function to assemble the global stiffness matrix
    Eigen::MatrixXd matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                std::vector<std::vector<unsigned int>> &velts);

}; // class c3d6

#endif // C3D6_HPP