#ifndef C2D4_HPP
#define C2D4_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element2d.hpp"

class C2D4 : public Element2D {

public:

    C2D4(double E, double nu);

    // Function to compute shape functions for a triangular prism
    Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta);

    // Function to compute the derivative of shape functions for a triangular prism
    Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta);

    // Function to compute the Jacobian matrix for a triangular prism
    Eigen::MatrixXd computeJacobian(const std::vector<Eigen::Vector3d>& nodes, Eigen::MatrixXd &dN);

    // Function to compute the inverse of the Jacobian matrix and its determinant
    std::pair<Eigen::MatrixXd, double> computeInverseJacobianAndDet(const Eigen::MatrixXd& J);

    // Function to compute the Strain-Displacement Matrix (B)
    Eigen::MatrixXd computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ);

    // Function to compute the stiffness matrix for a triangular prism
    Eigen::MatrixXd computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes);

    // Function to assemble the global stiffness matrix
    Eigen::MatrixXd matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                std::vector<std::vector<unsigned int>> &velts);

}; // class C2D4

#endif // C2D4_HPP