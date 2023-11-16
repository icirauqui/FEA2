#ifndef C3D6_HPP
#define C3D6_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

namespace c3d6 {

    // Function to compute shape functions for a triangular prism
    Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta);

    // Function to compute the derivative of shape functions for a triangular prism
    Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta);

    // Function to compute the Jacobian matrix for a triangular prism
    Eigen::MatrixXd computeJacobian(const std::array<Eigen::Vector3d, 6>& nodes, double xi, double eta, double zeta);

    // Function to compute the inverse of the Jacobian matrix and its determinant
    std::pair<Eigen::MatrixXd, double> computeInverseJacobianAndDet(const Eigen::MatrixXd& J);

    // Function to compute the Strain-Displacement Matrix (B)
    Eigen::MatrixXd computeStrainDisplacementMatrix(const std::array<Eigen::Vector3d, 6>& nodes, double xi, double eta, double zeta);

    // Function to compute the elasticity matrix
    Eigen::MatrixXd computeElasticityMatrix(double E, double nu);

    // Function to compute the stiffness matrix for a triangular prism
    Eigen::MatrixXd computeStiffnessMatrix(const std::array<Eigen::Vector3d, 6>& nodes, double E, double nu);

} // namespace c3d6

#endif // C3D6_HPP