#include "c3d6.hpp"



C3D6::C3D6(double E, double nu) : Element3D(E, nu) {
  _element_name = "C3D6";
}


Eigen::VectorXd C3D6::computeShapeFunctions(double xi, double eta, double zeta) {
    // Initialize a vector to store the shape functions
    Eigen::VectorXd N(6);

    // Compute the shape functions for each node
    N(0) = (1 - xi - eta) * (1 - zeta) / 2;
    N(1) = xi * (1 - zeta) / 2;
    N(2) = eta * (1 - zeta) / 2;
    N(3) = (1 - xi - eta) * (1 + zeta) / 2;
    N(4) = xi * (1 + zeta) / 2;
    N(5) = eta * (1 + zeta) / 2;

    return N;
}

Eigen::MatrixXd C3D6::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
    Eigen::MatrixXd dN(6, 3); // 6 nodes, 3 natural coordinates (xi, eta, zeta)

    // Derivatives of shape functions with respect to xi, eta, zeta
    dN(0, 0) = -(1 - zeta) / 2; // dN1/dxi
    dN(0, 1) = -(1 - zeta) / 2; // dN1/deta
    dN(0, 2) = -(1 - xi - eta) / 2; // dN1/dzeta

    dN(1, 0) = (1 - zeta) / 2;  // dN2/dxi
    dN(1, 1) = 0;               // dN2/deta
    dN(1, 2) = -xi / 2;         // dN2/dzeta

    dN(2, 0) = 0;               // dN3/dxi
    dN(2, 1) = (1 - zeta) / 2;  // dN3/deta
    dN(2, 2) = -eta / 2;        // dN3/dzeta

    dN(3, 0) = -(1 + zeta) / 2; // dN4/dxi
    dN(3, 1) = -(1 + zeta) / 2; // dN4/deta
    dN(3, 2) = (1 - xi - eta) / 2; // dN4/dzeta

    dN(4, 0) = (1 + zeta) / 2;  // dN5/dxi
    dN(4, 1) = 0;               // dN5/deta
    dN(4, 2) = xi / 2;          // dN5/dzeta

    dN(5, 0) = 0;               // dN6/dxi
    dN(5, 1) = (1 + zeta) / 2;  // dN6/deta
    dN(5, 2) = eta / 2;         // dN6/dzeta

    return dN;
}

Eigen::MatrixXd C3D6::computeJacobian(const std::vector<Eigen::Vector3d>& nodes, double xi, double eta, double zeta) {
    // Compute the derivatives of the shape functions
    Eigen::MatrixXd dN = C3D6::computeShapeFunctionDerivatives(xi, eta, zeta);

    // Initialize the Jacobian matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, 3);

    // Compute the Jacobian matrix
    for (int i = 0; i < 6; ++i) {
        J.col(0) += nodes[i] * dN(i, 0); // Contribution to the first column of J
        J.col(1) += nodes[i] * dN(i, 1); // Contribution to the second column of J
        J.col(2) += nodes[i] * dN(i, 2); // Contribution to the third column of J
    }

    return J;
}

// Function to compute the inverse of the Jacobian matrix and its determinant
std::pair<Eigen::MatrixXd, double> C3D6::computeInverseJacobianAndDet(const Eigen::MatrixXd& J) {
    double detJ = J.determinant();
    Eigen::MatrixXd invJ = J.inverse();
    return {invJ, detJ};
}

// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd C3D6::computeStrainDisplacementMatrix(const std::vector<Eigen::Vector3d>& nodes, double xi, double eta, double zeta) {
    // Compute the derivatives of the shape functions
    Eigen::MatrixXd dN = C3D6::computeShapeFunctionDerivatives(xi, eta, zeta);

    // Compute the Jacobian matrix and its inverse
    Eigen::MatrixXd J = C3D6::computeJacobian(nodes, xi, eta, zeta);
    auto [invJ, detJ] = C3D6::computeInverseJacobianAndDet(J);

    // Compute the derivatives of shape functions w.r.t. global coordinates
    Eigen::MatrixXd dNdXYZ = invJ * dN.transpose();

    // Initialize the B matrix
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 18); // 6 strain components, 18 displacement components (3 per node)

    // Fill the B matrix
    for (int i = 0; i < 6; ++i) {
        B(0, i * 3)     = dNdXYZ(0, i); // Strain εxx
        B(1, i * 3 + 1) = dNdXYZ(1, i); // Strain εyy
        B(2, i * 3 + 2) = dNdXYZ(2, i); // Strain εzz
        B(3, i * 3)     = dNdXYZ(1, i); // Shear γyx
        B(3, i * 3 + 1) = dNdXYZ(0, i);
        B(4, i * 3 + 1) = dNdXYZ(2, i); // Shear γzy
        B(4, i * 3 + 2) = dNdXYZ(1, i);
        B(5, i * 3)     = dNdXYZ(2, i); // Shear γzx
        B(5, i * 3 + 2) = dNdXYZ(0, i);
    }

    return B;
}


// Function to compute the stiffness matrix for a triangular prism
Eigen::MatrixXd C3D6::computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes) {
    // Initialize the stiffness matrix
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(18, 18); // 18x18 for 6 nodes, 3 DOF each

    // Gauss quadrature points and weights (2-point quadrature)
    std::array<double, 2> gaussPoints = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    std::array<double, 2> gaussWeights = {1.0, 1.0};

    // Integration (simplified example, replace with appropriate Gauss quadrature in real applications)
    for (int i = 0; i < gaussPoints.size(); ++i) {
        for (int j = 0; j < gaussPoints.size(); ++j) {
            for (int k = 0; k < gaussPoints.size(); ++k) {
                double xi = gaussPoints[i];
                double eta = gaussPoints[j];
                double zeta = gaussPoints[k];

                Eigen::MatrixXd J = C3D6::computeJacobian(nodes, xi, eta, zeta);
                auto [invJ, detJ] = C3D6::computeInverseJacobianAndDet(J);
                Eigen::MatrixXd B = C3D6::computeStrainDisplacementMatrix(nodes, xi, eta, zeta);

                // Weight calculation considering different weights
                double weight = gaussWeights[i] * gaussWeights[j] * gaussWeights[k];

                // K += B^T * D * B * detJ
                K += B.transpose() * _D * B * detJ;
            }
        }
    }

    return K;
}




Eigen::MatrixXd C3D6::matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                std::vector<std::vector<unsigned int>> &velts) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3*vpts.size(), 3*vpts.size());

  for (auto elt : velts) {
    std::vector<Eigen::Vector3d> xyzi(6);
    std::vector<int> mn(6);

    for (unsigned int i=0; i<elt.size(); i++) {
      xyzi[i] = vpts[elt[i]];
      mn[i] = elt[i]*3;
    }

    Eigen::MatrixXd Kei = C3D6::computeStiffnessMatrix(xyzi);

    for (unsigned int ni = 0; ni < mn.size(); ni++) {
      for (unsigned int nj = 0; nj < mn.size(); nj++) {
        for (unsigned int m = 0; m < 3; m++) {
          for (unsigned int n = 0; n < 3; n++) {
            K(mn[ni]+m, mn[nj]+n) += Kei(ni*3+m, nj*3+n);
          }
        }
      }
    }
  }

  return K;
}