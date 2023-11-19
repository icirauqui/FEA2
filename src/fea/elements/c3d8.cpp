#include "c3d8.hpp"



Eigen::VectorXd c3d8::computeShapeFunctions(double xi, double eta, double zeta) {
    // Initialize a vector to store the shape functions
    Eigen::VectorXd N(8);

    // Compute the shape functions for each node
    N(0) = (1 - xi) * (1 - eta) * (1 - zeta) / 8;
    N(1) = (1 + xi) * (1 - eta) * (1 - zeta) / 8;
    N(2) = (1 + xi) * (1 + eta) * (1 - zeta) / 8;
    N(3) = (1 - xi) * (1 + eta) * (1 - zeta) / 8;
    N(4) = (1 - xi) * (1 - eta) * (1 + zeta) / 8;
    N(5) = (1 + xi) * (1 - eta) * (1 + zeta) / 8;
    N(6) = (1 + xi) * (1 + eta) * (1 + zeta) / 8;
    N(7) = (1 - xi) * (1 + eta) * (1 + zeta) / 8;    

    return N;
}

Eigen::MatrixXd c3d8::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
    Eigen::MatrixXd dN(8, 3); // 6 nodes, 3 natural coordinates (xi, eta, zeta)

    // Derivatives of shape functions with respect to xi, eta, zeta
    dN(0, 0) = -(1 - eta) * (1 - zeta) / 8; // dN1/dxi
    dN(0, 1) = -(1 - xi) * (1 - zeta) / 8; // dN1/deta
    dN(0, 2) = -(1 - xi) * (1 - eta) / 8; // dN1/dzeta

    dN(1, 0) =  (1 - eta) * (1 - zeta) / 8; // dN2/dxi
    dN(1, 1) = -(1 + xi) * (1 - zeta) / 8; // dN2/deta
    dN(1, 2) = -(1 + xi) * (1 - eta) / 8; // dN2/dzeta

    dN(2, 0) =  (1 + eta) * (1 - zeta) / 8; // dN3/dxi
    dN(2, 1) =  (1 + xi) * (1 - zeta) / 8; // dN3/deta
    dN(2, 2) = -(1 + xi) * (1 + eta) / 8; // dN3/dzeta

    dN(3, 0) = -(1 + eta) * (1 - zeta) / 8; // dN4/dxi
    dN(3, 1) =  (1 - xi) * (1 - zeta) / 8; // dN4/deta
    dN(3, 2) = -(1 - xi) * (1 + eta) / 8; // dN4/dzeta

    dN(4, 0) = -(1 - eta) * (1 + zeta) / 8; // dN5/dxi
    dN(4, 1) = -(1 - xi) * (1 + zeta) / 8; // dN5/deta
    dN(4, 2) =  (1 - xi) * (1 - eta) / 8; // dN5/dzeta

    dN(5, 0) =  (1 - eta) * (1 + zeta) / 8; // dN6/dxi
    dN(5, 1) = -(1 + xi) * (1 + zeta) / 8; // dN6/deta
    dN(5, 2) =  (1 + xi) * (1 - eta) / 8; // dN6/dzeta

    dN(6, 0) =  (1 + eta) * (1 + zeta) / 8; // dN7/dxi
    dN(6, 1) =  (1 + xi) * (1 + zeta) / 8; // dN7/deta
    dN(6, 2) =  (1 + xi) * (1 + eta) / 8; // dN7/dzeta

    dN(7, 0) = -(1 + eta) * (1 + zeta) / 8; // dN8/dxi
    dN(7, 1) =  (1 - xi) * (1 + zeta) / 8; // dN8/deta
    dN(7, 2) =  (1 - xi) * (1 + eta) / 8; // dN8/dzeta

    return dN;
}

Eigen::MatrixXd c3d8::computeJacobian(const std::array<Eigen::Vector3d, 8>& nodes, double xi, double eta, double zeta) {
    // Compute the derivatives of the shape functions
    Eigen::MatrixXd dN = c3d8::computeShapeFunctionDerivatives(xi, eta, zeta);

    // Initialize the Jacobian matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, 3);

    // Compute the Jacobian matrix
    for (int i = 0; i < 8; ++i) {
        J.col(0) += nodes[i] * dN(i, 0); // Contribution to the first column of J
        J.col(1) += nodes[i] * dN(i, 1); // Contribution to the second column of J
        J.col(2) += nodes[i] * dN(i, 2); // Contribution to the third column of J
    }

    return J;
}

// Function to compute the inverse of the Jacobian matrix and its determinant
std::pair<Eigen::MatrixXd, double> c3d8::computeInverseJacobianAndDet(const Eigen::MatrixXd& J) {
    double detJ = J.determinant();
    Eigen::MatrixXd invJ = J.inverse();
    return {invJ, detJ};
}

// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd c3d8::computeStrainDisplacementMatrix(const std::array<Eigen::Vector3d, 8>& nodes, double xi, double eta, double zeta) {
    // Compute the derivatives of the shape functions
    Eigen::MatrixXd dN = c3d8::computeShapeFunctionDerivatives(xi, eta, zeta);

    // Compute the Jacobian matrix and its inverse
    Eigen::MatrixXd J = c3d8::computeJacobian(nodes, xi, eta, zeta);
    auto [invJ, detJ] = c3d8::computeInverseJacobianAndDet(J);

    // Compute the derivatives of shape functions w.r.t. global coordinates
    Eigen::MatrixXd dNdXYZ = invJ * dN.transpose();

    // Initialize the B matrix
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 24); // 6 strain components, 24 displacement components (3 per node)

    // Fill the B matrix
    for (int i = 0; i < 8; ++i) {
        B(0, 3 * i) = dNdXYZ(0, i);
        B(1, 3 * i + 1) = dNdXYZ(1, i);
        B(2, 3 * i + 2) = dNdXYZ(2, i);

        B(3, 3 * i) = dNdXYZ(1, i);
        B(3, 3 * i + 1) = dNdXYZ(0, i);

        B(4, 3 * i + 1) = dNdXYZ(2, i);
        B(4, 3 * i + 2) = dNdXYZ(1, i);

        B(5, 3 * i) = dNdXYZ(2, i);
        B(5, 3 * i + 2) = dNdXYZ(0, i);
    }

    return B;
}


Eigen::MatrixXd c3d8::computeElasticityMatrix(double E, double nu) {
    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(6, 6);
    D(0, 0) = D(1, 1) = D(2, 2) = lambda + 2 * mu;
    D(0, 1) = D(0, 2) = D(1, 0) = D(1, 2) = D(2, 0) = D(2, 1) = lambda;
    D(3, 3) = D(4, 4) = D(5, 5) = mu;

    return D;
}


// Function to compute the stiffness matrix for a triangular prism
Eigen::MatrixXd c3d8::computeStiffnessMatrix(const std::array<Eigen::Vector3d, 8>& nodes, double E, double nu) {
    // Define material properties
    Eigen::MatrixXd D = c3d8::computeElasticityMatrix(E, nu);

    // Initialize the stiffness matrix
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(24, 24); // 24x24 for 8 nodes, 3 DOF each

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

                Eigen::MatrixXd J = c3d8::computeJacobian(nodes, xi, eta, zeta);
                auto [invJ, detJ] = c3d8::computeInverseJacobianAndDet(J);
                Eigen::MatrixXd B = c3d8::computeStrainDisplacementMatrix(nodes, xi, eta, zeta);

                // Weight calculation considering different weights
                double weight = gaussWeights[i] * gaussWeights[j] * gaussWeights[k];

                // K += B^T * D * B * detJ
                K += B.transpose() * D * B * detJ;
            }
        }
    }

    return K;
}



Eigen::MatrixXd c3d8::matAssembly(std::vector<std::vector<float> > &vpts, 
                                  std::vector<std::vector<int> > &velts, 
                                  float E, float nu) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3*vpts.size(), 3*vpts.size());

  for (auto elt : velts) {
    std::array<Eigen::Vector3d, 8> xyzi;
    std::vector<int> mn;

    for (unsigned int i=0; i<elt.size(); i++) {
      std::vector<float> coords = vpts[elt[i]];
      xyzi[i] = Eigen::Vector3d(coords[0], coords[1], coords[2]);
      mn.push_back(elt[i]*3);
    }

    Eigen::MatrixXd Kei = c3d8::computeStiffnessMatrix(xyzi, E, nu);

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