#include "c3d8.hpp"


C3D8::C3D8(double E, double nu) : Element3D(E, nu) {
  _element_name = "C3D8";
  _num_nodes = 8;
  _dof_per_node = 3;
}

Eigen::VectorXd C3D8::computeShapeFunctions(double xi, double eta, double zeta) {
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

Eigen::MatrixXd C3D8::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
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

Eigen::MatrixXd C3D8::computeJacobian(const std::vector<Eigen::Vector3d>& nodes, double xi, double eta, double zeta) {
    // Compute the derivatives of the shape functions
    Eigen::MatrixXd dN = C3D8::computeShapeFunctionDerivatives(xi, eta, zeta);

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
std::pair<Eigen::MatrixXd, double> C3D8::computeInverseJacobianAndDet(const Eigen::MatrixXd& J) {
    double detJ = J.determinant();
    Eigen::MatrixXd invJ = J.inverse();
    return {invJ, detJ};
}

// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd C3D8::computeStrainDisplacementMatrix(const std::vector<Eigen::Vector3d>& nodes, double xi, double eta, double zeta) {
    // Compute the derivatives of the shape functions
    Eigen::MatrixXd dN = C3D8::computeShapeFunctionDerivatives(xi, eta, zeta);

    // Compute the Jacobian matrix and its inverse
    Eigen::MatrixXd J = C3D8::computeJacobian(nodes, xi, eta, zeta);
    auto [invJ, detJ] = C3D8::computeInverseJacobianAndDet(J);

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


// Function to compute the stiffness matrix for a triangular prism
Eigen::MatrixXd C3D8::computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes) {
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

                Eigen::MatrixXd J = C3D8::computeJacobian(nodes, xi, eta, zeta);
                auto [invJ, detJ] = C3D8::computeInverseJacobianAndDet(J);
                Eigen::MatrixXd B = C3D8::computeStrainDisplacementMatrix(nodes, xi, eta, zeta);

                // Weight calculation considering different weights
                double weight = gaussWeights[i] * gaussWeights[j] * gaussWeights[k];

                // K += B^T * D * B * detJ
                K += B.transpose() * _D * B * detJ;
            }
        }
    }

    return K;
}


// template a function to print std::vector
template <typename T>
void printVector(std::string title, std::vector<T> &v) {
  std::cout << title << ":";
  for (auto i : v) {
    std::cout << " " << i;
  }
  std::cout << std::endl;
}


Eigen::MatrixXd C3D8::matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                  std::vector<std::vector<unsigned int>> &velts) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3*vpts.size(), 3*vpts.size());

  int num_element = 0;
  for (auto elt : velts) {

    std::vector<Eigen::Vector3d> xyzi(8);
    std::vector<int> mn(8);


    for (unsigned int i=0; i<elt.size(); i++) {
      xyzi[i] = vpts[elt[i]];
      mn[i] = elt[i]*_dof_per_node;
    }

    //printVector("Element " + std::to_string(num_element), elt);
    //printVector("       mn", mn);
    //printVector("xyzi", xyzi);
    
    num_element++;

    Eigen::MatrixXd Kei = C3D8::computeStiffnessMatrix(xyzi);

    // For each node (8) in the element
    for (unsigned int ni = 0; ni < mn.size(); ni++) {
      for (unsigned int nj = 0; nj < mn.size(); nj++) {
        for (unsigned int m = 0; m < _dof_per_node; m++) {
          for (unsigned int n = 0; n < _dof_per_node; n++) {
            K(mn[ni]+m, mn[nj]+n) += Kei(ni*3+m, nj*3+n);
          }
        }
      }
    }
  }

  return K;
}
