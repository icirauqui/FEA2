#include "c2d4.hpp"


C2D4::C2D4(double E, double nu) : Element2D(E, nu) {
  _element_name = "C2D4";
  _num_nodes = 4;
  _dof_per_node = 2;
}

Eigen::VectorXd C2D4::computeShapeFunctions(double xi, double eta, double zeta) {
    // Initialize a vector to store the shape functions
    Eigen::VectorXd N(_num_nodes);

    // Compute the shape functions for each node
    N(0) = (1 - xi) * (1 - eta) / 4;
    N(1) = (1 + xi) * (1 - eta) / 4;
    N(2) = (1 + xi) * (1 + eta) / 4;
    N(3) = (1 - xi) * (1 + eta) / 4; 

    return N;
}

Eigen::MatrixXd C2D4::computeShapeFunctionDerivatives(double xi, double eta, double zeta) {
    Eigen::MatrixXd dN(4, 2); // 4 nodes, 2 natural coordinates (xi, eta)

    // Derivatives of shape functions with respect to xi, eta, zeta
    dN(0, 0) = -(1 - eta) / 4; // dN1/dxi
    dN(0, 1) = -(1 - xi)  / 4; // dN1/deta

    dN(1, 0) =  (1 - eta) / 4; // dN2/dxi
    dN(1, 1) = -(1 + xi)  / 4; // dN2/deta

    dN(2, 0) =  (1 + eta) / 4; // dN3/dxi
    dN(2, 1) =  (1 + xi)  / 4; // dN3/deta

    dN(3, 0) = -(1 + eta) / 4; // dN4/dxi
    dN(3, 1) =  (1 - xi)  / 4; // dN4/deta

    return dN;
}

Eigen::MatrixXd C2D4::computeJacobian(const std::vector<Eigen::Vector3d>& nodes, Eigen::MatrixXd &dN) {
    // Initialize the Jacobian matrix
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);

    // Compute the Jacobian matrix
    for (int i = 0; i < 4; ++i) {
      Eigen::Vector2d node = nodes[i].head(2);
      J.col(0) += node * dN(i, 0); // Contribution to the first column of J
      J.col(1) += node * dN(i, 1); // Contribution to the second column of J
    }

    return J;
}

// Function to compute the inverse of the Jacobian matrix and its determinant
std::pair<Eigen::MatrixXd, double> C2D4::computeInverseJacobianAndDet(const Eigen::MatrixXd& J) {
    double detJ = J.determinant();
    Eigen::MatrixXd invJ = J.inverse();
    return {invJ, detJ};
}

// Function to compute the Strain-Displacement Matrix (B)
Eigen::MatrixXd C2D4::computeStrainDisplacementMatrix(Eigen::MatrixXd &dN, Eigen::MatrixXd &invJ, double &detJ) {
    // Compute the derivatives of shape functions w.r.t. global coordinates
    Eigen::MatrixXd dNdXYZ = invJ * dN.transpose();

    // Initialize the B matrix
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8); // 3 strain components, 8 displacement components (2 per node)

    // Fill the B matrix
    for (int i = 0; i < _num_nodes; ++i) {
        B(0, _dof_per_node * i)     = dNdXYZ(0, i);
        B(1, _dof_per_node * i + 1) = dNdXYZ(1, i);

        B(2, _dof_per_node * i)     = dNdXYZ(1, i);
        B(2, _dof_per_node * i + 1) = dNdXYZ(0, i);
    }

    return B;
}


// Function to compute the stiffness matrix for a triangular prism
Eigen::MatrixXd C2D4::computeStiffnessMatrix(const std::vector<Eigen::Vector3d>& nodes) {
  // Initialize the stiffness matrix: 24x24 for 4 nodes, 2 DOF each
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(_num_nodes*_dof_per_node, _num_nodes*_dof_per_node);

  // Gauss quadrature points and weights (2-point quadrature)
  std::array<double, 2> gaussPoints = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  std::array<double, 2> gaussWeights = {1.0, 1.0};

  // Integration (simplified example, replace with appropriate Gauss quadrature in real applications)
  for (int i = 0; i < gaussPoints.size(); ++i) {
    for (int j = 0; j < gaussPoints.size(); ++j) {
      double xi = gaussPoints[i];
      double eta = gaussPoints[j];
      double zeta = 0.0;

      Eigen::MatrixXd dN = C2D4::computeShapeFunctionDerivatives(xi, eta, zeta);
                
      Eigen::MatrixXd J = C2D4::computeJacobian(nodes, dN);
      auto [invJ, detJ] = C2D4::computeInverseJacobianAndDet(J);
      Eigen::MatrixXd B = C2D4::computeStrainDisplacementMatrix(dN, invJ, detJ);

      // Weight calculation considering different weights
      double weight = gaussWeights[i] * gaussWeights[j];

      // K += B^T * D * B * detJ
      K += B.transpose() * _D * B * detJ;
    }
  }

  return K;
}


Eigen::MatrixXd C2D4::matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                  std::vector<std::vector<unsigned int>> &velts) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(_dof_per_node*vpts.size(), _dof_per_node*vpts.size());

  for (auto elt : velts) {

    std::vector<Eigen::Vector3d> xyzi(_num_nodes);
    std::vector<int> mn(_num_nodes);

    for (unsigned int i=0; i<elt.size(); i++) {
      xyzi[i] = vpts[elt[i]];
      mn[i] = elt[i]*_dof_per_node;
    }

    Eigen::MatrixXd Kei = C2D4::computeStiffnessMatrix(xyzi);

    for (unsigned int ni = 0; ni < mn.size(); ni++) {
      for (unsigned int nj = 0; nj < mn.size(); nj++) {
        for (unsigned int m = 0; m < _dof_per_node; m++) {
          for (unsigned int n = 0; n < _dof_per_node; n++) {
            K(mn[ni]+m, mn[nj]+n) += Kei(ni*_dof_per_node+m, nj*_dof_per_node+n);
          }
        }
      }
    }
  }

  return K;
}
