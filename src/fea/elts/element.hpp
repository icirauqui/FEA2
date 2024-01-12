#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

template<size_t N>
class Element {
public:
  Element(double E, double nu) {
    _E = E;
    _nu = nu;
  }
  //virtual ~Element() = 0;

  // Function to compute shape functions for a triangular prism
  virtual Eigen::VectorXd computeShapeFunctions(double xi, double eta, double zeta) = 0;

  // Function to compute the derivative of shape functions for a triangular prism
  virtual Eigen::MatrixXd computeShapeFunctionDerivatives(double xi, double eta, double zeta) = 0;

  // Function to compute the Jacobian matrix for a triangular prism
  virtual Eigen::MatrixXd computeJacobian(const std::array<Eigen::Vector3d, N>& nodes, double xi, double eta, double zeta) = 0;

  // Function to compute the inverse of the Jacobian matrix and its determinant
  virtual std::pair<Eigen::MatrixXd, double> computeInverseJacobianAndDet(const Eigen::MatrixXd& J) = 0;

  // Function to compute the Strain-Displacement Matrix (B)
  virtual Eigen::MatrixXd computeStrainDisplacementMatrix(const std::array<Eigen::Vector3d, N>& nodes, double xi, double eta, double zeta) = 0;

  // Function to compute the elasticity matrix
  virtual void computeElasticityMatrix() = 0;

  // Function to compute the stiffness matrix for a triangular prism
  virtual Eigen::MatrixXd computeStiffnessMatrix(const std::array<Eigen::Vector3d, N>& nodes) = 0;

  // Function to assemble the global stiffness matrix
  virtual Eigen::MatrixXd matAssembly(std::vector<Eigen::Vector3d> &vpts, 
                                      std::vector<std::array<unsigned int, N>> &velts) = 0;


  // Accessors
  std::string getElementCode() const { return _element_name; }
  double getE() const { return _E; }
  double getNu() const { return _nu; }
  Eigen::MatrixXd getD() const { return _D; }

protected:
  std::string _element_name;
  double _E;
  double _nu;
  Eigen::MatrixXd _D;
};


#endif // ELEMENT_HPP