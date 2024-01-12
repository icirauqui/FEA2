#ifndef ELEMENT3D_HPP
#define ELEMENT3D_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element.hpp"

template<size_t N>
class Element3D : public Element<N> {
public:
  using Element<N>::_E;
  using Element<N>::_nu;
  using Element<N>::_D;

  Element3D(double E, double nu) : Element<N>(E, nu) {
    computeElasticityMatrix();
  }

  void computeElasticityMatrix() override {
    double lambda = _E * _nu / ((1 + _nu) * (1 - 2 * _nu));
    double mu = _E / (2 * (1 + _nu));

    _D = Eigen::MatrixXd::Zero(6, 6);
    _D(0, 0) = _D(1, 1) = _D(2, 2) = lambda + 2 * mu;
    _D(0, 1) = _D(0, 2) = _D(1, 0) = _D(1, 2) = _D(2, 0) = _D(2, 1) = lambda;
    _D(3, 3) = _D(4, 4) = _D(5, 5) = mu;
  }
};


#endif // ELEMENT3D_HPP