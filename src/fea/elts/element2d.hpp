#ifndef ELEMENT2D_HPP
#define ELEMENT2D_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "element.hpp"


class Element2D : public Element {
public:
  Element2D(double E, double nu) : Element(E, nu) {
    computeElasticityMatrix();
  }

  void computeElasticityMatrix() override {
    double lambda = _E * _nu / ((1 + _nu) * (1 - 2 * _nu));
    double mu = _E / (2 * (1 + _nu));

    _D = Eigen::MatrixXd::Zero(3, 3);
    _D(0, 0) = _D(1, 1) = lambda + 2 * mu;
    _D(2, 2) = mu;
    _D(0, 1) = _D(1, 0) = lambda;
  }
};


#endif // ELEMENT2D_HPP