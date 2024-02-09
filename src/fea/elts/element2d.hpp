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

    // [1.0-v,     v, 0.0]
    // [v,     1.0-v, 0.0]
    // [0.0,     0.0, 0.5-v]

    double mult = _E/(1.0+_nu)/(1.0-2.0*_nu);
    _D = Eigen::MatrixXd::Zero(3, 3);
    _D(0, 0) = _D(1, 1) = (1-_nu);
    _D(0, 1) = _D(1, 0) = _nu;
    _D(2, 2) = 0.5 - _nu;
    _D *= mult;

    std::cout << "D: " << std::endl << _D << std::endl;

    double L = _E * _nu / ((1 + _nu) * (1 - 2 * _nu));
    double G = _E / (2 * (1 + _nu));
    
    double Ls = 2*L*G / (L + 2*G);
    Ls = L;
    double Gs = G;
    
    _D = Eigen::MatrixXd::Zero(3, 3);
    _D(0, 0) = _D(1, 1) = Ls + 2 * Gs;
    _D(2, 2) = Gs;
    _D(0, 1) = _D(1, 0) = Ls;

    std::cout << "D: " << std::endl << _D << std::endl;
  }
};


#endif // ELEMENT2D_HPP