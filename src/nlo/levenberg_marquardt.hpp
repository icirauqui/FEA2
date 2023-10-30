#ifndef NLO_HPP
#define NLO_HPP

#include <eigen3/Eigen/Dense>
#include "../fea/fea.hpp"
#include "../fea/pos.hpp"

//#include "node.hpp"
//#include "edge.hpp"


class LevenbergMarquardt 
{

public:

  //LevenbergMarquardt(
  //  POS* pos, FEA* fea,
  //  int maxIterations, double lambda, 
  //  double tolerance, double dampingFactor);

  LevenbergMarquardt(
    POS* pos, FEA* fea,
    int maxIterations = 1000, double lambda = 0.001, 
    double tolerance = 0.000001, double damping_factor = 2.0);

  std::pair<int, Eigen::VectorXd> Optimize(int max_iters, const Eigen::VectorXd params0);

  std::pair<int, Eigen::MatrixXd> OptimizeStep(const Eigen::VectorXd& params0);

  double GetResidual();

  //void Optimize();

private:

  double ComputeResidual(const Eigen::VectorXd& params);
  Eigen::MatrixXd ComputeJacobian(const Eigen::VectorXd& params, double residual_original, 
                                  double delta_t = 0.00001, double delta_q = 0.000001, double delta_s = 0.00001);

  POS* pos_;
  FEA* fea_;

  double residual_;

  int max_iterations_;
  double lambda_;
  double tolerance_;
  double damping_factor_;




};


#endif