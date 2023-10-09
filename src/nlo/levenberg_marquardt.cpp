#include "levenberg_marquardt.hpp"

LevenbergMarquardt::LevenbergMarquardt(
  POS* pos, FEA* fea,
  int maxIterations = 1000, double lambda = 0.001, double tolerance = 0.000001, double damping_factor = 2.0)
  : pos_(pos), fea_(fea),
    max_iterations_(maxIterations), lambda_(lambda), tolerance_(tolerance), damping_factor_(damping_factor) {}

double LevenbergMarquardt::ComputeResidual(const Eigen::VectorXd& params) {
  Eigen::Vector3d tvec(params(0), params(1), params(2));
  Eigen::Vector4d qvec(params(3), params(4), params(5), params(6));
  double scale = params(7);

  std::vector<Eigen::Vector3d> nodes_k1 = pos_->GetPoints();

  pos_->Transform(qvec, tvec, scale);
  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k = pos_->GetPose();

  std::vector<Eigen::Vector3d> nodes_k = pos_->GetPoints();
  std::vector<Eigen::Vector3d> nodes_k_front, nodes_k_back;
  for (unsigned int i=0; i<nodes_k.size(); i++) {
    if (i < nodes_k.size()/2) {
      nodes_k_front.push_back(nodes_k[i]);
    } else {
      nodes_k_back.push_back(nodes_k[i]);
    }
  }

  double sE = fea_->ComputeStrainEnergy(nodes_k1, nodes_k);

  return sE;
}

std::pair<int, Eigen::MatrixXd> LevenbergMarquardt::OptimizeStep(const Eigen::VectorXd& params0) {
  double residual = ComputeResidual(params0);
  Eigen::MatrixXd jacobian = ComputeJacobian(params0, residual);

  // Compute the approximate Hessian
  Eigen::VectorXd hessian = jacobian.transpose() * jacobian;
  
  // Compute the gradient
  double gradient = jacobian(0,0) * residual;

  // Updthe diagonal of the hessian with the damping factor
  for (int i = 0; i < hessian.size(); i++) {
    hessian(i) += lambda_ * hessian(i);
  }

  // Solve for the update (delta)
  Eigen::VectorXd delta = hessian.array().inverse() * (-gradient);

  // Update the parameters
  Eigen::VectorXd params = params0 + delta;

  // Compute the new residual
  double residual_new = ComputeResidual(params);

  // If the new residual is smaller, accept the update
  int result = 0;
  if (residual_new < residual) {
    lambda_ /= damping_factor_;
    
    if (delta.norm() < tolerance_) {
      return std::make_pair(0, params); // Stop
    } else {
      return std::make_pair(1, params); // Continue
    }
  }
  // If the new residual is larger, reject the update
  else {
    lambda_ *= damping_factor_;
    return std::make_pair(-1, params0);
  }
}

Eigen::MatrixXd  LevenbergMarquardt::ComputeJacobian(const Eigen::VectorXd& params, double residual_original) {
  double delta = 0.000001;
  Eigen::MatrixXd jacobian(1, params.size());

  for (int i = 0; i < params.size(); i++) {
    Eigen::VectorXd params_delta = params;
    params_delta(i) += delta;

    double residual_delta = ComputeResidual(params_delta);
    jacobian(0, i) = (residual_delta - residual_original) / delta;
  }

  return jacobian;
}