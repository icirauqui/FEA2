#include "levenberg_marquardt.hpp"

LevenbergMarquardt::LevenbergMarquardt(
  POS* pos, FEA* fea,
  int maxIterations, double lambda, double tolerance, double damping_factor)
  : pos_(pos), fea_(fea),
    max_iterations_(maxIterations), lambda_(lambda), tolerance_(tolerance), damping_factor_(damping_factor) {}

double LevenbergMarquardt::ComputeResidual(const Eigen::VectorXd& params) {
  Eigen::Vector3d tvec(params(0), params(1), params(2));
  Eigen::Vector4d qvec(params(3), params(4), params(5), params(6));
  double scale = params(7);

  std::vector<Eigen::Vector3d> nodes_k0 = pos_->GetTarget();
  //std::vector<Eigen::Vector3d> nodes_k1 = pos_->GetPoints();

  std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k1 = std::make_pair(qvec, tvec);
  //pos_->Transform(qvec, tvec, scale);
  std::pair<std::pair<Eigen::Vector4d, Eigen::Vector3d>, std::vector<Eigen::Vector3d>> sim = pos_->SimulateTransformToPose(pose_k1, scale);
  //std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k = pos_->GetPose();

  //std::vector<Eigen::Vector3d> nodes_k = pos_->GetPoints();
  std::vector<Eigen::Vector3d> nodes_k_front, nodes_k_back;
  for (unsigned int i=0; i<sim.second.size(); i++) {
    if (i < sim.second.size()/2) {
      nodes_k_front.push_back(sim.second[i]);
    } else {
      nodes_k_back.push_back(sim.second[i]);
    }
  }

  double sE = fea_->ComputeStrainEnergy(nodes_k0, sim.second);

  return sE;
}

std::pair<int, Eigen::MatrixXd> LevenbergMarquardt::OptimizeStep(const Eigen::VectorXd& params0) {

  //std::cout << " params0 = " << params0.transpose() << std::endl;

  double residual0 = ComputeResidual(params0);

  //std::cout << " residual0 = " << residual0 << std::endl;

  Eigen::MatrixXd jacobian = ComputeJacobian(params0, residual0, 0.0000001);

  // Compute the approximate Hessian
  Eigen::VectorXd hessian = jacobian.transpose() * jacobian;
  
  // Compute the gradient
  double gradient = jacobian(0,0) * residual0;

  // Update the diagonal of the hessian with the damping factor
  for (int i = 0; i < hessian.size(); i++) {
    hessian(i) += lambda_ * hessian(i);
  }

  // Solve for the update (delta)
  Eigen::VectorXd delta = hessian.array().inverse() * (-gradient);
  //std::cout << "   delta = " << delta.transpose() << std::endl;

  //std::cout << " delta = " << delta.transpose() << std::endl;

  // Update the parameters
  Eigen::VectorXd params = params0 + delta;

  //std::cout << " params1 = " << params.transpose() << std::endl;

  // Compute the new residual
  double residual1 = ComputeResidual(params);

  //std::cout << " residual1 = " << residual1 << std::endl;

  // If the new residual is smaller, accept the update, else reject it
  int result = 0;
  if (residual1 < residual0) {
    lambda_ /= damping_factor_;

    //std::cout << std::endl;
    //std::cout << " residual1 < residual0 = " << residual1 << " < " << residual0 << std::endl;

    pos_->TransformToPose(std::make_pair(params.segment(3,4), params.segment(0,3)), params(7));

    std::pair<Eigen::Vector4d, Eigen::Vector3d> pose_k = pos_->GetPose();
    //std::cout << "  pose_k = " << pose_k.second.transpose() << " " << pose_k.first.transpose() << std::endl;
    
    if (delta.norm() < tolerance_) {
      return std::make_pair(0, params); // Stop
    } else {
      return std::make_pair(1, params); // Continue
    }
  }
  else {
    //std::cout << " residual1 >= residual0 = " << residual1 << " >= " << residual0 << std::endl;
    lambda_ *= damping_factor_;
    return std::make_pair(-1, params0);
  }
}

Eigen::MatrixXd  LevenbergMarquardt::ComputeJacobian(const Eigen::VectorXd& params, double residual_original, double delta) {
  Eigen::MatrixXd jacobian(1, params.size());

  //std::cout << "  ComputeJacobian = ";
  for (int i = 0; i < params.size(); i++) {
    Eigen::VectorXd params_delta = params;
    params_delta(i) += delta;

    double residual_delta = ComputeResidual(params_delta);
    //std::cout << residual_delta << " ";
    jacobian(0, i) = (residual_delta - residual_original) / delta;
  }

  //std::cout << std::endl;

  return jacobian;
}