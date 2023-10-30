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


std::pair<int, Eigen::VectorXd> LevenbergMarquardt::Optimize(int max_iters, const Eigen::VectorXd params0) {

  Eigen::VectorXd params = params0;

  for (int it = 0; it < max_iters; it++) {
    std::pair<int, Eigen::MatrixXd> params1 = OptimizeStep(params);
    params = params1.second;
    if (params1.first == 0 || params1.first == 1) {
      std::cout << "Step " << it+1 << " / " << max_iters << " : "
                << params1.first << " : "
                //<< params1.second.transpose() << " : "
                << ComputeResidual(params1.second) << std::endl;
      if (params1.first == 0) {
        break;
      } 
    } else {
      std::cout << "Step " << it+1 << " / " << max_iters << " : "
                << params1.first << " : "
                //<< params1.second.transpose() << " : "
                << ComputeResidual(params1.second) << std::endl;
    }
  }

  return std::make_pair(1, params);
}


std::pair<int, Eigen::MatrixXd> LevenbergMarquardt::OptimizeStep(const Eigen::VectorXd& params0) {

  bool verbose = false;
  //std::cout << " params0 = " << params0.transpose() << std::endl;

  double residual0 = ComputeResidual(params0);

  Eigen::MatrixXd jacobian = ComputeJacobian(params0, residual0, 0.00001, 0.001, 0.0000001);

  // Compute the approximate Hessian
  Eigen::MatrixXd hessian = jacobian.transpose() * jacobian;
  
  // Compute the gradient
  Eigen::VectorXd gradient = jacobian.transpose() * residual0;
  
  // Update the diagonal of the Hessian with the damping factor
  for (int i = 0; i < hessian.rows(); i++) {
    hessian(i,i) += lambda_ * hessian(i,i);
  }

  // Solve for the update (delta)
  Eigen::VectorXd delta = hessian.llt().solve(-gradient);

  // Update the parameters
  Eigen::VectorXd params = params0 + delta;

  // Compute the new residual
  double residual1 = ComputeResidual(params);
  std::cout << " residual = " << residual0 << " -> " << residual1 << std::endl;

  // Compute delta norm
  double delta_norm = delta.norm();

  
  // If the new residual is smaller, accept the update, else reject it
  if (residual1 < residual0) {
    lambda_ /= damping_factor_;
    residual_ = residual1;

    pos_->TransformToPose(std::make_pair(params.segment(3,4), params.segment(0,3)), params(7));
    
    if (delta.norm() < tolerance_) {
      return std::make_pair(0, params); // Stop
    } else {
      return std::make_pair(1, params); // Continue
    }
  }
  else {
    lambda_ *= damping_factor_;
    return std::make_pair(-1, params0);
  }
  
}

double LevenbergMarquardt::GetResidual() {
  return residual_;
}

Eigen::MatrixXd  LevenbergMarquardt::ComputeJacobian(const Eigen::VectorXd& params, double residual_original, 
                                                     double delta_t, double delta_q, double delta_s) {
  Eigen::MatrixXd jacobian(1, params.size());

  Eigen::VectorXd deltas = params;
  deltas.segment(0,3) *= delta_t;
  deltas(3) = delta_q;
  deltas(4) = delta_q;
  deltas(5) = delta_q;
  deltas(6) = delta_q;
  deltas(7) = delta_s;

  for (int i = 0; i < params.size(); i++) {
    Eigen::VectorXd params_delta = params;
    params_delta(i) += deltas(i);
    double residual_delta = ComputeResidual(params_delta);
    jacobian(0, i) = (residual_delta - residual_original) / deltas(i);
  }

  return jacobian;
}