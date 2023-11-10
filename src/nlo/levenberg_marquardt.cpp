#include "levenberg_marquardt.hpp"

LevenbergMarquardt::LevenbergMarquardt(
  POS* pos, FEA* fea,
  int maxIterations, double lambda, double damping_factor, double tolerance)
  : pos_(pos), fea_(fea),
    max_iterations_(maxIterations), lambda_(lambda), damping_factor_(damping_factor), tolerance_(tolerance) {}

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


//std::pair<int, Eigen::VectorXd> LevenbergMarquardt::Optimize(const Eigen::VectorXd params0) {
//
//  Eigen::VectorXd params = params0;
//
//  for (int it = 0; it < max_iterations_; it++) {
//    std::pair<int, Eigen::MatrixXd> params1 = OptimizeStep(params);
//    if (params1.first == 0 || params1.first == 1) {
//      //params = params1.second;
//      params = pos_->GetPoseVector();
//      Eigen::VectorXd pose_i = pos_->GetPoseVector();
//      std::cout << "Step " << it+1 << " / " << max_iterations_ << " : "
//                << params1.first << " : "
//                //<< params1.second.transpose() << " : "
//                << ComputeResidual(params) << std::endl
//                << "lambda = " << lambda_ << std::endl;
//
//      std::cout << "  " << ComputeResidual(params1.second) << " " << params1.second.transpose() << std::endl;
//      std::cout << "  " << ComputeResidual(params) << " " << params.transpose() << std::endl;
//      std::cout << "  " << ComputeResidual(pos_->GetPoseVector()) << " " << pos_->GetPoseVector().transpose() << std::endl;
//      std::cout << std::endl;
//
//      if (params1.first == 0) {
//        break;
//      } 
//    } 
//    //else {
//    //  std::cout << "Step " << it+1 << " / " << max_iterations_ << " : "
//    //            << params1.first << " : "
//    //            //<< params1.second.transpose() << " : "
//    //            << ComputeResidual(params1.second) << std::endl;
//    //}
//  }
//
//  return std::make_pair(1, params);
//}




std::pair<int, Eigen::VectorXd> LevenbergMarquardt::Optimize(const Eigen::VectorXd params0) {

  Eigen::VectorXd params = params0;

  double error = std::numeric_limits<double>::max();
  int it = 0;

  while (error > tolerance_ && it < max_iterations_) {
    it++;

    std::pair<int, double> params1 = OptimizeStep(params);
    error = params1.second;

    //if ( it == 3)
    //  break;

    if (params1.first == 1) {
      params = pos_->GetPoseVector();

      std::cout << "Step " << it+1 << " / " << max_iterations_ << " : "
                << params1.first << " : "
                << ComputeResidual(params) << std::endl
                << "lambda = " << lambda_ << std::endl;
    } 
  }

  std::cout << "LM terminated at iteration " << it << " with error " << error << std::endl;

  return std::make_pair(1, params);
}


std::pair<int, double> LevenbergMarquardt::OptimizeStep(const Eigen::VectorXd& params0) {

  //std::cout << " params0 = " << params0.transpose() << std::endl;

  double residual0 = ComputeResidual(params0);

  Eigen::MatrixXd jacobian = ComputeJacobian(params0, residual0, 0.0001, 0.001, 0.0001);

  // Compute the approximate Hessian
  //Eigen::MatrixXd hessian = jacobian.transpose() * jacobian;
  Eigen::MatrixXd hessian = jacobian.transpose() * jacobian;
  //std::cout << "hessian = " << std::endl << hessian << std::endl << std::endl;

  
  //hessian += lambda_ * Eigen::MatrixXd::Identity(hessian.rows(), hessian.cols());  
  //std::cout << "hessian = " << std::endl << hessian << std::endl << std::endl;
  
  // Update the diagonal of the Hessian with the damping factor
  //for (int i = 0; i < hessian.rows(); i++) {
  //  hessian(i,i) += lambda_ * hessian(i,i);
  //}
  hessian += lambda_ * Eigen::MatrixXd::Identity(hessian.rows(), hessian.cols());


  // Compute the gradient
  Eigen::VectorXd gradient = jacobian.transpose() * residual0;

  // Solve for the update (delta)
  //Eigen::VectorXd delta = hessian.llt().solve(-gradient);
  Eigen::VectorXd delta = hessian.ldlt().solve(-gradient);

  std::cout << " parm0 = " << params0.transpose() << std::endl;
  std::cout << " delta = " << delta.transpose() << std::endl;

  // Update the parameters
  Eigen::VectorXd params = params0 + delta;

  // Compute the new residual
  double residual1 = ComputeResidual(params);
  
  // Compute delta norm
  double delta_norm = delta.norm();

  
  // If the new residual is smaller, accept the update, else reject it
  if (residual1 < residual0) {
    std::cout << " residual = " << residual0 << " -> " << residual1 << std::endl;

    lambda_ /= damping_factor_;
    residual_ = residual1;

    pos_->TransformToPose(std::make_pair(params.segment(3,4), params.segment(0,3)), params(7));
    
    return std::make_pair(1, residual1);
  }
  else {
    lambda_ *= damping_factor_;
    return std::make_pair(0, residual0);
  }
  
}

double LevenbergMarquardt::GetResidual() {
  return residual_;
}

Eigen::MatrixXd  LevenbergMarquardt::ComputeJacobian(const Eigen::VectorXd& params, double residual_original, 
                                                     double delta_t, double delta_q, double delta_s) {
  Eigen::MatrixXd jacobian(1, params.size());

  //Eigen::VectorXd deltas = params;
  //deltas.segment(0,3) *= delta_t;
  //deltas(0) = deltas(0) * delta_t;
  //deltas(1) = deltas(1) * delta_t;
  //deltas(2) = deltas(2) * delta_t;
  //deltas(3) = delta_q;
  //deltas(4) = delta_q;
  //deltas(5) = delta_q;
  //deltas(6) = delta_q;
  //deltas(7) = delta_s;

  // vector deltas equal to 0.0001
  Eigen::VectorXd deltas = Eigen::VectorXd::Constant(params.size(), 0.0001);

  for (int i = 0; i < params.size(); i++) {
    Eigen::VectorXd params_delta = params;
    params_delta(i) += deltas(i);
    double residual_delta = ComputeResidual(params_delta);
    jacobian(0, i) = (residual_delta - residual_original) / deltas(i);
  }

  //std::cout << "jacobian = " << std::endl << jacobian << std::endl << std::endl;
  
  return jacobian;
}