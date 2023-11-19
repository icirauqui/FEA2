#include "fea.hpp"


FEA::FEA(int frame_id, 
         std::string element_type, 
         float young_modulus, float poisson_coefficient, float element_depth, 
         float gauss_point,
         bool debug_mode) {
  frame_id_ = frame_id;
  element_ = element_type;
  E_ = young_modulus;
  nu_ = poisson_coefficient;
  h_ = element_depth;
  debug_mode_ = debug_mode;

  lambda_ = ( nu_*E_ ) / ( (1+nu_) * (1-2*nu_) );
  G_ = E_ / ( 2 * (1+nu_) );

  if (element_ == "C3D6") {
    InitC3D6();
  } else if (element_ == "C3D8") {
    InitC3D8();
  } else {
    std::cout << "Element not supported" << std::endl;
  }

  InitGaussPoints(gauss_point);

  if (debug_mode_){
    std::cout << "D_ = " << std::endl << D_ << std::endl;
    std::cout << "gs_ = " << std::endl << gs_ << std::endl;
  }
}



void FEA::MatAssembly(std::vector<std::vector<float> > &vpts, 
                      std::vector<std::vector<int> > &velts) {

  if (element_ == "C3D6") {
    K_ = c3d6::matAssembly(vpts, velts, E_, nu_);
  } else if (element_ == "C3D8") {
    K_ = c3d8::matAssembly(vpts, velts, E_, nu_);
  } else {
    std::cout << "Element not supported" << std::endl;
  }

  Print_K();
  Eigenvalues_K();
}

void FEA::Eigenvalues_K() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_);
  std::cout << std::endl
            << "The eigenvalues of K_ are:" << std::endl
            << es.eigenvalues().transpose() << std::endl;

}

void FEA::Print_K() {
  std::cout << std::endl
            << "K_ = " << std::endl
            << K_ << std::endl;
}


void FEA::ComputeForces() {
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
  F_ = K_ * U_;
}

void FEA::SetForces(std::vector<std::vector<float>> &vF) {
  F_ = Eigen::MatrixXd::Zero(3*vF.size(), 1);
  for (unsigned int i = 0; i < vF.size(); i++) {
    F_(i*3, 0) = vF[i][0];
    F_(i*3+1, 0) = vF[i][1];
    F_(i*3+2, 0) = vF[i][2];
  }
}



void FEA::EncastreBackLayer(float k_large) {
  std::vector<int> dir;
  for (unsigned int i = K_.rows()/2; i < K_.rows(); i++) {
    // Set row and column to zero
    for (unsigned int j = 0; j < K_.cols(); j++) {
      K_(i,j) = 0.0;
      K_(j,i) = 0.0;
    }
    // Set diagonal to k_large
    K_(i,i) = k_large;
  }
}

void FEA::ImposeDirichletEncastre(std::vector<int> &dir, float k_large) {
  for (auto d : dir) {
    int mp0 = 3*(d - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    K_(mp0,mp0) = k_large;
    K_(mp1,mp1) = k_large;
    K_(mp2,mp2) = k_large;

    F_(mp0,0) = 0.0;
    F_(mp1,0) = 0.0;
    F_(mp2,0) = 0.0;
  }
}


void FEA::ImposeDirichletEncastre(std::vector<std::vector<int>> &dir, float k_large) {
  for (auto d : dir) {
    int mp0 = 3*(d[0] - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    K_(mp0,mp0) = k_large;
    K_(mp1,mp1) = k_large;
    K_(mp2,mp2) = k_large;

    F_(mp0,0) = 0.0;
    F_(mp1,0) = 0.0;
    F_(mp2,0) = 0.0;
  }

}


void FEA::ComputeDisplacements() {
  K1_ = K_.inverse();
  U_ = K1_ * F_;
}


void FEA::ComputeStrainEnergy() {
  sE_ = (U_.transpose() * F_)(0,0);
  sE_ = std::abs(sE_);
}

double FEA::ComputeStrainEnergy(std::vector<Eigen::Vector3d> &u0,
                                std::vector<Eigen::Vector3d> &u1) {
  int dim_in = u0.size() * 3;
  int dim_k = K_.rows();

  if (dim_in != dim_k) {
    std::cout << "Error: dim_in != dim_k" << std::endl;
    return -1.0;
  }

  //std::cout << std::endl;
  U_ = Eigen::MatrixXd::Zero(dim_in, 1);
  for (unsigned int i = 0; i < u0.size(); i++) {
    U_(i*3, 0) = u1[i][0] - u0[i][0];
    U_(i*3+1, 0) = u1[i][1] - u0[i][1];
    U_(i*3+2, 0) = u1[i][2] - u0[i][2];
  }

  //for (unsigned int i=0; i<u0.size(); i++) {
  //  std::cout << i << ": " << u0[i].transpose() << " -> " << u1[i].transpose() << " = ";
  //  std::cout << U_(i*3, 0) << " " << U_(i*3+1, 0) << " " << U_(i*3+2, 0) << std::endl;
  //}



  ComputeForces();

  //std::cout << std::endl;
  //std::cout << "K_ = " << std::endl;
  //for (unsigned int i=0; i<K_.rows(); i++) {
  //  for (unsigned int j=0; j<K_.cols(); j++) {
  //    std::cout << K_(i,j) << " ";
  //  }
  //  std::cout << std::endl;
  //}
  
  //std::cout << std::endl;
  //std::cout << "F_ = " << std::endl;
  //for (unsigned int i=0; i<u0.size(); i++) {
  //  std::cout << i << ": " << F_(i*3, 0) << " " << F_(i*3+1, 0) << " " << F_(i*3+2, 0) << std::endl;
  //}

  ComputeStrainEnergy();

  return sE_;
}

Eigen::MatrixXd FEA::K() {
  return K_;
}


Eigen::MatrixXd FEA::F() {
  return F_;
}


Eigen::MatrixXd FEA::U() {
  return U_;
}


float FEA::StrainEnergy() {
  return sE_;
}


void FEA::InitC3D6() {
  D_ <<  lambda_+2*G_ ,   lambda_    ,   lambda_    ,     0.0     ,     0.0     ,     0.0     ,
            lambda_   , lambda_+2*G_ ,   lambda_    ,     0.0     ,     0.0     ,     0.0     ,
            lambda_   ,   lambda_    , lambda_+2*G_ ,     0.0     ,     0.0     ,     0.0     ,
              0.0     ,     0.0      ,     0.0      ,      G_     ,     0.0     ,     0.0     ,
              0.0     ,     0.0      ,     0.0      ,     0.0     ,      G_     ,     0.0     ,
              0.0     ,     0.0      ,     0.0      ,     0.0     ,     0.0     ,      G_     ;
  Kei_ = Eigen::MatrixXd::Zero(18, 18);
  dndgs_ = Eigen::MatrixXd::Zero(6, 3);
  B_ = Eigen::MatrixXd::Zero(6, 18);
  base_size_ = 6;
}


void FEA::InitC3D8() {
  D_ <<  lambda_+2*G_ ,   lambda_    ,   lambda_    ,     0.0     ,     0.0     ,     0.0     ,
            lambda_   , lambda_+2*G_ ,   lambda_    ,     0.0     ,     0.0     ,     0.0     ,
            lambda_   ,   lambda_    , lambda_+2*G_ ,     0.0     ,     0.0     ,     0.0     ,
              0.0     ,     0.0      ,     0.0      ,      G_     ,     0.0     ,     0.0     ,
              0.0     ,     0.0      ,     0.0      ,     0.0     ,      G_     ,     0.0     ,
              0.0     ,     0.0      ,     0.0      ,     0.0     ,     0.0     ,      G_     ;
  Kei_ = Eigen::MatrixXd::Zero(24, 24);
  dndgs_ = Eigen::MatrixXd::Zero(8, 3);
  B_ = Eigen::MatrixXd::Zero(6, 24);
  base_size_ = 8;
}


void FEA::InitGaussPoints(float fg) {
  gs_ << -fg, -fg, -fg,
         +fg, -fg, -fg,
         +fg, +fg, -fg,
         -fg, +fg, -fg,
         -fg, -fg, +fg,
         +fg, -fg, +fg,
         +fg, +fg, +fg,
         -fg, +fg, +fg;
}





//void FEA::ComputeKei(std::vector<std::vector<float>> &vfPts) {
//  // vFpts to Eigen::MatrixXd
//  Eigen::MatrixXd vFpts(base_size_,3);
//  for (unsigned int i=0; i<base_size_; i++) {
//    vFpts(i,0) = vfPts[i][0];
//    vFpts(i,1) = vfPts[i][1];
//    vFpts(i,2) = vfPts[i][2];
//  }
//
//  std::cout << "Element node coordinates:" << std::endl
//            << vFpts << std::endl;
//
//  for (unsigned int ops=0; ops<gs_.rows(); ops++) {
//    float xi   = gs_(ops, 0);
//    float eta  = gs_(ops, 1);
//    float zeta = gs_(ops, 2);
//
//    dNdgs(xi, eta, zeta, base_size_);
//
//    // Compute Jacobian
//    J_ = dndgs_.transpose() * vFpts;
//
//    // Compute Jacobian determinant
//    float Jdet = J_.determinant();
//
//    // Compute inverse of Jacobian
//    J1_ = J_.inverse();
//
//    // Compute dNdxyz
//    Eigen::MatrixXd dNdxyz = (J1_ * dndgs_.transpose()).transpose();
//
//    // Compute B matrix
//    for (unsigned int i=0; i<base_size_; i++) {
//      B_(0,3*i+0) = dNdxyz(i,0); B_(0,3*i+1) = 0;           B_(0,3*i+2) = 0;
//      B_(1,3*i+0) = 0;           B_(1,3*i+1) = dNdxyz(i,1); B_(1,3*i+2) = 0;
//      B_(2,3*i+0) = 0;           B_(2,3*i+1) = 0;           B_(2,3*i+2) = dNdxyz(i,2);
//      B_(3,3*i+0) = dNdxyz(i,1); B_(3,3*i+1) = dNdxyz(i,0); B_(3,3*i+2) = 0;
//      B_(4,3*i+0) = 0;           B_(4,3*i+1) = dNdxyz(i,2); B_(4,3*i+2) = dNdxyz(i,1);
//      B_(5,3*i+0) = dNdxyz(i,2); B_(5,3*i+1) = 0;           B_(5,3*i+2) = dNdxyz(i,0);
//    }
//
//    // Compute Bt * D * B * Jdet and accumulate to vBtDB
//    if (ops == 0) {
//      Kei_ = B_.transpose() * D_ * B_ * Jdet;
//    } else {
//      Kei_ += B_.transpose() * D_ * B_ * Jdet;
//    }
//  }
//
//  std::cout << "Kei = " << std::endl
//            << Kei_ << std::endl;
//}


void FEA::dNdgs(float xi, float eta, float zeta, int dim) {
  // Col 0 = dN/dgs(0,0), Col 1 = dN/dgs(0,1), Col 2 = dN/dgs(0,2)

  float fact = 1.0 / dim;

  if (element_ == "C3D6") {
    dndgs_(0,0) = -(1 + zeta)/2;    dndgs_(0,1) = -(1 + zeta)/2;    dndgs_(0,2) =  (1-xi-eta)/2;
    dndgs_(1,0) =  (1 + zeta)/2;    dndgs_(1,1) =   0.0;            dndgs_(1,2) =  xi/2;
    dndgs_(2,0) =  0.0;             dndgs_(2,1) =  (1 + zeta)/2;    dndgs_(2,2) =  eta/2;
    dndgs_(3,0) = -(1 - zeta)/2;    dndgs_(3,1) = -(1 - zeta)/2;    dndgs_(3,2) = -(1-xi-eta)/2;
    dndgs_(4,0) =  (1 - zeta)/2;    dndgs_(4,1) =   0.0;            dndgs_(4,2) = -xi/2;
    dndgs_(5,0) =  0.0;             dndgs_(5,1) =  (1 - zeta)/2;    dndgs_(5,2) = -eta/2;
  } else if (element_ == "C3D8") {
    dndgs_(0,0) = -0.125 * ( (1-eta) * (1-zeta) );     dndgs_(0,1) = -0.125 * ( (1-xi) * (1-zeta) );     dndgs_(0,2) = -0.125 * ( (1-xi) * (1-eta) ); 
    dndgs_(1,0) = +0.125 * ( (1-eta) * (1-zeta) );     dndgs_(1,1) = -0.125 * ( (1+xi) * (1-zeta) );     dndgs_(1,2) = -0.125 * ( (1+xi) * (1-eta) ); 
    dndgs_(2,0) = +0.125 * ( (1+eta) * (1-zeta) );     dndgs_(2,1) = +0.125 * ( (1+xi) * (1-zeta) );     dndgs_(2,2) = -0.125 * ( (1+xi) * (1+eta) ); 
    dndgs_(3,0) = -0.125 * ( (1+eta) * (1-zeta) );     dndgs_(3,1) = +0.125 * ( (1-xi) * (1-zeta) );     dndgs_(3,2) = -0.125 * ( (1-xi) * (1+eta) ); 
    dndgs_(4,0) = -0.125 * ( (1-eta) * (1+zeta) );     dndgs_(4,1) = -0.125 * ( (1-xi) * (1+zeta) );     dndgs_(4,2) = +0.125 * ( (1-xi) * (1-eta) ); 
    dndgs_(5,0) = +0.125 * ( (1-eta) * (1+zeta) );     dndgs_(5,1) = -0.125 * ( (1+xi) * (1+zeta) );     dndgs_(5,2) = +0.125 * ( (1+xi) * (1-eta) ); 
    dndgs_(6,0) = +0.125 * ( (1+eta) * (1+zeta) );     dndgs_(6,1) = +0.125 * ( (1+xi) * (1+zeta) );     dndgs_(6,2) = +0.125 * ( (1+xi) * (1+eta) ); 
    dndgs_(7,0) = -0.125 * ( (1+eta) * (1+zeta) );     dndgs_(7,1) = +0.125 * ( (1-xi) * (1+zeta) );     dndgs_(7,2) = +0.125 * ( (1-xi) * (1+eta) ); 
  }
}