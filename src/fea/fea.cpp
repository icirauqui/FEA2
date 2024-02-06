#include "fea.hpp"

FEA::FEA(std::string element_type,
         float young_modulus, float poisson_coefficient, 
         bool debug_mode) : E_(young_modulus), nu_(poisson_coefficient), debug_mode_(debug_mode) {

  if (element_type == "C2D4") {
    element_ = new C2D4(young_modulus, poisson_coefficient);
    //bcs_ = new BoundaryConditions2d(2, nullptr);
    //loads_ = new Loads2d(2, nullptr);
  } else if (element_type == "C3D6") {
    //setElement(new C3D6(young_modulus, poisson_coefficient));
    element_ = new C3D6(young_modulus, poisson_coefficient);
    //bcs_ = new BoundaryConditions3d(3, nullptr);
    //loads_ = new Loads3d(3, nullptr);
  } else if (element_type == "C3D8") {
    //setElement(new C3D8(young_modulus, poisson_coefficient));
    element_ = new C3D8(young_modulus, poisson_coefficient);
    //bcs_ = new BoundaryConditions3d(3, nullptr);
  } else {
    std::cout << "Element not supported" << std::endl;
  }

  if (debug_mode_) {
    std::cout << std::endl;
    std::cout << "FEA: " << element_type << std::endl;
    std::cout << "element pointer address: " << element_ << std::endl;
    std::cout << "element name: " << element_->getElementName() << std::endl;
    std::cout << "element dim: " << element_->getNumNodes() << std::endl;
    std::cout << "element dof: " << element_->getDofPerNode() << std::endl;
    std::cout << "element D: " << std::endl << element_->getD() << std::endl;
  }
}


void FEA::MatAssembly(std::vector<Eigen::Vector2d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {

  std::vector<Eigen::Vector3d> vpts3d;
  for (auto v : vpts) {
    vpts3d.push_back(Eigen::Vector3d(v[0], v[1], 0.0));
  }
  
  K_ = element_->matAssembly(vpts3d, velts);

  // Reset F_ and U_ to zero
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
  U_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
}


void FEA::MatAssembly(std::vector<Eigen::Vector3d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {
  K_ = element_->matAssembly(vpts, velts);

  // Reset F_ and U_ to zero
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
  U_ = Eigen::MatrixXd::Zero(K_.rows(), 1);
}




void FEA::ApplyBoundaryConditions(BoundaryConditions &bc) {

  int bc_encastre = 0;
  int bc_displ = 0;

  int num_dof = bc.NumDof();

  std::cout << " K Size = " << K_.rows() << " x " << K_.cols() << std::endl;
  std::cout << " F Size = " << F_.rows() << " x " << F_.cols() << std::endl;

  bc_.resize(bc.NodeIds().size());

  for (unsigned int node = 0; node < bc.NodeIds().size(); node++) {
    bc_[node] = bc.NodeIds()[node];

    if (!bc.NodeIds()[node])
      continue;
    
    std::vector<double> values = bc.Values(node);

    for (unsigned int i=0; i<values.size(); i++) {
      unsigned int m = num_dof*node + i;

      if (values[i] == 0) {
        for (unsigned int j=0; j<K_.cols(); j++) {
          K_(m,j) = 0.0;
          K_(j,m) = 0.0;
        }
        K_(m,m) = 1.0;
        F_(m,0) = 0.0;
        bc_encastre++;
      }
    }
  }

  std::cout << " - Encastre conditions: " << bc_encastre << std::endl;
  std::cout << " - Displacement conditions: " << bc_displ << std::endl;
}





void FEA::ApplyLoads(Loads &loads) {

  int loads_node = 0;

  int num_dof = loads.NumDof();

  loads_.resize(loads.NodeIds().size());

  for (unsigned int node = 0; node < loads.NodeIds().size(); node++) {
    loads_[node] = loads.NodeIds()[node];

    if (!loads.NodeIds()[node])
      continue;
    
    std::vector<double> values = loads.Values(node);

    for (unsigned int i=0; i<values.size(); i++) {
      unsigned int m = num_dof*node + i;
      F_(m,0) = values[i];
      loads_node++;
    }
  }

  std::cout << " - Node loads: " << loads_node << std::endl;
}



void FEA::Solve(std::string method) {
  // Check if K is singular or ill-conditioned
  if (solver::IsSingularOrIllConditioned2(K_) == 1) {
    return;
  }

  std::cout << " - Solving system with " << method << " method" << std::endl;
  U_ = solver::SolveSystemWithPreconditioning(K_, F_, method);

  Fi_ = K_ * U_;
}








// 
// LEGACY FEA
// 


void FEA::EncastreBackLayer() {
  std::vector<int> dir;
  for (unsigned int i = K_.rows()/2; i < K_.rows(); i++) {
    // Set row and column to zero
    for (unsigned int j = 0; j < K_.cols(); j++) {
      K_(i,j) = 0.0;
      K_(j,i) = 0.0;
    }
    // Set diagonal to k_large_
    K_(i,i) = k_large_;
  }
}


void FEA::ImposeDirichletEncastre(std::vector<int> &dir) {
  for (auto d : dir) {
    int mp0 = 3*(d - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    K_(mp0,mp0) = k_large_;
    K_(mp1,mp1) = k_large_;
    K_(mp2,mp2) = k_large_;

    F_(mp0,0) = 0.0;
    F_(mp1,0) = 0.0;
    F_(mp2,0) = 0.0;
  }
}


void FEA::ImposeDirichletEncastre(std::vector<std::vector<int>> &dir) {
  for (auto d : dir) {
    int mp0 = 3*(d[0] - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    K_(mp0,mp0) = k_large_;
    K_(mp1,mp1) = k_large_;
    K_(mp2,mp2) = k_large_;

    F_(mp0,0) = 0.0;
    F_(mp1,0) = 0.0;
    F_(mp2,0) = 0.0;
  }

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

  ComputeForces();

  ComputeStrainEnergy();

  return sE_;
}





// 
// REPORT FUNCTIONS
// 

std::string buildHeader() {
  int len_field = 12;
  std::vector<std::string> hdrs = {"n", "f", "f.x", "f.y", "f.z", "fi", "fi.x", "fi.y", "fi.z", "u", "u.x", "u.y", "u.z"};

  std::string header;
  for (auto h : hdrs) {
    for (unsigned int i = 0; i < len_field - h.size(); i++) {
      header += " ";
    }
    header += h + " ";
  }

  return header;
}


void FEA::ReportNodes(std::string filename) {

  // Eigen Matrix of size U_.rows()/3 x 7
  // n, u.x, u.y, u.z, f.x, f.y, f.z
  Eigen::MatrixXd exp = Eigen::MatrixXd::Zero(U_.rows()/3, 13);

  for (unsigned int n=0; n<exp.rows(); n++) {
    exp(n,0) = n;

    exp(n,2) = F_(n*3, 0);
    exp(n,3) = F_(n*3+1, 0);
    exp(n,4) = F_(n*3+2, 0);
    exp(n,1) = sqrt(exp(n,2)*exp(n,2) + exp(n,3)*exp(n,3) + exp(n,4)*exp(n,4));

    exp(n,6) = Fi_(n*3, 0);
    exp(n,7) = Fi_(n*3+1, 0);
    exp(n,8) = Fi_(n*3+2, 0);
    exp(n,5) = sqrt(exp(n,6)*exp(n,6) + exp(n,7)*exp(n,7) + exp(n,8)*exp(n,8));

    exp(n,10) = U_(n*3, 0);
    exp(n,11) = U_(n*3+1, 0);
    exp(n,12) = U_(n*3+2, 0);
    exp(n,9) = sqrt(exp(n,10)*exp(n,10) + exp(n,11)*exp(n,11) + exp(n,12)*exp(n,12));
  }

  std::ofstream file;
  file.open(filename);
  file << buildHeader() << std::endl;
  file << exp << std::endl;

  file.close();

  std::cout << "ReportNodes [U,F]: " << filename << std::endl;
}

void FEA::ExportAll(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << K_ << std::endl;
  file << std::endl << "K ( " << K_.rows() << ", " << K_.cols() << " ) " << std::endl << std::endl;
  file << F_ << std::endl;
  file << std::endl << "F ( " << F_.rows() << ", " << F_.cols() << " ) " << std::endl << std::endl;
  file << U_ << std::endl;
  file << std::endl << "U ( " << U_.rows() << ", " << U_.cols() << " ) " << std::endl << std::endl;
  file.close();

  std::cout << "ExportAll [K,F,U]: " << filename << std::endl;
}

void FEA::ExportK(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << K_ << std::endl;
  file << std::endl << "K ( " << K_.rows() << ", " << K_.cols() << " ) " << std::endl;
  file.close();

  std::cout << "ExportK [K]: " << filename << std::endl;
}

void FEA::ExportF(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << F_ << std::endl;
  file << std::endl << "F ( " << F_.rows() << ", " << F_.cols() << " ) " << std::endl;
  file.close();

  std::cout << "ExportF [F]: " << filename << std::endl;
}

void FEA::ExportU(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << U_ << std::endl;
  file << std::endl << "U ( " << U_.rows() << ", " << U_.cols() << " ) " << std::endl;
  file.close();

  std::cout << "ExportU [U]: " << filename << std::endl;
}


void FEA::PrintEigenvaluesK() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_);
  std::cout << std::endl
            << "The eigenvalues of K_ are:" << std::endl
            << es.eigenvalues().transpose() << std::endl;

}

void FEA::PrintK() {
  std::cout << std::endl
            << "K_ = " << std::endl
            << K_ << std::endl;
}