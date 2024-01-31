#include "fea.hpp"

FEA::FEA(std::string element_type,
         float young_modulus, float poisson_coefficient, 
         bool debug_mode) : E_(young_modulus), nu_(poisson_coefficient), debug_mode_(debug_mode) {

  if (element_type == "C3D6") {
    //setElement(new C3D6(young_modulus, poisson_coefficient));
    element_ = new C3D6(young_modulus, poisson_coefficient);
  } else if (element_type == "C3D8") {
    //setElement(new C3D8(young_modulus, poisson_coefficient));
    element_ = new C3D8(young_modulus, poisson_coefficient);
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
    std::cout << "element D: " << element_->getD() << std::endl;
  }
}


void FEA::MatAssembly(std::vector<Eigen::Vector3d> &vpts, 
                      std::vector<std::vector<unsigned int>> &velts) {
  K_ = element_->matAssembly(vpts, velts);
}


void FEA::ApplyBoundaryConditions(BoundaryConditions &bc) {
  // Reset F_ to zero
  F_ = Eigen::MatrixXd::Zero(K_.rows(), 1);

  int bc_encastre = 0;
  int bc_force = 0;
  int bc_displ = 0;

  bc_.resize(bc.NodeIds().size());

  for (unsigned int node = 0; node < bc.NodeIds().size(); node++) {
    bc_[node] = bc.NodeIds()[node];

    if (!bc.NodeIds()[node])
      continue;
    
    std::vector<unsigned int>* dof = bc.Dof(node);
    std::vector<double> values = bc.Values(node);

    for (unsigned int i=0; i<dof->size(); i++) {
      unsigned int m = 3*node + i;
      if ((*dof)[i] == 0) {
        // set row and col m to zero
        for (unsigned int j=0; j<K_.cols(); j++) {
          K_(m,j) = 0.0;
          K_(j,m) = 0.0;
        }
        K_(m,m) = 1.0;

        if (values[i] != 0.0) {
          F_(m,0) = values[i];
          bc_displ++;
        } else {
          bc_encastre++;
        }
      }
      else {
        F_(m,0) = values[i];
        bc_force++;
      }  
    }
  }

  std::cout << " - Encastre conditions: " << bc_encastre << std::endl;
  std::cout << " - Force conditions: " << bc_force << std::endl;
  std::cout << " - Displacement conditions: " << bc_displ << std::endl;
}



void FEA::Solve() {
  // Check if K is singular or ill-conditioned
  if (IsSingularOrIllConditioned2(K_)) {
    std::cout << "K is singular or ill-conditioned" << std::endl;
    return;
  }
  std::cout << "K is OK: not singular or ill-conditioned" << std::endl;

  // Rearrange the matrix for efficiency
  U_ = SolveSystemWithPreconditioning(K_, F_, "CG");

//  // Solve using the appropiate method for efficiency
//  if (UseDirectSolver(K)) {
//    U_ = DirectSolver(K, F_);
//  } else {
//    U_ = IterativeSolver(K, F_);
//  }
}

bool FEA::IsSingularOrIllConditioned1(Eigen::MatrixXd &A) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
  double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
  if (cond > 1e10) {
    std::cout << "Condition number: " << cond << std::endl;
    return true;
  }
  return false;
}

bool FEA::IsSingularOrIllConditioned2(Eigen::MatrixXd &A) {
    // Threshold for determining if the determinant is effectively zero
    const double detThreshold = 1e-10;
    // Threshold for numerical instability in solutions, indicative of ill-conditioning
    const double solutionThreshold = 1e-10;

    // Check if the matrix is singular by determinant
    double det = A.determinant();
    if (std::fabs(det) < detThreshold) {
        return true; // Matrix is singular or near-singular
    }

    // Attempt to solve A*x = b for a unit vector b to check for ill-conditioning
    Eigen::VectorXd b = Eigen::VectorXd::Ones(A.rows());
    Eigen::VectorXd x = A.fullPivLu().solve(b);

    // Check if the solution is significantly different from what is expected
    Eigen::VectorXd b_check = A * x; // Recompute b to check solution
    if ((b - b_check).norm() / b.norm() > solutionThreshold) {
        return true; // Solution instability suggests ill-conditioning
    }

    // If neither singularity nor significant solution instability was detected
    return false;
}



Eigen::VectorXd FEA::SolveSystemWithPreconditioning(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b, std::string method) {

  Eigen::VectorXd b_vec = MatrixToVector(b);

  if (method == "CG") {
    return SolveSystemWithPreconditioningCG(A, b_vec);
  } else if (method == "BiCGSTAB") {
    return SolveSystemWithPreconditioningBiCGSTAB(A, b_vec);
  } else if (method == "Lu") {
    return SolveSystemWithLU(A, b_vec);
  } else if (method == "LuFull") {
    return SolveSystemWithFullPivLU(A, b_vec);
  } else {
    std::cout << "Method not supported" << std::endl;
    return Eigen::VectorXd::Zero(b.size());
  }

}

// Function to solve a system of linear equations A*x = b using an iterative solver with diagonal preconditioning
Eigen::VectorXd FEA::SolveSystemWithPreconditioningCG(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Convert the dense matrix A to a sparse matrix for efficiency with iterative solvers
    Eigen::SparseMatrix<double> Asparse = A.sparseView();

    // Initialize the iterative solver with diagonal preconditioning
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    cg.compute(Asparse);

    // Check if the factorization was successful
    if(cg.info() != Eigen::Success) {
        std::cerr << "Matrix decomposition failed. Returning zero vector." << std::endl;
        return Eigen::VectorXd::Zero(b.size());
    }

    // Solve the system A*x = b
    Eigen::VectorXd x = cg.solve(b);

    // Return the solution vector
    return x;
}

// Function to solve a system of linear equations A*x = b using Eigen::BiCGSTAB with diagonal preconditioning
Eigen::VectorXd FEA::SolveSystemWithPreconditioningBiCGSTAB(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Convert the dense matrix A to a sparse matrix format for efficiency
    Eigen::SparseMatrix<double> Asparse = A.sparseView();

    // Initialize the BiCGSTAB solver with diagonal preconditioning for the sparse matrix
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> bicgstab;
    bicgstab.compute(Asparse);

    // Check if the decomposition was successful
    if (bicgstab.info() != Eigen::Success) {
        std::cerr << "Matrix decomposition failed. Returning zero vector." << std::endl;
        return Eigen::VectorXd::Zero(b.size());
    }

    // Solve the system A*x = b
    Eigen::VectorXd x = bicgstab.solve(b);

    return x;
}

Eigen::VectorXd FEA::SolveSystemWithLU(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Compute the LU decomposition of A
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);

    // Use the LU decomposition to solve for x
    Eigen::VectorXd x = lu.solve(b);

    return x;
}


Eigen::VectorXd FEA::SolveSystemWithFullPivLU(const Eigen::MatrixXd &A, const Eigen::VectorXd &b) {
    // Compute the LU decomposition of A with full pivoting
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);

    // Before solving, check if the solution is possible and the system is not singular
    if (!lu.isInvertible()) {
        std::cerr << "Matrix is singular or nearly singular, solution may not be accurate or possible." << std::endl;
        return Eigen::VectorXd::Zero(b.size()); // Return a zero vector if not invertible
    }

    // Use the LU decomposition to solve for x
    Eigen::VectorXd x = lu.solve(b);

    return x;
}

// Function to convert a column matrix (Eigen::MatrixXd) to Eigen::VectorXd
Eigen::VectorXd FEA::MatrixToVector(const Eigen::MatrixXd &columnMatrix) {
    // Ensure the matrix is indeed a column matrix
    if (columnMatrix.cols() != 1) {
        std::cerr << "Error: The input matrix is not a column matrix." << std::endl;
        return Eigen::VectorXd(); // Return an empty vector in case of error
    }

    // Directly convert the column matrix to a vector
    Eigen::VectorXd vector = columnMatrix.col(0);

    return vector;
}


//Eigen::MatrixXd PreconditionMatrix(Eigen::MatrixXd &A) {
//    // Convert the input dense matrix A to a sparse matrix for efficient preconditioning
//    Eigen::SparseMatrix<double> Asparse = A.sparseView();
//
//    // Initialize the diagonal preconditioner
//    Eigen::DiagonalPreconditioner<double> preconditioner;
//
//    // Compute the preconditioner based on the sparse matrix
//    preconditioner.compute(Asparse);
//
//    // Although directly applying the preconditioner to the matrix is not typical,
//    // for demonstration, we show how the preconditioned matrix might be formed.
//    // In practical applications, you would use the preconditioner directly with an iterative solver.
//    Eigen::MatrixXd Apreconditioned(Asparse.rows(), Asparse.cols());
//    
//    // The following loop applies the diagonal preconditioning effect manually for demonstration.
//    // This mimics the preconditioner's effect by scaling each row by the inverse of its diagonal element.
//    for (int i = 0; i < Asparse.rows(); ++i) {
//        double diagVal = preconditioner.diag()(i);
//        if (diagVal != 0) {  // Avoid division by zero
//            for (int j = 0; j < Asparse.cols(); ++j) {
//                Apreconditioned(i, j) = A(i, j) / diagVal;
//            }
//        }
//    }
//
//    return Apreconditioned;
//}


//bool FEA::UseDirectSolver(Eigen::MatrixXd &A);
//Eigen::MatrixXd FEA::DirectSolver(Eigen::MatrixXd &A, Eigen::MatrixXd &b);
//Eigen::MatrixXd FEA::IterativeSolver(Eigen::MatrixXd &A, Eigen::MatrixXd &b);

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

void FEA::ReportNodes(std::string filename) {
  // Transform U_ to size nx3
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(U_.rows()/3, 3);
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(F_.rows()/3, 3);
  for (unsigned int n=0; n<U.rows(); n++) {
    U(n/3, n%3) = U_(n,0);
    F(n/3, n%3) = F_(n,0);
  }

  std::ofstream file;
  file.open(filename);
  file << "n\tu.x\tu.y\tu.z\tf.x\tf.y\tf.z" << std::endl;

  for (unsigned int n=0; n<U.rows(); n++) {
      file << n << "\t" << U(n,0) << "\t" << U(n,1) << "\t" << U(n,2) << "\t" << F(n,0) << "\t" << F(n,1) << "\t" << F(n,2) << std::endl;
  }

  file.close();
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
}

void FEA::ExportK(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << K_ << std::endl;
  file << std::endl << "K ( " << K_.rows() << ", " << K_.cols() << " ) " << std::endl;
  file.close();
}

void FEA::ExportF(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << F_ << std::endl;
  file << std::endl << "F ( " << F_.rows() << ", " << F_.cols() << " ) " << std::endl;
  file.close();
}

void FEA::ExportU(std::string filename) {
  std::ofstream file;
  file.open(filename);
  file << U_ << std::endl;
  file << std::endl << "U ( " << U_.rows() << ", " << U_.cols() << " ) " << std::endl;
  file.close();
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