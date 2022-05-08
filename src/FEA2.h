#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>
#include <bits/stdc++.h> 
#include <eigen3/Eigen/Dense>

using namespace std;

class FEA2{

    public:
        
        // Frame
        int nFrameId;

        // Young Modulus [Pa]
        float E;

        // Poisson Coefficient
        float nu;

        // Behaviour Matrix
        vector<vector<float> > D;

        // Lamé parameters
        float lambda = 0.0;
        float G = 0.0;

        // Gauss Points
        float fg = 0.0;
        vector<vector<float> > gs;

        // Element depth
        float h = 0.0;

        // Stiffness matrix size
        int Ksize = 0;
        
        // Stiffness matrix
        vector<vector<float> > K;
        vector<vector<float> > K1;

        // Displacements
        vector<vector<float> > vU;

        // Forces
        vector<vector<float> > vF;

        // Stiffness matrix determinant
        float DetK = 0.0;

        // Force vector
        vector<vector<float> > f;

        // Debug mode
        bool bDebugMode = false;


        // Constructor & Destructor
        FEA2(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug);
        ~FEA2();

        // Compute elemental stiffness matrix
        vector<vector<float> > ComputeKeiC3D6(vector<vector<float> > vfPts);
        vector<vector<float> > ComputeKeiC3D8(vector<vector<float> > vfPts);

        // Assembly stiffness matrix
        vector<vector<float> > MatrixAssemblyC3D6(vector<vector<float> > vpts, vector<vector<int> > vElts);
        vector<vector<float> > MatrixAssemblyC3D8(vector<vector<float> > vpts, vector<vector<int> > vElts); 

        // Impose Dirichlet conditions (Encastre)
        void ImposeDirichletEncastre(vector<vector<int> > vD, float Klarge);
        
        // Multiply matrices
        vector<vector<float> > MultiplyMatricesEigen(vector<vector<float> > m1, vector<vector<float> > m2);
        
        // Invert matrix
        vector<vector<float> > InvertMatrixEigen(vector<vector<float> > m1);


        // Aux
        void put_to_file_vvf (string outputPath, string delimiter, vector<vector<float> > vvoutput, bool append);
        vector<vector<float> > vector_resize_cols(vector<vector<float> > v1, int n);

};

