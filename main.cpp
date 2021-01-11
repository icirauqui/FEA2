#include <iostream>
#include <vector>
#include <math.h>

#include "aux.h"
#include "FEA2.h"


using namespace std;

float E = 3500.0;
float nu = 0.495;
float fg = 0.577350269;
float Klarge = 100000000.0;

vector<vector<float> > vpts;
vector<vector<int> > vElts;
vector<vector<float> > vvF;
vector<vector<int> > vDir;




int main() {

    vpts = get_from_file_vvf("input_C3D6/input_points.csv",",");
    vElts = get_from_file_vvn("input_C3D6/input_elements.csv",",");
    vvF = get_from_file_vvf("input_C3D6/input_forces.csv",",");
    vDir = get_from_file_vvn("input_C3D6/input_Dirichlet.csv",",");

    // Node numbers in Abaqus start at 1, change to 0 for cpp
    for (unsigned int i=0; i<vElts.size(); i++){
        for (unsigned int j=0; j<vElts[i].size(); j++){
            vElts[i][j] -= 1;
        }
    }

    FEA2 feaC3D6(1,E,nu,1,fg,false);

    feaC3D6.K = feaC3D6.MatrixAssemblyC3D6(vpts,vElts);
    feaC3D6.vF = vector_resize_cols(vvF,1);
    
    feaC3D6.ImposeDirichletEncastre(vDir,Klarge);
    
    feaC3D6.K1 = feaC3D6.InvertMatrixEigen(feaC3D6.K);
    feaC3D6.vU = feaC3D6.MultiplyMatricesEigen(feaC3D6.K1,feaC3D6.vF);

    vector<vector<float> > vvU = vector_resize_cols(feaC3D6.vU,3);

    vector<vector<float> > vUt;
    vector<float> vUti;
    for (unsigned int i=0; i<feaC3D6.vU.size(); i++){
        vUti.push_back(feaC3D6.vU[i][0]);
    }
    vUt.push_back(vUti);

    vector<vector<float> > sE = feaC3D6.MultiplyMatricesEigen(vUt,feaC3D6.vF);
    cout << "Strain Energy C3D6 = " << sE[0][0] << " Jules" << endl;

    put_to_file_vvf("output_C3D6/vvK.csv",",",feaC3D6.K, false);
    put_to_file_vvf("output_C3D6/vvK1.csv",",",feaC3D6.K1,false);
    put_to_file_vvf("output_C3D6/vvU.csv",",",vvU,false);






    vpts = get_from_file_vvf("input_C3D8/input_points.csv",",");
    vElts = get_from_file_vvn("input_C3D8/input_elements.csv",",");
    vvF = get_from_file_vvf("input_C3D8/input_forces.csv",",");
    vDir = get_from_file_vvn("input_C3D8/input_Dirichlet.csv",",");

    // Node numbers in Abaqus start at 1, change to 0 for cpp
    for (unsigned int i=0; i<vElts.size(); i++){
        for (unsigned int j=0; j<vElts[i].size(); j++){
            vElts[i][j] -= 1;
        }
    }

    FEA2 fea(1,E,nu,1,fg,false);

    //fea.MatrixAssembly(vpts,vElts);
    fea.K = fea.MatrixAssemblyC3D8(vpts,vElts);
    fea.vF = vector_resize_cols(vvF,1);
    
    fea.ImposeDirichletEncastre(vDir,Klarge);

    fea.K1 = fea.InvertMatrixEigen(fea.K);
    fea.vU = fea.MultiplyMatricesEigen(fea.K1,fea.vF);

    vector<vector<float> > vvU2 = vector_resize_cols(fea.vU,3);

    vector<vector<float> > vUt2;
    vector<float> vUti2;
    for (unsigned int i=0; i<fea.vU.size(); i++){
        vUti2.push_back(fea.vU[i][0]);
    }
    vUt2.push_back(vUti2);

    vector<vector<float> > sE2 = fea.MultiplyMatricesEigen(vUt2,fea.vF);
    cout << "Strain Energy C3D8 = " << sE2[0][0] << " Jules" << endl;

    put_to_file_vvf("output_C3D8/vvK.csv",",",fea.K, false);
    put_to_file_vvf("output_C3D8/vvK1.csv",",",fea.K1,false);
    put_to_file_vvf("output_C3D8/vvU.csv",",",vvU2,false);

    return 0;
}
