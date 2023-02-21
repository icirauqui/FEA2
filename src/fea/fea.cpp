#include "fea.hpp"



FEA2::FEA2(unsigned int input_nFrameId, unsigned int input_E, float input_nu, float input_h, float input_fg1, bool bSetDebug):
    nFrameId(input_nFrameId), E(input_E), nu(input_nu), h(input_h), fg(input_fg1), bDebugMode(bSetDebug) {
    
    lambda = ( nu*E ) / ( (1+nu) * (1-2*nu) );
    G = E / ( 2 * (1+nu) );

    D.clear();
    D.push_back(vector<float>({ (lambda+2*G) ,    lambda    ,    lambda    ,     0.0     ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({    lambda    , (lambda+2*G) ,    lambda    ,     0.0     ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({    lambda    ,    lambda    , (lambda+2*G) ,     0.0     ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({      0.0     ,      0.0     ,      0.0     ,      G      ,     0.0     ,     0.0     }));
    D.push_back(vector<float>({      0.0     ,      0.0     ,      0.0     ,     0.0     ,      G      ,     0.0     }));
    D.push_back(vector<float>({      0.0     ,      0.0     ,      0.0     ,     0.0     ,     0.0     ,      G      }));
 
    gs.clear();
    gs.push_back(vector<float>({-fg,-fg,-fg}));
    gs.push_back(vector<float>({+fg,-fg,-fg}));
    gs.push_back(vector<float>({+fg,+fg,-fg}));
    gs.push_back(vector<float>({-fg,+fg,-fg}));
    gs.push_back(vector<float>({-fg,-fg,+fg}));
    gs.push_back(vector<float>({+fg,-fg,+fg}));
    gs.push_back(vector<float>({+fg,+fg,+fg}));
    gs.push_back(vector<float>({-fg,+fg,+fg}));
}


FEA2::~FEA2() {
}






vector<vector<float> > FEA2::ComputeKeiC3D6(vector<vector<float> > vfPts){

    vector<vector<float> > vBtDB = vector<vector<float> >(24,vector<float>(24,0.0));

    for (unsigned int ops=0; ops<gs.size(); ops++)
    {
        float xi = gs[ops][0];
        float eta = gs[ops][1];
        float zeta = gs[ops][2];

        float dN1dxi = -(1 + zeta)/2;    float dN1deta = -(1 + zeta)/2;    float dN1dzeta =  (1-xi-eta)/2;
        float dN2dxi =  (1 + zeta)/2;    float dN2deta =   0.0;            float dN2dzeta =  xi/2;
        float dN3dxi =  0.0;             float dN3deta =  (1 + zeta)/2;    float dN3dzeta =  eta/2;
        float dN4dxi = -(1 - zeta)/2;    float dN4deta = -(1 - zeta)/2;    float dN4dzeta = -(1-xi-eta)/2;
        float dN5dxi =  (1 - zeta)/2;    float dN5deta =   0.0;            float dN5dzeta = -xi/2;
        float dN6dxi =  0.0;             float dN6deta =  (1 - zeta)/2;    float dN6dzeta = -eta/2;

        /*dxdxi*/   float J_00 = dN1dxi*vfPts[0][0]   + dN2dxi*vfPts[1][0]   + dN3dxi*vfPts[2][0]   + dN4dxi*vfPts[3][0]   + dN5dxi*vfPts[4][0]   + dN6dxi*vfPts[5][0];
        /*dydxi*/   float J_01 = dN1dxi*vfPts[0][1]   + dN2dxi*vfPts[1][1]   + dN3dxi*vfPts[2][1]   + dN4dxi*vfPts[3][1]   + dN5dxi*vfPts[4][1]   + dN6dxi*vfPts[5][1];
        /*dzdxi*/   float J_02 = dN1dxi*vfPts[0][2]   + dN2dxi*vfPts[1][2]   + dN3dxi*vfPts[2][2]   + dN4dxi*vfPts[3][2]   + dN5dxi*vfPts[4][2]   + dN6dxi*vfPts[5][2];
        /*dxdeta*/  float J_10 = dN1deta*vfPts[0][0]  + dN2deta*vfPts[1][0]  + dN3deta*vfPts[2][0]  + dN4deta*vfPts[3][0]  + dN5deta*vfPts[4][0]  + dN6deta*vfPts[5][0];
        /*dydeta*/  float J_11 = dN1deta*vfPts[0][1]  + dN2deta*vfPts[1][1]  + dN3deta*vfPts[2][1]  + dN4deta*vfPts[3][1]  + dN5deta*vfPts[4][1]  + dN6deta*vfPts[5][1];
        /*dzdeta*/  float J_12 = dN1deta*vfPts[0][2]  + dN2deta*vfPts[1][2]  + dN3deta*vfPts[2][2]  + dN4deta*vfPts[3][2]  + dN5deta*vfPts[4][2]  + dN6deta*vfPts[5][2];
        /*dxdzeta*/ float J_20 = dN1dzeta*vfPts[0][0] + dN2dzeta*vfPts[1][0] + dN3dzeta*vfPts[2][0] + dN4dzeta*vfPts[3][0] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0];
        /*dydzeta*/ float J_21 = dN1dzeta*vfPts[0][1] + dN2dzeta*vfPts[1][1] + dN3dzeta*vfPts[2][1] + dN4dzeta*vfPts[3][1] + dN5dzeta*vfPts[4][1] + dN6dzeta*vfPts[5][1];
        /*dzdzeta*/ float J_22 = dN1dzeta*vfPts[0][2] + dN2dzeta*vfPts[1][2] + dN3dzeta*vfPts[2][2] + dN4dzeta*vfPts[3][2] + dN5dzeta*vfPts[4][2] + dN6dzeta*vfPts[5][2];

        float Jac = J_00*J_11*J_22 + J_01*J_12*J_20 + J_10*J_21*J_02 - J_20*J_11*J_02 - J_10*J_01*J_22 - J_21*J_12*J_00;

        //cout << "Jac = " << Jac << endl;

        float J1_00 = (+1) * ( (J_11*J_22) - (J_21*J_12) ) / Jac;     float J1_01 = (-1) * ( (J_01*J_22) - (J_21*J_02) ) / Jac;     float J1_02 = (-1) * ( (J_01*J_12) - (J_11*J_02) ) / Jac;
        float J1_10 = (-1) * ( (J_10*J_22) - (J_20*J_12) ) / Jac;     float J1_11 = (-1) * ( (J_00*J_22) - (J_20*J_02) ) / Jac;     float J1_12 = (-1) * ( (J_00*J_12) - (J_10*J_02) ) / Jac;
        float J1_20 = (+1) * ( (J_10*J_21) - (J_20*J_11) ) / Jac;     float J1_21 = (-1) * ( (J_00*J_21) - (J_20*J_01) ) / Jac;     float J1_22 = (-1) * ( (J_00*J_11) - (J_10*J_01) ) / Jac;

        float dN1dx = J1_00*dN1dxi + J1_01*dN1deta + J1_02*dN1dzeta;     float dN1dy = J1_10*dN1dxi + J1_11*dN1deta + J1_12*dN1dzeta;     float dN1dz = J1_20*dN1dxi + J1_21*dN1deta + J1_22*dN1dzeta;
        float dN2dx = J1_00*dN2dxi + J1_01*dN2deta + J1_02*dN2dzeta;     float dN2dy = J1_10*dN2dxi + J1_11*dN2deta + J1_12*dN2dzeta;     float dN2dz = J1_20*dN2dxi + J1_21*dN2deta + J1_22*dN2dzeta;
        float dN3dx = J1_00*dN3dxi + J1_01*dN3deta + J1_02*dN3dzeta;     float dN3dy = J1_10*dN3dxi + J1_11*dN3deta + J1_12*dN3dzeta;     float dN3dz = J1_20*dN3dxi + J1_21*dN3deta + J1_22*dN3dzeta;
        float dN4dx = J1_00*dN4dxi + J1_01*dN4deta + J1_02*dN4dzeta;     float dN4dy = J1_10*dN4dxi + J1_11*dN4deta + J1_12*dN4dzeta;     float dN4dz = J1_20*dN4dxi + J1_21*dN4deta + J1_22*dN4dzeta;
        float dN5dx = J1_00*dN5dxi + J1_01*dN5deta + J1_02*dN5dzeta;     float dN5dy = J1_10*dN5dxi + J1_11*dN5deta + J1_12*dN5dzeta;     float dN5dz = J1_20*dN5dxi + J1_21*dN5deta + J1_22*dN5dzeta;
        float dN6dx = J1_00*dN6dxi + J1_01*dN6deta + J1_02*dN6dzeta;     float dN6dy = J1_10*dN6dxi + J1_11*dN6deta + J1_12*dN6dzeta;     float dN6dz = J1_20*dN6dxi + J1_21*dN6deta + J1_22*dN6dzeta;

        float B[6][18] = {  dN1dx ,  0.0  ,  0.0  , dN2dx ,  0.0  ,  0.0  , dN3dx ,  0.0  ,  0.0  , dN4dx ,  0.0  ,  0.0  , dN5dx ,  0.0  ,  0.0  , dN6dx ,  0.0  ,  0.0  ,
                            0.0  , dN1dy ,  0.0  ,  0.0  , dN2dy ,  0.0  ,  0.0  , dN3dy ,  0.0  ,  0.0  , dN4dy ,  0.0  ,  0.0  , dN5dy ,  0.0  ,  0.0  , dN6dy ,  0.0  ,
                            0.0  ,  0.0  , dN1dz ,  0.0  ,  0.0  , dN2dz ,  0.0  ,  0.0  , dN3dz ,  0.0  ,  0.0  , dN4dz ,  0.0  ,  0.0  , dN5dz ,  0.0  ,  0.0  , dN6dz ,
                            dN1dy , dN1dx ,  0.0  , dN2dy , dN2dx ,  0.0  , dN3dy , dN3dx ,  0.0  , dN4dy , dN4dx ,  0.0  , dN5dy , dN5dx ,  0.0  , dN6dy , dN6dx ,  0.0  ,
                            dN1dz ,  0.0  , dN1dx , dN2dz ,  0.0  , dN2dx , dN3dz ,  0.0  , dN3dx , dN4dz ,  0.0  , dN4dx , dN5dz ,  0.0  , dN5dx , dN6dz ,  0.0  , dN6dx ,
                            0.0  , dN1dz , dN1dy ,  0.0  , dN2dz , dN2dy ,  0.0  , dN3dz , dN3dy ,  0.0  , dN4dz , dN4dy ,  0.0  , dN5dz , dN5dy ,  0.0  , dN6dz , dN6dy  };

        float BtD[18][6] = {0.0};

        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<6; j++)
                BtD[i][j] = B[0][i]*D[0][j] + B[1][i]*D[1][j] + B[2][i]*D[2][j] + B[3][i]*D[3][j] + B[4][i]*D[4][j] + B[5][i]*D[5][j];

        for (unsigned int i=0; i<18; i++)
            for (unsigned int j=0; j<18; j++){
                float vBtDBaux = BtD[i][0]*B[0][j] + BtD[i][1]*B[1][j] + BtD[i][2]*B[2][j] + BtD[i][3]*B[3][j] + BtD[i][4]*B[4][j] + BtD[i][5]*B[5][j];
                vBtDB[i][j] += vBtDBaux * Jac;
            }
    
    }

    return vBtDB;
}


vector<vector<float> > FEA2::ComputeKeiC3D8(vector<vector<float> > vfPts){

    vector<vector<float> > vBtDB = vector<vector<float> >(24,vector<float>(24,0.0));

    for (unsigned int ops=0; ops<gs.size(); ops++)
    {
        float xi = gs[ops][0];
        float eta = gs[ops][1];
        float zeta = gs[ops][2];

        float dN1dxi = -0.125 * ( (1-eta) * (1-zeta) );     float dN1deta = -0.125 * ( (1-xi) * (1-zeta) );     float dN1dzeta = -0.125 * ( (1-xi) * (1-eta) ); 
        float dN2dxi = +0.125 * ( (1-eta) * (1-zeta) );     float dN2deta = -0.125 * ( (1+xi) * (1-zeta) );     float dN2dzeta = -0.125 * ( (1+xi) * (1-eta) ); 
        float dN3dxi = +0.125 * ( (1+eta) * (1-zeta) );     float dN3deta = +0.125 * ( (1+xi) * (1-zeta) );     float dN3dzeta = -0.125 * ( (1+xi) * (1+eta) ); 
        float dN4dxi = -0.125 * ( (1+eta) * (1-zeta) );     float dN4deta = +0.125 * ( (1-xi) * (1-zeta) );     float dN4dzeta = -0.125 * ( (1-xi) * (1+eta) ); 
        float dN5dxi = -0.125 * ( (1-eta) * (1+zeta) );     float dN5deta = -0.125 * ( (1-xi) * (1+zeta) );     float dN5dzeta = +0.125 * ( (1-xi) * (1-eta) ); 
        float dN6dxi = +0.125 * ( (1-eta) * (1+zeta) );     float dN6deta = -0.125 * ( (1+xi) * (1+zeta) );     float dN6dzeta = +0.125 * ( (1+xi) * (1-eta) ); 
        float dN7dxi = +0.125 * ( (1+eta) * (1+zeta) );     float dN7deta = +0.125 * ( (1+xi) * (1+zeta) );     float dN7dzeta = +0.125 * ( (1+xi) * (1+eta) ); 
        float dN8dxi = -0.125 * ( (1+eta) * (1+zeta) );     float dN8deta = +0.125 * ( (1-xi) * (1+zeta) );     float dN8dzeta = +0.125 * ( (1-xi) * (1+eta) ); 
                
        /*dxdxi*/   float J_00 = dN1dxi*vfPts[0][0] + dN2dxi*vfPts[1][0] + dN3dxi*vfPts[2][0] + dN4dxi*vfPts[3][0] + dN5dxi*vfPts[4][0] + dN6dxi*vfPts[5][0] + dN7dxi*vfPts[6][0] + dN8dxi*vfPts[7][0];
        /*dydxi*/   float J_01 = dN1dxi*vfPts[0][1] + dN2dxi*vfPts[1][1] + dN3dxi*vfPts[2][1] + dN4dxi*vfPts[3][1] + dN5dxi*vfPts[4][1] + dN6dxi*vfPts[5][1] + dN7dxi*vfPts[6][1] + dN8dxi*vfPts[7][1];
        /*dzdxi*/   float J_02 = dN1dxi*vfPts[0][2] + dN2dxi*vfPts[1][2] + dN3dxi*vfPts[2][2] + dN4dxi*vfPts[3][2] + dN5dxi*vfPts[4][2] + dN6dxi*vfPts[5][2] + dN7dxi*vfPts[6][2] + dN8dxi*vfPts[7][2];
        /*dxdeta*/  float J_10 = dN1deta*vfPts[0][0] + dN2deta*vfPts[1][0] + dN3deta*vfPts[2][0] + dN4deta*vfPts[3][0] + dN5deta*vfPts[4][0] + dN6deta*vfPts[5][0] + dN7deta*vfPts[6][0] + dN8deta*vfPts[7][0];
        /*dydeta*/  float J_11 = dN1deta*vfPts[0][1] + dN2deta*vfPts[1][1] + dN3deta*vfPts[2][1] + dN4deta*vfPts[3][1] + dN5deta*vfPts[4][1] + dN6deta*vfPts[5][1] + dN7deta*vfPts[6][1] + dN8deta*vfPts[7][1];
        /*dzdeta*/  float J_12 = dN1deta*vfPts[0][2] + dN2deta*vfPts[1][2] + dN3deta*vfPts[2][2] + dN4deta*vfPts[3][2] + dN5deta*vfPts[4][2] + dN6deta*vfPts[5][2] + dN7deta*vfPts[6][2] + dN8deta*vfPts[7][2];
        /*dxdzeta*/ float J_20 = dN1dzeta*vfPts[0][0] + dN2dzeta*vfPts[1][0] + dN3dzeta*vfPts[2][0] + dN4dzeta*vfPts[3][0] + dN5dzeta*vfPts[4][0] + dN6dzeta*vfPts[5][0] + dN7dzeta*vfPts[6][0] + dN8dzeta*vfPts[7][0];
        /*dydzeta*/ float J_21 = dN1dzeta*vfPts[0][1] + dN2dzeta*vfPts[1][1] + dN3dzeta*vfPts[2][1] + dN4dzeta*vfPts[3][1] + dN5dzeta*vfPts[4][1] + dN6dzeta*vfPts[5][1] + dN7dzeta*vfPts[6][1] + dN8dzeta*vfPts[7][1];
        /*dzdzeta*/ float J_22 = dN1dzeta*vfPts[0][2] + dN2dzeta*vfPts[1][2] + dN3dzeta*vfPts[2][2] + dN4dzeta*vfPts[3][2] + dN5dzeta*vfPts[4][2] + dN6dzeta*vfPts[5][2] + dN7dzeta*vfPts[6][2] + dN8dzeta*vfPts[7][2];

        float Jac = J_00*J_11*J_22 + J_01*J_12*J_20 + J_10*J_21*J_02 - J_20*J_11*J_02 - J_10*J_01*J_22 - J_21*J_12*J_00;

        float J1_00 = (+1) * ( (J_11*J_22) - (J_21*J_12) ) / Jac;     float J1_01 = (-1) * ( (J_01*J_22) - (J_21*J_02) ) / Jac;     float J1_02 = (-1) * ( (J_01*J_12) - (J_11*J_02) ) / Jac;
        float J1_10 = (-1) * ( (J_10*J_22) - (J_20*J_12) ) / Jac;     float J1_11 = (-1) * ( (J_00*J_22) - (J_20*J_02) ) / Jac;     float J1_12 = (-1) * ( (J_00*J_12) - (J_10*J_02) ) / Jac;
        float J1_20 = (+1) * ( (J_10*J_21) - (J_20*J_11) ) / Jac;     float J1_21 = (-1) * ( (J_00*J_21) - (J_20*J_01) ) / Jac;     float J1_22 = (-1) * ( (J_00*J_11) - (J_10*J_01) ) / Jac;

        float dN1dx = J1_00*dN1dxi + J1_01*dN1deta + J1_02*dN1dzeta;     float dN1dy = J1_10*dN1dxi + J1_11*dN1deta + J1_12*dN1dzeta;     float dN1dz = J1_20*dN1dxi + J1_21*dN1deta + J1_22*dN1dzeta;
        float dN2dx = J1_00*dN2dxi + J1_01*dN2deta + J1_02*dN2dzeta;     float dN2dy = J1_10*dN2dxi + J1_11*dN2deta + J1_12*dN2dzeta;     float dN2dz = J1_20*dN2dxi + J1_21*dN2deta + J1_22*dN2dzeta;
        float dN3dx = J1_00*dN3dxi + J1_01*dN3deta + J1_02*dN3dzeta;     float dN3dy = J1_10*dN3dxi + J1_11*dN3deta + J1_12*dN3dzeta;     float dN3dz = J1_20*dN3dxi + J1_21*dN3deta + J1_22*dN3dzeta;
        float dN4dx = J1_00*dN4dxi + J1_01*dN4deta + J1_02*dN4dzeta;     float dN4dy = J1_10*dN4dxi + J1_11*dN4deta + J1_12*dN4dzeta;     float dN4dz = J1_20*dN4dxi + J1_21*dN4deta + J1_22*dN4dzeta;
        float dN5dx = J1_00*dN5dxi + J1_01*dN5deta + J1_02*dN5dzeta;     float dN5dy = J1_10*dN5dxi + J1_11*dN5deta + J1_12*dN5dzeta;     float dN5dz = J1_20*dN5dxi + J1_21*dN5deta + J1_22*dN5dzeta;
        float dN6dx = J1_00*dN6dxi + J1_01*dN6deta + J1_02*dN6dzeta;     float dN6dy = J1_10*dN6dxi + J1_11*dN6deta + J1_12*dN6dzeta;     float dN6dz = J1_20*dN6dxi + J1_21*dN6deta + J1_22*dN6dzeta;
        float dN7dx = J1_00*dN7dxi + J1_01*dN7deta + J1_02*dN7dzeta;     float dN7dy = J1_10*dN7dxi + J1_11*dN7deta + J1_12*dN7dzeta;     float dN7dz = J1_20*dN7dxi + J1_21*dN7deta + J1_22*dN7dzeta;
        float dN8dx = J1_00*dN8dxi + J1_01*dN8deta + J1_02*dN8dzeta;     float dN8dy = J1_10*dN8dxi + J1_11*dN8deta + J1_12*dN8dzeta;     float dN8dz = J1_20*dN8dxi + J1_21*dN8deta + J1_22*dN8dzeta;

        float B[6][24] = {  dN1dx ,  0.0  ,  0.0  , dN2dx ,  0.0  ,  0.0  , dN3dx ,  0.0  ,  0.0  , dN4dx ,  0.0  ,  0.0  , dN5dx ,  0.0  ,  0.0  , dN6dx ,  0.0  ,  0.0  , dN7dx ,  0.0  ,  0.0  , dN8dx ,  0.0  ,  0.0  ,
                            0.0  , dN1dy ,  0.0  ,  0.0  , dN2dy ,  0.0  ,  0.0  , dN3dy ,  0.0  ,  0.0  , dN4dy ,  0.0  ,  0.0  , dN5dy ,  0.0  ,  0.0  , dN6dy ,  0.0  ,  0.0  , dN7dy ,  0.0  ,  0.0  , dN8dy ,  0.0  ,
                            0.0  ,  0.0  , dN1dz ,  0.0  ,  0.0  , dN2dz ,  0.0  ,  0.0  , dN3dz ,  0.0  ,  0.0  , dN4dz ,  0.0  ,  0.0  , dN5dz ,  0.0  ,  0.0  , dN6dz ,  0.0  ,  0.0  , dN7dz ,  0.0  ,  0.0  , dN8dz ,
                            dN1dy , dN1dx ,  0.0  , dN2dy , dN2dx ,  0.0  , dN3dy , dN3dx ,  0.0  , dN4dy , dN4dx ,  0.0  , dN5dy , dN5dx ,  0.0  , dN6dy , dN6dx ,  0.0  , dN7dy , dN7dx ,  0.0  , dN8dy , dN8dx ,  0.0  ,
                            dN1dz ,  0.0  , dN1dx , dN2dz ,  0.0  , dN2dx , dN3dz ,  0.0  , dN3dx , dN4dz ,  0.0  , dN4dx , dN5dz ,  0.0  , dN5dx , dN6dz ,  0.0  , dN6dx , dN7dz ,  0.0  , dN7dx , dN8dz ,  0.0  , dN8dx ,
                            0.0  , dN1dz , dN1dy ,  0.0  , dN2dz , dN2dy ,  0.0  , dN3dz , dN3dy ,  0.0  , dN4dz , dN4dy ,  0.0  , dN5dz , dN5dy ,  0.0  , dN6dz , dN6dy ,  0.0  , dN7dz , dN7dy ,  0.0  , dN8dz , dN8dy  };

        float BtD[24][6] = {0.0};

        for (unsigned int i=0; i<24; i++)
            for (unsigned int j=0; j<6; j++)
                BtD[i][j] = B[0][i]*D[0][j] + B[1][i]*D[1][j] + B[2][i]*D[2][j] + B[3][i]*D[3][j] + B[4][i]*D[4][j] + B[5][i]*D[5][j];

        for (unsigned int i=0; i<24; i++)
            for (unsigned int j=0; j<24; j++){
                float vBtDBaux = BtD[i][0]*B[0][j] + BtD[i][1]*B[1][j] + BtD[i][2]*B[2][j] + BtD[i][3]*B[3][j] + BtD[i][4]*B[4][j] + BtD[i][5]*B[5][j];
                vBtDB[i][j] += vBtDBaux * Jac;
            }
    }

    return vBtDB;
}


vector<vector<float> > FEA2::MatrixAssemblyC3D6(vector<vector<float> > vpts, vector<vector<int> > vElts) {
    vector<vector<float> > Kt = vector<vector<float> >(3*vpts.size(),vector<float>(3*vpts.size(),0.0));

    for (unsigned int i=0; i<vElts.size(); i++){
        vector<vector<float> > xyzi;
        vector<int> nodesi;
        for (unsigned int j=0; j<vElts[i].size(); j++){
            int node = vElts[i][j];
            nodesi.push_back(node);
            xyzi.push_back(vpts[node]);
        }

        vector<vector<float> > Kei = ComputeKeiC3D6(xyzi);

        /*
        for (unsigned int j=0; j<Kei.size(); j++){
            for (unsigned int k=0; k<Kei[j].size(); k++){
                cout << Kei[j][k] << " ";
            }
            cout << endl;
        }
        */

        vector<int> mn;
        for (unsigned j=0; j<nodesi.size(); j++) {
            mn.push_back(nodesi[j]*3);
        }

        for (unsigned int ni=0; ni<nodesi.size(); ni++){
            for (unsigned int nj=0; nj<nodesi.size(); nj++){
                for (unsigned int m=0; m<3; m++){
                    for (unsigned int n=0; n<3; n++){
                        Kt[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                    }
                }
            }
        }
    }
    
    return Kt;
}


vector<vector<float> > FEA2::MatrixAssemblyC3D8(vector<vector<float> > vpts, vector<vector<int> > vElts) {
    vector<vector<float> > Kt = vector<vector<float> >(3*vpts.size(),vector<float>(3*vpts.size(),0.0));

    for (unsigned int i=0; i<vElts.size(); i++){
        vector<vector<float> > xyzi;
        vector<int> nodesi;
        for (unsigned int j=0; j<vElts[i].size(); j++){
            int node = vElts[i][j];
            nodesi.push_back(node);
            xyzi.push_back(vpts[node]);
        }

        vector<vector<float> > Kei = ComputeKeiC3D8(xyzi);

        vector<int> mn;
        for (unsigned j=0; j<nodesi.size(); j++) {
            mn.push_back(nodesi[j]*3);
        }

        for (unsigned int ni=0; ni<nodesi.size(); ni++){
            for (unsigned int nj=0; nj<nodesi.size(); nj++){
                for (unsigned int m=0; m<3; m++){
                    for (unsigned int n=0; n<3; n++){
                        Kt[mn[ni]+m][mn[nj]+n] += Kei[3*ni+m][3*nj+n];
                    }
                }
            }
        }
    }
    
    return Kt;
}


void FEA2::ImposeDirichletEncastre(vector<vector<int> > vD, float Klarge){
    for (unsigned int i=0; i<vD.size(); i++){
        int mp0 = 3*(vD[i][0] - 1);
        int mp1 = mp0 + 1;
        int mp2 = mp0 + 2;

        K[mp0][mp0] = Klarge;
        K[mp1][mp1] = Klarge;
        K[mp2][mp2] = Klarge;

        vF[mp0][0] = 0.0;
        vF[mp1][0] = 0.0;
        vF[mp2][0] = 0.0;
    }
}


vector<vector<float> > FEA2::MultiplyMatricesEigen(vector<vector<float> > m1, vector<vector<float> > m2) {
    vector<vector<float> > pr = vector<vector<float> >(m1.size(),vector<float>(m2[0].size(),0.0));

    if (m1[0].size() != m2.size()) {
        pr = vector<vector<float> >(1,vector<float>(1,0.0));
        return pr;
    }

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> A;
	A.resize(m1.size(),m1[0].size());
	for (unsigned int i=0; i<m1.size(); i++){
		for (unsigned int j=0; j<m1[i].size(); j++){
			A(i,j) = m1[i][j];
		}
	}

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> B;
	B.resize(m2.size(),m2[0].size());
	for (unsigned int i=0; i<m2.size(); i++){
		for (unsigned int j=0; j<m2[i].size(); j++){
			B(i,j) = m2[i][j];
		}
	}

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> C;
	C.resize(pr.size(),pr[0].size());
    C = A*B;

    for (unsigned int i=0; i<pr.size(); i++){
        for (unsigned int j=0; j<pr[i].size(); j++){
            pr[i][j] = C(i,j);
        }
    }

    return pr;
}


vector<vector<float> > FEA2::InvertMatrixEigen(vector<vector<float> > m1) { 
	Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> A;
	A.resize(m1.size(),m1.size());

	for (unsigned int i=0; i<m1.size(); i++){
		for (unsigned int j=0; j<m1[i].size(); j++){
			A(i,j) = m1[i][j];
		}
	}

    float detA = A.determinant();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> A1 = A.inverse();
    
    vector<vector<float> > Ainv;
    if (detA!=0.0){
        for (unsigned int i=0; i<m1.size(); i++){
            vector<float> Ainvi;
            for (unsigned int j=0; j<m1[i].size(); j++){
                Ainvi.push_back(A1(i,j));
            }
            Ainv.push_back(Ainvi);
        }
    }
    else{
        vector<float> Ainvi;
        Ainvi.push_back(0.0);
        Ainv.push_back(Ainvi);
    }

    return Ainv;
} 


void FEA2::put_to_file_vvf (string outputPath, string delimiter, vector<vector<float> > vvoutput, bool append){
    
    ofstream os;
    if (append)
        os.open(outputPath.c_str(), ios::out | ios::app );
    else
        os.open(outputPath.c_str(), ios::out);
    
    for (unsigned int i=0; i<vvoutput.size(); i++) {
        for (unsigned int j=0; j<vvoutput[i].size(); j++) {
            os << vvoutput[i][j] << " ";
        }
        os << endl;
    }
    os << endl;
    
    os.close();
}


vector<vector<float> > FEA2::vector_resize_cols(vector<vector<float> > v1, int n){
    vector<vector<float> > v2;
    vector<float> v2i;
    for (unsigned int i=0; i<v1.size(); i++){
        for (unsigned int j=0; j<v1[i].size(); j++){
            v2i.push_back(v1[i][j]);
            if (v2i.size() == n){
                v2.push_back(v2i);
                v2i.clear();
            }
        }
    }
    return v2;
}


