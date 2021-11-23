#include <blockInfo.h>

void gBlockInfo::formCmat()
{
    // cmat  {};

    double YoungMod = materialProperty[1];
    double PoissonRatio = materialProperty[2];
    double LamdaM = YoungMod * PoissonRatio / ((1 + PoissonRatio) * (1 - 2 * PoissonRatio));
    double MuM = 0.5 * YoungMod / (1 + PoissonRatio);
    double  LamdaPlus2Mu = 2 * MuM + LamdaM;
    materialProperty[18]=LamdaM;
    materialProperty[19]=MuM;


    cmat.val[0][0] = LamdaPlus2Mu;
    cmat.val[1][1] = LamdaPlus2Mu;
    cmat.val[2][2] = LamdaPlus2Mu;
    cmat.val[3][3] = MuM;
    cmat.val[4][4] = MuM;
    cmat.val[5][5] = MuM;
    cmat.val[0][1] = LamdaM;
    cmat.val[1][0] = LamdaM;
    cmat.val[0][2] = LamdaM;
    cmat.val[2][0] = LamdaM;
    cmat.val[1][2] = LamdaM;
    cmat.val[2][1] = LamdaM;

}