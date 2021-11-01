#include "rkpmCudaMethod.h"

__global__ void
findBin(const cellPosition *positionDev, int *nodeBinId, int nc, const gmbinInfo binInfo, const domainBound modelBound)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < nc)
    {
        for (int dm = 0; dm < 3; dm++)
        {
            double localCoo = positionDev[index].coo[dm];
            if (localCoo > modelBound.max[dm])
            {
                nodeBinId[3 * index + dm] = binInfo.num[dm] - 1;
            }
            else if (localCoo < modelBound.min[dm])
            {
                nodeBinId[3 * index + dm] = 0;
            }
            else
            {
                nodeBinId[3 * index + dm] = __double2int_rd((localCoo - modelBound.min[dm]) / binInfo.size);
            }
        }

        // nodeBinId[index]=nid[0]+nid[1]*binInfo.num[0]+nid[2]*binInfo.num[0]*binInfo.num[1];
    }
}

__global__ void cleanBin(int nb, gmBinNode *binNode)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < nb)
    {
        binNode[index].nodeNum = 0;
    }
}
__global__ void countBin(const int *nodeBinId, const int nc, gmBinNode *binNode, const gmbinInfo binInfo)
{
    for (int index = 0; index < nc; index++)
    {
        int i, j, k, binId;
        i = nodeBinId[3 * index];
        j = nodeBinId[3 * index + 1];
        k = nodeBinId[3 * index + 2];
        binId = i + j * binInfo.num[0] + k * binInfo.num[0] * binInfo.num[1];
        binNode[binId].nodeNum += 1;
        binNode[binId].nodeId[binNode[binId].nodeNum - 1] = index;
    }
}

__device__ void L2Norm(const double vec[3], double *norm)
{
    *norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

__device__ void hardSearch(const cellPosition *localCellPosition, gmNodeNeighbor *localNodeNeighbor, const cellPosition *positionDev, const gmBinNode *localBinNode)
{
    for (int index = 0; index < (*localBinNode).nodeNum; index++)
    {
        int comId = (*localBinNode).nodeId[index];
        double vec[3]{}, dist{};
        double *dptr = &dist;
        for (int i = 0; i < 3; i++)
            vec[i] = positionDev[comId].coo[i] - (*localCellPosition).coo[i];
        // L2Norm(vec,dptr);

        if ((abs(vec[0]) < positionDev[comId].win) && (abs(vec[1]) < positionDev[comId].win) && (abs(vec[2]) < positionDev[comId].win))
        {
            (*localNodeNeighbor).neighborId[(*localNodeNeighbor).neighborNum] = comId;
            (*localNodeNeighbor).neighborNum++;
        }
    }
}

__global__ void nodeNeighborSearch(const gmBinNode *binNode, const int nc, gmNodeNeighbor *nodeNeighbor, const gmbinInfo binInfo, int *nodeBinId, const cellPosition *positionDev)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < nc)
    {

        cellPosition localCellPosition;
        localCellPosition = positionDev[index];

        int XId[3], YId[3], ZId[3];
        for (int i = 0; i < 3; i++)
        {
            XId[i] = nodeBinId[3 * index] - 1 + i;
            YId[i] = nodeBinId[3 * index + 1] - 1 + i;
            ZId[i] = nodeBinId[3 * index + 2] - 1 + i;
        }
        for (int i = 0; i < 3; i++)
        {
            if ((XId[i] >= 0) && (XId[i] < binInfo.num[0]))
            {
                for (int j = 0; j < 3; j++)
                {
                    if ((YId[j] >= 0) && (YId[j] < binInfo.num[1]))
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if ((ZId[k] >= 0) && (ZId[k] < binInfo.num[2]))
                            {
                                int localBinId = XId[i] + YId[j] * binInfo.num[0] + ZId[k] * binInfo.num[1] * binInfo.num[0];
                                hardSearch(positionDev + index, nodeNeighbor + index, positionDev, binNode + localBinId);
                            }
                        }
                    }
                }
            }
        }
    }
}

__global__ void updateRK(int nc, gmNodeNeighbor *nodeNeighbor, const cellPosition *positionDev, gmShape *shape, gmShapeGradient *shapeGradient)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < nc)
    {
        double localCoo[3];
        for (int i = 0; i < 3; i++)
            localCoo[i] = positionDev[index].coo[i];

        mat4 M, dm1, dm2, dm3;
        // mat4 *Mptr=&M;
        // mat4 *dm1ptr=&dm1;
        // mat4 *dm2ptr=&dm2;
        // mat4 *dm3ptr=&dm3;
        // mat4Clear(Mptr);
        // mat4Clear(dm1ptr);
        // mat4Clear(dm2ptr);
        // mat4Clear(dm3ptr);

        for (int i = 0; i < nodeNeighbor[index].neighborNum; i++)
        {
            vec4 hx;

            int localNeighborId = nodeNeighbor[index].neighborId[i];
            for (int j = 0; j < 3; j++)
                hx.val[j + 1] = localCoo[j] - positionDev[localNeighborId].coo[j];
            hx.val[0] = 1;

            vec4 phi;
            getPhi(&hx, &phi, positionDev[localNeighborId].win);

            mat4 mtemp4;
            vec4Expand(&hx, &hx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&M, &mtemp4);

            vec4 dhx{};
            dhx.val[1] = 1;
            vec4Expand(&dhx, &hx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&dm1, &mtemp4);
            vec4Expand(&hx, &dhx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&dm1, &mtemp4);
            vec4Expand(&hx, &hx, &mtemp4);
            mat4Scale(phi.val[1], &mtemp4);
            mat4Increase(&dm1, &mtemp4);

            dhx.val[1] = 0;
            dhx.val[2] = 1;
            vec4Expand(&dhx, &hx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&dm2, &mtemp4);
            vec4Expand(&hx, &dhx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&dm2, &mtemp4);
            vec4Expand(&hx, &hx, &mtemp4);
            mat4Scale(phi.val[2], &mtemp4);
            mat4Increase(&dm2, &mtemp4);

            dhx.val[2] = 0;
            dhx.val[3] = 1;
            vec4Expand(&dhx, &hx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&dm3, &mtemp4);
            vec4Expand(&hx, &dhx, &mtemp4);
            mat4Scale(phi.val[0], &mtemp4);
            mat4Increase(&dm3, &mtemp4);
            vec4Expand(&hx, &hx, &mtemp4);
            mat4Scale(phi.val[3], &mtemp4);
            mat4Increase(&dm3, &mtemp4);
        }
        mat4 invM, invdm1, invdm2, invdm3;
        mat4 mtemp4;

        mat4Invert(&M, &invM);

        mat4Prod(&invM, &dm1, &mtemp4);
        mat4Scale(-1, &mtemp4);
        mat4Prod(&mtemp4, &invM, &invdm1);

        mat4Prod(&invM, &dm2, &mtemp4);
        mat4Scale(-1, &mtemp4);
        mat4Prod(&mtemp4, &invM, &invdm2);


        mat4Prod(&invM, &dm3, &mtemp4);
        mat4Scale(-1, &mtemp4);
        mat4Prod(&mtemp4, &invM, &invdm3);

        vec4 b0, b1, b2, b3;

        for (int i = 0; i < 4; i++)
        {
            b0.val[i] = invM.val[0][i];
            b1.val[i] = invdm1.val[0][i];
            b2.val[i] = invdm2.val[0][i];
            b3.val[i] = invdm3.val[0][i];
        }
        for (int i = 0; i < nodeNeighbor[index].neighborNum; i++)
        {
            vec4 hx;

            int localNeighborId = nodeNeighbor[index].neighborId[i];
            for (int j = 0; j < 3; j++)
                hx.val[j + 1] = localCoo[j] - positionDev[localNeighborId].coo[j];
            hx.val[0] = 1;

            vec4 phi;
            getPhi(&hx, &phi, positionDev[localNeighborId].win);
            double temp, temp1, temp2, temp3;

            vec4Dot(&b0,&hx,&temp);
            shape[index].val[i]=phi.val[0]*temp;

            vec4Dot(&b1,&hx,&temp1);
            vec4Dot(&b2,&hx,&temp2);
            vec4Dot(&b3,&hx,&temp3);
            shapeGradient[index].val[0][i]=phi.val[0]*temp1 + phi.val[0] * b0.val[1] + phi.val[1] * temp;
            shapeGradient[index].val[1][i]=phi.val[0]*temp2 + phi.val[0] * b0.val[2] + phi.val[2] * temp;
            shapeGradient[index].val[2][i]=phi.val[0]*temp3 + phi.val[0] * b0.val[3] + phi.val[3] * temp;

            // shape[index].val[i]=positionDev[localNeighborId].coo[0];            
            // shape[index].val[i]=localCoo[2];           
            // shapeGradient[index].val[0][i]=positionDev[localNeighborId].coo[2];
            // shapeGradient[index].val[0][i]=positionDev[localNeighborId].win;
            // shapeGradient[index].val[0][i]=localCoo[2];
            // shape[index].val[i]=hx.val[1];
            // shapeGradient[index].val[0][i]=hx.val[2];
            // shapeGradient[index].val[1][i]=hx.val[3];
        }

        // for (int i = 0; i < 4; i++)
        //     for (int j = 0; j < 4; j++)
        //     {
        //         shape[index].val[i*4+j]=M.val[i][j];
        //     }
        
    }
}
__device__ void getPhi(vec4 *hxPtr, vec4 *phiPtr, double win)
{
    double localPhi[3], localPhiD[3];
    double zl[3];

    for (int i = 0; i < 3; i++)
    {
        zl[i] = abs((*hxPtr).val[i + 1])/win;
        if (zl[i] <= 0.5)
        {
            localPhi[i] = 2.0 / 3.0 - 4 * zl[i] * zl[i] + 4 * zl[i] * zl[i] * zl[i];
            localPhiD[i] = -8 * zl[i] + 12 * zl[i] * zl[i];
            if ((*hxPtr).val[i + 1] < 0)
                localPhiD[i] = -localPhiD[i];
        }
        else if ((zl[i] >= 0.5) && (zl[i] <= 1))
        {
            localPhi[i] = 4.0 / 3.0 - 4 * zl[i] + 4 * zl[i] * zl[i] - 4.0 / 3.0 * zl[i] * zl[i] * zl[i];
            localPhiD[i] = -4 + 8 * zl[i] - 4 * zl[i] * zl[i];
            if ((*hxPtr).val[i + 1] < 0)
                localPhiD[i] = -localPhiD[i];
        }
        else
        {
            localPhi[i] = 0;
            localPhiD[0] = 0;
            localPhiD[1] = 0;
            localPhiD[2] = 0;
        }
    }
    (*phiPtr).val[0] = localPhi[0] * localPhi[1] * localPhi[2];
    (*phiPtr).val[1] = localPhi[1] * localPhi[2] * localPhiD[0] / win;
    (*phiPtr).val[2] = localPhi[0] * localPhi[2] * localPhiD[1] / win;
    (*phiPtr).val[3] = localPhi[0] * localPhi[1] * localPhiD[2] / win;
}

__device__ void vec4Expand(vec4 *va, vec4 *vb, mat4 *ans)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            (*ans).val[i][j] = (*va).val[i] * (*vb).val[j];
}

__device__ void vec4Dot(vec4 *va, vec4 *vb, double *ans)
{
    double temp {0};
    for (int i = 0; i < 4; i++) temp+=(*va).val[i]*(*vb).val[i];
    (*ans)=temp;
}

__device__ void mat4Scale(const double scaler, mat4 *ans)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            (*ans).val[i][j] = scaler * (*ans).val[i][j];
}
// __device__ void mat4Clear(mat4 *ans)
// {
//     for (int i = 0; i < 4; i++)
//         for (int j = 0; j < 4; j++)
//             (*ans).val[i][j] = 0;
// }

__device__ void mat4Increase(mat4 *ans, mat4 *increment)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            (*ans).val[i][j] = (*ans).val[i][j] + (*increment).val[i][j];
}

__device__ void mat4Prod(mat4 *ma, mat4 *mb, mat4 *ans)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            (*ans).val[i][j] = 0;
            for (int k = 0; k < 4; k++)
                (*ans).val[i][j] += (*ma).val[i][k] * (*mb).val[k][j];
        }
}

// __device__ void vec4Tmat4(vec4 *va, mat4 *ma, vec4 *ans)
// {
//     for (int i = 0; i < 4; i++)
//     {
//         (*ans).val[i]=0;
//         for (int j = 0; j < 4; j++)
//             (*ans).val[i]=0;

//     }

// }

__device__ void mat4Invert(mat4 *ma, mat4 *ans)
{
    double det;
    mat4 ta;
    det = (*ma).val[0][0] * ((*ma).val[1][1] * ((*ma).val[2][2] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][2]) + (*ma).val[1][2] * ((*ma).val[2][3] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][3]) + (*ma).val[1][3] * ((*ma).val[2][1] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][1])) - (*ma).val[0][1] * ((*ma).val[1][0] * ((*ma).val[2][2] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][2]) + (*ma).val[1][2] * ((*ma).val[2][3] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][3]) + (*ma).val[1][3] * ((*ma).val[2][0] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][0])) + (*ma).val[0][2] * ((*ma).val[1][0] * ((*ma).val[2][1] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][1]) + (*ma).val[1][1] * ((*ma).val[2][3] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][3]) + (*ma).val[1][3] * ((*ma).val[2][0] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][0])) - (*ma).val[0][3] * ((*ma).val[1][0] * ((*ma).val[2][1] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][1]) + (*ma).val[1][1] * ((*ma).val[2][2] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][2]) + (*ma).val[1][2] * ((*ma).val[2][0] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][0]));
    ta.val[0][0] = (*ma).val[1][1] * ((*ma).val[2][2] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][2]) + (*ma).val[1][2] * ((*ma).val[2][3] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][3]) + (*ma).val[1][3] * ((*ma).val[2][1] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][1]);
    ta.val[1][0] = (*ma).val[1][0] * ((*ma).val[2][3] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][3]) + (*ma).val[1][2] * ((*ma).val[2][0] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][0]) + (*ma).val[1][3] * ((*ma).val[2][2] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][2]);
    ta.val[2][0] = (*ma).val[1][0] * ((*ma).val[2][1] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][1]) + (*ma).val[1][1] * ((*ma).val[2][3] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][3]) + (*ma).val[1][3] * ((*ma).val[2][0] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][0]);
    ta.val[3][0] = (*ma).val[1][0] * ((*ma).val[2][2] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][2]) + (*ma).val[1][1] * ((*ma).val[2][0] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][0]) + (*ma).val[1][2] * ((*ma).val[2][1] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][1]);
    ta.val[0][1] = (*ma).val[0][1] * ((*ma).val[2][3] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][3]) + (*ma).val[0][2] * ((*ma).val[2][1] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][1]) + (*ma).val[0][3] * ((*ma).val[2][2] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][2]);
    ta.val[1][1] = (*ma).val[0][0] * ((*ma).val[2][2] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][2]) + (*ma).val[0][2] * ((*ma).val[2][3] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][3]) + (*ma).val[0][3] * ((*ma).val[2][0] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][0]);
    ta.val[2][1] = (*ma).val[0][0] * ((*ma).val[2][3] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][3]) + (*ma).val[0][1] * ((*ma).val[2][0] * (*ma).val[3][3] - (*ma).val[2][3] * (*ma).val[3][0]) + (*ma).val[0][3] * ((*ma).val[2][1] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][1]);
    ta.val[3][1] = (*ma).val[0][0] * ((*ma).val[2][1] * (*ma).val[3][2] - (*ma).val[2][2] * (*ma).val[3][1]) + (*ma).val[0][1] * ((*ma).val[2][2] * (*ma).val[3][0] - (*ma).val[2][0] * (*ma).val[3][2]) + (*ma).val[0][2] * ((*ma).val[2][0] * (*ma).val[3][1] - (*ma).val[2][1] * (*ma).val[3][0]);
    ta.val[0][2] = (*ma).val[0][1] * ((*ma).val[1][2] * (*ma).val[3][3] - (*ma).val[1][3] * (*ma).val[3][2]) + (*ma).val[0][2] * ((*ma).val[1][3] * (*ma).val[3][1] - (*ma).val[1][1] * (*ma).val[3][3]) + (*ma).val[0][3] * ((*ma).val[1][1] * (*ma).val[3][2] - (*ma).val[1][2] * (*ma).val[3][1]);
    ta.val[1][2] = (*ma).val[0][0] * ((*ma).val[1][3] * (*ma).val[3][2] - (*ma).val[1][2] * (*ma).val[3][3]) + (*ma).val[0][2] * ((*ma).val[1][0] * (*ma).val[3][3] - (*ma).val[1][3] * (*ma).val[3][0]) + (*ma).val[0][3] * ((*ma).val[1][2] * (*ma).val[3][0] - (*ma).val[1][0] * (*ma).val[3][2]);
    ta.val[2][2] = (*ma).val[0][0] * ((*ma).val[1][1] * (*ma).val[3][3] - (*ma).val[1][3] * (*ma).val[3][1]) + (*ma).val[0][1] * ((*ma).val[1][3] * (*ma).val[3][0] - (*ma).val[1][0] * (*ma).val[3][3]) + (*ma).val[0][3] * ((*ma).val[1][0] * (*ma).val[3][1] - (*ma).val[1][1] * (*ma).val[3][0]);
    ta.val[3][2] = (*ma).val[0][0] * ((*ma).val[1][2] * (*ma).val[3][1] - (*ma).val[1][1] * (*ma).val[3][2]) + (*ma).val[0][1] * ((*ma).val[1][0] * (*ma).val[3][2] - (*ma).val[1][2] * (*ma).val[3][0]) + (*ma).val[0][2] * ((*ma).val[1][1] * (*ma).val[3][0] - (*ma).val[1][0] * (*ma).val[3][1]);
    ta.val[0][3] = (*ma).val[0][1] * ((*ma).val[1][3] * (*ma).val[2][2] - (*ma).val[1][2] * (*ma).val[2][3]) + (*ma).val[0][2] * ((*ma).val[1][1] * (*ma).val[2][3] - (*ma).val[1][3] * (*ma).val[2][1]) + (*ma).val[0][3] * ((*ma).val[1][2] * (*ma).val[2][1] - (*ma).val[1][1] * (*ma).val[2][2]);
    ta.val[1][3] = (*ma).val[0][0] * ((*ma).val[1][2] * (*ma).val[2][3] - (*ma).val[1][3] * (*ma).val[2][2]) + (*ma).val[0][2] * ((*ma).val[1][3] * (*ma).val[2][0] - (*ma).val[1][0] * (*ma).val[2][3]) + (*ma).val[0][3] * ((*ma).val[1][0] * (*ma).val[2][2] - (*ma).val[1][2] * (*ma).val[2][0]);
    ta.val[2][3] = (*ma).val[0][0] * ((*ma).val[1][3] * (*ma).val[2][1] - (*ma).val[1][1] * (*ma).val[2][3]) + (*ma).val[0][1] * ((*ma).val[1][0] * (*ma).val[2][3] - (*ma).val[1][3] * (*ma).val[2][0]) + (*ma).val[0][3] * ((*ma).val[1][1] * (*ma).val[2][0] - (*ma).val[1][0] * (*ma).val[2][1]);
    ta.val[3][3] = (*ma).val[0][0] * ((*ma).val[1][1] * (*ma).val[2][2] - (*ma).val[1][2] * (*ma).val[2][1]) + (*ma).val[0][1] * ((*ma).val[1][2] * (*ma).val[2][0] - (*ma).val[1][0] * (*ma).val[2][2]) + (*ma).val[0][2] * ((*ma).val[1][0] * (*ma).val[2][1] - (*ma).val[1][1] * (*ma).val[2][0]);

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            (*ans).val[i][j] = ta.val[i][j] / det;
}
