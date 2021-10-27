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
                nodeBinId[3*index+dm] = binInfo.num[dm] - 1;
            }
            else if (localCoo < modelBound.min[dm])
            {
                nodeBinId[3*index+dm] = 0;
            }
            else
            {
                nodeBinId[3*index+dm]= __double2int_rd((localCoo - modelBound.min[dm]) / binInfo.size);
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
        i =nodeBinId[3*index];
        j =nodeBinId[3*index+1];
        k =nodeBinId[3*index+2];
        binId=i+j*binInfo.num[0]+k*binInfo.num[0]*binInfo.num[1];
        binNode[binId].nodeNum += 1;
        binNode[binId].nodeId[binNode[binId].nodeNum-1] = index;
    }
}

__device__ void L2Norm(const double  vec[3], double *norm)
{
    *norm= sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

__device__ void hardSearch(const cellPosition *localCellPosition, gmNodeNeighbor *localNodeNeighbor, const cellPosition *positionDev, const gmBinNode *localBinNode)
{
    for (int index = 0; index < (*localBinNode).nodeNum; index++)
    {
        int comId=(*localBinNode).nodeId[index];
        double vec[3] {}, dist {};
        double *dptr=&dist;
        for (int i = 0; i < 3; i++)
            vec[i]=positionDev[comId].coo[i]-(*localCellPosition).coo[i];
        L2Norm(vec,dptr);

        if (dist<positionDev[comId].win)
        {
            (*localNodeNeighbor).neighborId[(*localNodeNeighbor).neighborNum]=comId;
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
        localCellPosition=positionDev[index];

        int XId[3], YId[3], ZId[3];
        for (int i = 0; i < 3; i++)
        {
            XId[i] = nodeBinId[3*index] - 1 + i;
            YId[i] = nodeBinId[3*index + 1] - 1 + i;
            ZId[i] = nodeBinId[3*index + 2] - 1 + i;
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
                                int localBinId=XId[i]+YId[i]*binInfo.num[0]+ZId[i]*binInfo.num[1]*binInfo.num[0];
                                hardSearch(positionDev+index, nodeNeighbor+index,positionDev,binNode+localBinId);
                            }
                        }
                    }
                }
            }
        }
    }
}

