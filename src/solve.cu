#include "rkpmProlblem.h"

void grkpm::neighborSearch()
{
    binInfo.size=SimulationParameter.winMax*1.5;
    for (int i =0; i<3; i++) binInfo.num[i]=floor((modelBound.max[i]-modelBound.min[i])/binInfo.size)+1;

    int *nodeBinId;
    err = cudaMalloc(&nodeBinId,  3*nc * sizeof(int));
    except(err,"Fail to allocate device memory");
    findBin<<<blocksPerGrid,threadsPerBlock>>>(positionDev, nodeBinId, nc, binInfo, modelBound);

    // int *hnodeBinId=(int*)malloc(3*nc * sizeof(int));
    // err = cudaMemcpy(hnodeBinId, nodeBinId, 3*nc * sizeof(int), cudaMemcpyDeviceToHost);
    // except(err,"Fail to tansfer data to device, force");
    // for (int i =0; i<nc; i++)
    // {
    //     EchoVarDebug(hnodeBinId[3*i]);
    //     EchoVarDebug(hnodeBinId[3*i+1]);
    //     EchoVarDebug(hnodeBinId[3*i+2]);
    // }

    int totalBinNum=(binInfo.num[0])*(binInfo.num[1])*(binInfo.num[2]);
    gmBinNode* binNode;
    err = cudaMalloc(&binNode, totalBinNum * sizeof(gmBinNode));
    except(err,"Fail to allocate device memory");

    blocksPerGrid = (totalBinNum + threadsPerBlock - 1) / threadsPerBlock;
    cleanBin<<<blocksPerGrid,threadsPerBlock>>>(totalBinNum, binNode);
    
    blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;
    countBin<<<1,1>>>(nodeBinId,nc,binNode,binInfo);

    // gmBinNode* hbn=(gmBinNode*)malloc(totalBinNum * sizeof(gmBinNode));
    // err = cudaMemcpy(hbn, binNode,totalBinNum * sizeof(gmBinNode), cudaMemcpyDeviceToHost);
    // except(err,"Fail to tansfer data to device, force");
    // for (int i =0; i<totalBinNum; i++)
    // {
    //     EchoVar("bin id", i);
    //     for (int j =0; j<hbn[i].nodeNum; j++) std::cout<<hbn[i].nodeId[j]<<", ";
    //     std::cout<<std::endl;
    // }
    nodeNeighborSearch<<<blocksPerGrid,threadsPerBlock>>>(binNode,nc,nodeNeighbor, binInfo, nodeBinId, positionDev);
    // gmNodeNeighbor* hnn=(gmNodeNeighbor*)malloc(nc * sizeof(gmNodeNeighbor));
    // err = cudaMemcpy(hnn, nodeNeighbor,nc * sizeof(gmNodeNeighbor), cudaMemcpyDeviceToHost);
    // except(err,"Fail to tansfer data to device, force");    
    // for (int i =0; i<nc; i++)
    // {
    //     EchoVar("neighborNum", i);
    //     for (int j =0; j<hnn[i].neighborNum; j++) std::cout<<hnn[i].neighborId[j]<<", ";
    //     std::cout<<std::endl;
    // }
    err = cudaFree(nodeBinId);
    err = cudaFree(binNode);
}

void grkpm::solve()
{

    blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;

    EchoVar("blocksPerGrid",blocksPerGrid);

    err = cudaMalloc(&dspDev, nc * sizeof(cellDsp));
    except(err,"Fail to allocate device memory");
    err = cudaMalloc(&positionDev, nc * sizeof(cellPosition));
    except(err,"Fail to allocate device memory");
    err = cudaMalloc(&forceDev, nc * sizeof(cellForce));
    except(err,"Fail to allocate device memory");
    err = cudaMalloc(&nodeNeighbor, nc * sizeof(gmNodeNeighbor));
    except(err,"Fail to allocate device memory");

    err = cudaMemcpy(dspDev, hostDsp, nc * sizeof(cellDsp), cudaMemcpyHostToDevice);
    except(err,"Fail to tansfer data to device, dsp");
    err = cudaMemcpy(positionDev, hostPosition, nc * sizeof(cellPosition), cudaMemcpyHostToDevice);
    except(err,"Fail to tansfer data to device, position");
    err = cudaMemcpy(forceDev, hostForce, nc * sizeof(cellForce), cudaMemcpyHostToDevice);
    except(err,"Fail to tansfer data to device, force");

    neighborSearch();

    blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;

    if (shapeFlag)
    {
        cudaFree(shape);
        cudaFree(shapeGradient);
    }
    err = cudaMalloc(&shape, nc * sizeof(gmShape));
    except(err,"Fail to allocate device memory");
    err = cudaMalloc(&shapeGradient, nc * sizeof(gmShapeGradient));
    except(err,"Fail to allocate device memory");
    updateRK<<<blocksPerGrid,threadsPerBlock>>>(nc, nodeNeighbor, positionDev, shape, shapeGradient);

    gmShape* hshape=(gmShape*)malloc(nc*sizeof(gmShape));
    err = cudaMemcpy(hshape, shape,nc*sizeof(gmShape), cudaMemcpyDeviceToHost);
    except(err,"Fail to tansfer data to device, force");
    for (int i =0; i<nc; i++)
    {
        EchoVar("node id", i);
        for (int j =0; j<14; j++) std::cout<<hshape[i].val[j]<<", ";
        std::cout<<std::endl;
    }
}