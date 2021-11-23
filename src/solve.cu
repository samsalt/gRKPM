#include "rkpmProlblem.h"
#include "rkpmCudaMethod.h"

void grkpm::neighborSearch()
{
    binInfo.size = simulationParameter.winMax * 1.5;
    for (int i = 0; i < 3; i++)
        binInfo.num[i] = floor((modelBound.max[i] - modelBound.min[i]) / binInfo.size) + 1;

    int *nodeBinId;
    err = cudaMalloc(&nodeBinId, 3 * nc * sizeof(int));
    except(err, "Fail to allocate device memory,  nodeBinId");
    findBin<<<blocksPerGrid, threadsPerBlock>>>(positionDev, nodeBinId, nc, binInfo, modelBound);

    int *hnodeBinId = (int *)malloc(3 * nc * sizeof(int));
    err = cudaMemcpy(hnodeBinId, nodeBinId, 3 * nc * sizeof(int), cudaMemcpyDeviceToHost);
    except(err, "Fail to tansfer data to device, force");
    for (int i = 0; i < nc; i++)
    {
        EchoVarDebug(hnodeBinId[3 * i]);
        EchoVarDebug(hnodeBinId[3 * i + 1]);
        EchoVarDebug(hnodeBinId[3 * i + 2]);
    }

    int totalBinNum = (binInfo.num[0]) * (binInfo.num[1]) * (binInfo.num[2]);
    gmBinNode *binNode;
    err = cudaMalloc(&binNode, totalBinNum * sizeof(gmBinNode));
    except(err, "Fail to allocate device memory");

    blocksPerGrid = (totalBinNum + threadsPerBlock - 1) / threadsPerBlock;
    cleanBin<<<blocksPerGrid, threadsPerBlock>>>(totalBinNum, binNode);

    blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;
    countBin<<<1, 1>>>(nodeBinId, nc, binNode, binInfo);

    // gmBinNode* hbn=(gmBinNode*)malloc(totalBinNum * sizeof(gmBinNode));
    // err = cudaMemcpy(hbn, binNode,totalBinNum * sizeof(gmBinNode), cudaMemcpyDeviceToHost);
    // except(err,"Fail to tansfer data to device, force");
    // for (int i =0; i<totalBinNum; i++)
    // {
    //     EchoVar("bin id", i);
    //     for (int j =0; j<hbn[i].nodeNum; j++) std::cout<<hbn[i].nodeId[j]<<", ";
    //     std::cout<<std::endl;
    // }
    // free(hbn);

    nodeNeighborSearch<<<blocksPerGrid, threadsPerBlock>>>(binNode, nc, nodeNeighbor, binInfo, nodeBinId, positionDev);

    // gmNodeNeighbor* hnn=(gmNodeNeighbor*)malloc(nc * sizeof(gmNodeNeighbor));
    // err = cudaMemcpy(hnn, nodeNeighbor,nc * sizeof(gmNodeNeighbor), cudaMemcpyDeviceToHost);
    // except(err,"Fail to tansfer data to device, force");
    // for (int i =0; i<nc; i++)
    // {
    //     EchoVar("neighborNum", i);
    //     for (int j =0; j<hnn[i].neighborNum; j++) std::cout<<hnn[i].neighborId[j]<<", ";
    //     std::cout<<std::endl;
    // }
    // free(hnn)
    err = cudaFree(nodeBinId);
    except(err, "Fail to free device memory, nodeBinId");
    err = cudaFree(binNode);
    except(err, "Fail to free device memory, binnode");
}

void grkpm::solve()
{
    neighborSearch();

    blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;

    if (shapeFlag)
    {
        cudaFree(shape);
        except(err, "Fail to free device memory, shape");
        cudaFree(shapeGradient);
        except(err, "Fail to free device memory, shapeGradient");
    }
    EchoVar("gmshape", sizeof(gmShape));
    err = cudaMalloc(&shape, nc * sizeof(gmShape));
    except(err, "Fail to allocate device memory, shape");
    err = cudaMalloc(&shapeGradient, nc * sizeof(gmShapeGradient));
    except(err, "Fail to allocate device memory, shapeGradient");
    updateRK<<<blocksPerGrid, threadsPerBlock>>>(nc, nodeNeighbor, positionDev, shape, shapeGradient);
    cudaDeviceSynchronize();
    shapeFlag = true;

    EchoStr("Start time integrating");
    // chrono::steady_clock sc;
    // auto start = sc.now();

    while (simulationParameter.timeCurrent < simulationParameter.timeEnd)
    {
        for (unsigned i = 0; i < essentialNodeSet.size(); i++)
            essentialNodeSet[i].updateCurrentBoundaryCondition(simulationParameter);

        simulationParameter.timeCurrent+=simulationParameter.dlt;
    }
    // auto end = sc.now();
    // auto simulationTime = static_cast<chrono::duration<double>>(end - start);
    // EchoVar("Simulation duration (sec): ", simulationTime.count());

    // gmShape* hshape=(gmShape*)malloc(nc*sizeof(gmShape));
    // gmShapeGradient* hgrashape=(gmShapeGradient*)malloc(nc*sizeof(gmShapeGradient));
    // err = cudaMemcpy(hshape, shape,nc*sizeof(gmShape), cudaMemcpyDeviceToHost);
    // err = cudaMemcpy(hgrashape, shapeGradient,nc*sizeof(gmShapeGradient), cudaMemcpyDeviceToHost);
    // except(err,"Fail to tansfer data to device, force");
    // for (int i =0; i<nc; i++)
    // {
    //     EchoVar("node id", i);
    //     for (int j =0; j<16; j++) std::cout<<hshape[i].val[j]<<", ";
    //     std::cout<<std::endl;
    //     EchoVar("neighbor", i);
    //     for (int j =0; j<14; j++) std::cout<<hgrashape[i].val[0][j]<<", ";
    //     std::cout<<std::endl;
    //     EchoVar("neighbor", i);
    //     for (int j =0; j<14; j++) std::cout<<hgrashape[i].val[1][j]<<", ";
    //     std::cout<<std::endl;
    // }
}