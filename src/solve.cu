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

    nodeNeighborSearch<<<blocksPerGrid, threadsPerBlock>>>(binNode, nc, nodeNeighbor, binInfo, nodeBinId, positionDev);

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

    EchoStr("Start time integration");
    // chrono::steady_clock sc;
    // auto start = sc.now();

    while (simulationParameter.timeCurrent < simulationParameter.timeEnd)
    {
        for (unsigned i = 0; i < essentialNodeSet.size(); i++)
            essentialNodeSet[i].updateCurrentBoundaryCondition(simulationParameter);

        predictor<<<blocksPerGrid, threadsPerBlock>>>(nc, simulationParameter.dlt, dspDev);

        for (int i = 0; i < simulationParameter.sideSetNum; i++)
        {
            blocksPerGrid = (nodeNumInSet[i] + threadsPerBlock - 1) / threadsPerBlock;
            
            esstialBoundaryEnforce<<<blocksPerGrid, threadsPerBlock>>>(nodeNumInSet[i], essentialNodeDev[i], essentialNodeSet[i].currentEssentialBoundaryCondition, dspDev);

            err = cudaMalloc(&(essentialNodeDev[i]), nodeNumInSet[i] * sizeof(int));
            except(err, "Fail to allocate device memory, essential node set");
            err = cudaMemcpy(essentialNodeDev[i], hostEssentialNode[i], nodeNumInSet[i] * sizeof(int), cudaMemcpyHostToDevice);

        }
        blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;

        fintCal<<<blocksPerGrid, threadsPerBlock>>>(nc, dspDev, nodeNeighbor, shapeGradient, forceDev);

        assemble<<<blocksPerGrid, threadsPerBlock>>>(nc, dspDev, nodeNeighbor, forceDev);

        corrector<<<blocksPerGrid, threadsPerBlock>>>(nc, simulationParameter.dlt, dspDev);

        if (simulationParameter.timeCurrent > simulationParameter.timetoOutput)
        {
            simulationParameter.timetoOutput += simulationParameter.timeOutputPeriod;

            err = cudaMemcpy(hostDsp, dspDev, nc * sizeof(cellDsp), cudaMemcpyDeviceToHost);

            except(err, "Fail to tansfer data to device, dsp in output");

            EchoVar("dsp at the end", hostDsp[nc - 1].dsp[0]);
        }

        simulationParameter.timeCurrent += simulationParameter.dlt;
    }
    EchoStr("The simulation is finished!");
    // auto end = sc.now();
    // auto simulationTime = static_cast<chrono::duration<double>>(end - start);
    // EchoVar("Simulation duration (sec): ", simulationTime.count());

}