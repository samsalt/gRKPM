#pragma once
#include <iostream>
#include <array>
#include <string.h>
#include <exodusII.h>
#include "control.h"
#include "gxcEcho.h"
#include "gxcMath.h"
// #include "rkpmCudaMethod.h"
#include "gType.h"
#include <cuda_runtime.h>
#include "gxcNodeSet.h"
#include <chrono>

class grkpm
{
    public:
    //method
    void readControl();
    void preprocess();
    void stateInitiate();
    void solve();
    void hexVolumePosition(const std::array<std::array<double, 3>, 8> &xyzel, double &volume, std::array<double, 3> &position, double &win);

    void neighborSearch();
    //data
    control simulationParameter;

    cellDsp* hostDsp;
    cellPosition* hostPosition;
    cellForce* hostForce;

    cellDsp* dspDev;
    cellPosition* positionDev;
    cellForce* forceDev;
    gmNodeNeighbor* nodeNeighbor;
    gmShape* shape;
    gmShapeGradient* shapeGradient;

    cudaError_t err = cudaSuccess;

    domainBound modelBound {};
    gmbinInfo binInfo {};
    double normalWin {1.3};
    bool shapeFlag {false};

    int threadsPerBlock {64};
    int blocksPerGrid {};
    int nc {};

    int** hostEssentialNode;
    int** essentialNodeDev;
    std::vector<int> nodeNumInSet;


    std::vector<gxcEssentialNodeSet> essentialNodeSet;

    // destructor
    ~grkpm()
    {
        free(hostDsp);
        free(hostPosition);
        free(hostForce);
        err = cudaFree(dspDev);
        err = cudaFree(forceDev);
        err = cudaFree(positionDev);
        err = cudaFree(nodeNeighbor);
        err = cudaFree(shape);
        err = cudaFree(shapeGradient);
        for (int i=0; i<simulationParameter.sideSetNum;i++)
        {
            free(hostEssentialNode[i]);
            cudaFree(essentialNodeDev[i]);
        }
    }

};