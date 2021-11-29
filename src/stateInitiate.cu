#include "rkpmProlblem.h"

void grkpm::stateInitiate()
{
    int i=0;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);


    blocksPerGrid = (nc + threadsPerBlock - 1) / threadsPerBlock;

    EchoVar("blocksPerGrid",blocksPerGrid);

    for (int i =0; i< simulationParameter.blockNum; i++)
    {
           simulationParameter.blockInfo[i].formCmat();
    }

    err = cudaMalloc(&dspDev, nc * sizeof(cellDsp));
    except(err,"Fail to allocate device memory, dspDev");
    err = cudaMalloc(&positionDev, nc * sizeof(cellPosition));
    except(err,"Fail to allocate device memory, position");
    err = cudaMalloc(&forceDev, nc * sizeof(cellForce));
    except(err,"Fail to allocate device memory, force");
    err = cudaMalloc(&nodeNeighbor, nc * sizeof(gmNodeNeighbor));
    except(err,"Fail to allocate device memory, nodeNeighbor");

    err = cudaMemcpy(dspDev, hostDsp, nc * sizeof(cellDsp), cudaMemcpyHostToDevice);
    except(err,"Fail to tansfer data to device, dsp");
    err = cudaMemcpy(positionDev, hostPosition, nc * sizeof(cellPosition), cudaMemcpyHostToDevice);
    except(err,"Fail to tansfer data to device, position");
    err = cudaMemcpy(forceDev, hostForce, nc * sizeof(cellForce), cudaMemcpyHostToDevice);
    except(err,"Fail to tansfer data to device, force");

    err = cudaMalloc(&(essentialNodeDev), simulationParameter.sideSetNum * sizeof(int*));
    except(err,"Fail to allocate device memory, essential");

    for (int i=0; i<simulationParameter.sideSetNum;i++)
    {
        err = cudaMalloc(&(essentialNodeDev[i]), nodeNumInSet[i] * sizeof(int));
        except(err,"Fail to allocate device memory, essential node set");
        err = cudaMemcpy(essentialNodeDev[i], hostEssentialNode[i], nodeNumInSet[i]  * sizeof(int), cudaMemcpyHostToDevice);
        except(err,"Fail to tansfer data to device, essential");
        free(hostEssentialNode[i]);
    }
    free(hostEssentialNode);

}
