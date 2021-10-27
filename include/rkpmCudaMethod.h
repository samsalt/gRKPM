#pragma once
#define MAXNODEINBIN 100
#define MAXNODENEIGHBOR 100

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h> 

struct cellDsp
{
    double vel[3] {};
    double acl[3] {};
    double dsp[3] {};
};
struct cellForce
{
    double strain[6] {};
    double stress[6] {};
    double volume {};
};
struct cellPosition
{
    double coo[3] {};
    double win {};
};

struct domainBound
{
    double min[3] {1e10, 1e10, 1e10};
    double max[3] {-1e10, -1e10, -1e10};
};

struct gmbinInfo
{
    double size {};
    int num[3] {};
};

struct gmBinNode
{
    int nodeNum {};
    int nodeId [MAXNODEINBIN];
};

struct gmNodeNeighbor
{
    int neighborNum {};
    int neighborId[MAXNODENEIGHBOR];
};
struct gmShape
{
    double val[MAXNODENEIGHBOR];
};
struct gmShapeGradient
{
    double val[3][MAXNODENEIGHBOR];
};

__global__ void findBin(const cellPosition* positionDev,int *nodeBinId,int nc, const gmbinInfo binInfo, const domainBound modelBound);

__global__ void cleanBin(int nb, gmBinNode* binNode);

__global__ void countBin(const int *nodeBinId, int nc, gmBinNode* binNode,const gmbinInfo binInfo);

__global__ void nodeNeighborSearch(const gmBinNode *binNode, const int nc, gmNodeNeighbor *nodeNeighbor, const gmbinInfo binInfo, int *nodeBinId, const cellPosition *positionDev);