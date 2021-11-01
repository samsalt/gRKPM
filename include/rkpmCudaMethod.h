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
    double mass {};
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
struct vec4
{
    double val[4] {};
};
struct mat4
{
    double val[4][4] {};
};

__global__ void findBin(const cellPosition* positionDev,int *nodeBinId,int nc, const gmbinInfo binInfo, const domainBound modelBound);

__global__ void cleanBin(int nb, gmBinNode* binNode);

__global__ void countBin(const int *nodeBinId, int nc, gmBinNode* binNode,const gmbinInfo binInfo);

__global__ void nodeNeighborSearch(const gmBinNode *binNode, const int nc, gmNodeNeighbor *nodeNeighbor, const gmbinInfo binInfo, int *nodeBinId, const cellPosition *positionDev);

__global__ void updateRK(int nc, gmNodeNeighbor *nodeNeighbor, const cellPosition *positionDev, gmShape *shape, gmShapeGradient *shapeGradient);

__device__ void getPhi(vec4* hxPtr, vec4* phiPtr, double win);

// __device__ void cleanMat4(mat4* mat);

__device__ void vec4Expand(vec4 *va, vec4 *vb, mat4 *ans);
__device__ void vec4Dot(vec4 *va, vec4 *vb, double *ans);
__device__ void mat4Scale(const double scaler, mat4 *ans);
__device__ void mat4Clear(mat4 *ans);
__device__ void mat4Increase(mat4 *ans, mat4 *increment);
__device__ void mat4Prod(mat4 *ma, mat4 *mb, mat4 *ans);
__device__ void mat4Invert(mat4 *ma, mat4 *ans);
