#pragma once
#include "gType.h"

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h> 



__global__ void findBin(const cellPosition* positionDev,int *nodeBinId,int nc, const gmbinInfo binInfo, const domainBound modelBound);

__global__ void cleanBin(int nb, gmBinNode* binNode);

__global__ void countBin(const int *nodeBinId, int nc, gmBinNode* binNode,const gmbinInfo binInfo);

__global__ void nodeNeighborSearch(const gmBinNode *binNode, const int nc, gmNodeNeighbor *nodeNeighbor, const gmbinInfo binInfo, int *nodeBinId, const cellPosition *positionDev);

__global__ void updateRK(int nc, gmNodeNeighbor *nodeNeighbor, const cellPosition *positionDev, gmShape *shape, gmShapeGradient *shapeGradient);

__global__ void predictor(int nc, double dlt, cellDsp *dspDev);

__global__ void fintCal(int nc, cellDsp *dspDev, gmNodeNeighbor* nodeNeighbor, gmShapeGradient* shapeGradient, cellForce* forceDev);

__global__ void corrector(int nc, double dlt, cellDsp *dspDev);

__device__ void getPhi(vec4* hxPtr, vec4* phiPtr, double win);

__device__ int findThisId(gmNodeNeighbor* nodeNeighbor, int neiborId, int index);

__global__ void assemble (int nc,cellDsp* dspDev,gmNodeNeighbor* nodeNeighbor,cellForce* forceDev);

__global__ void esstialBoundaryEnforce (int nodeNum, int* essentialNodeList, double currentEssentialBoundaryCondition[3], cellDsp* dspDev);
// __device__ void cleanMat4(mat4* mat);

__device__ void vec4Expand(vec4 *va, vec4 *vb, mat4 *ans);
__device__ void vec4Dot(vec4 *va, vec4 *vb, double *ans);
__device__ void mat4Scale(const double scaler, mat4 *ans);
__device__ void mat4Clear(mat4 *ans);
__device__ void mat4Increase(mat4 *ans, mat4 *increment);
__device__ void mat4Prod(mat4 *ma, mat4 *mb, mat4 *ans);
__device__ void mat4Invert(mat4 *ma, mat4 *ans);
