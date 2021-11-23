#pragma once
#define MAXNODEINBIN 100
#define MAXNODENEIGHBOR 100

struct vec4
{
    double val[4] {};
};
struct mat4
{
    double val[4][4] {};
};
struct mat6
{
    double val[6][6] {};
};

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
    mat6 cmat {};
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
