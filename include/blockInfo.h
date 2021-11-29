#pragma once

#include "gType.h"



class gBlockInfo
{
    public:
    // data
    double cmat [6][6];
    double materialProperty[20] {};
    double velocityInitial[3]  {};
    int materialType {};
    // method

    void formCmat();

};