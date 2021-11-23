#pragma once

#include "gType.h"



class gBlockInfo
{
    public:
    // data
    mat6 cmat {};
    double materialProperty[20] {};
    double velocityInitial[3]  {};
    int materialType {};
    // method

    void formCmat();

};