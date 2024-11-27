#ifndef gxc_TYPE_H
#define gxc_TYPE_H

#include <cinttypes>

#ifdef ID_RANGE
#if ID_RANGE == 64
typedef unsigned long long gxcId;
#elif ID_RANGE == 32
typedef unsigned long gxcId;
#endif
#else
typedef int gxcId;
#endif

#ifdef DATA_PRECISION
#if DATA_PRECISION == 64
typedef double gxcData;
#elif DATA_PRECISION == 32
typedef float gxcData;
#endif
#else
typedef double gxcData;
#endif

#endif