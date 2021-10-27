#pragma once

#include "gxcType.h"
 #include <array>
#include <vector>
#include <cstddef>
 #include "math.h"

gxcData determinant_3(const std::vector<std::vector<gxcData>> &m);
void InvMat3(const std::vector<std::vector<gxcData>> &m, std::vector<std::vector<gxcData>> &Invm);
gxcData LinearInterpolation(gxcData x,gxcData xa,gxcData xb,gxcData ya,gxcData yb);
void DeviatoricStress(const std::vector<gxcData> Stress,double *DevStress, double &SphereStress);
gxcData StressNorm(const double *Stress);
gxcData VecNorm3(const double *Vec);
std::array<std::array<double, 3>, 3> prodA33(const std::vector<std::vector<double>>& a1, const std::vector<std::vector<double>>& a2);

template <typename T, std::size_t nr, std::size_t nc>
inline std::array<std::array<T, nc>, nr> addA(const std::array<std::array<T, nc>, nr>& a1, const std::array<std::array<T, nc>, nr>& a2)
{
    std::array<std::array<T, nc>, nr> ans;

    for (auto i=0; i<nr;i++)
        for (auto j=0; j<nc;j++)
            ans[i][j]=a1[i][j]+a2[i][j];

    return ans;
}

template <typename T, std::size_t nr, std::size_t nc>
inline std::array<std::array<T, nc>, nr> subA(const std::array<std::array<T, nc>, nr>& a1, const std::array<std::array<T, nc>, nr>& a2)
{
    std::array<std::array<T, nc>, nr> ans;

    for (auto i=0; i<nr;i++)
        for (auto j=0; j<nc;j++)
            ans[i][j]=a1[i][j]-a2[i][j];

    return ans;
}

template <typename T, std::size_t nr, std::size_t nc, std::size_t ni>
inline std::array<std::array<T, nc>, nr> prodA(const std::array<std::array<T, ni>, nr>& a1, const std::array<std::array<T, nc>, ni>& a2)
{
    std::array<std::array<T, nc>, nr> ans {};

    for (auto i=0; i<nr;i++)
        for (auto j=0; j<nc;j++)
            for (auto p=0; p<ni;p++)
                ans[i][j]+=a1[i][p]*a2[p][j];

    return ans;
}


template <typename T, std::size_t nr, std::size_t nc>
inline std::array<std::array<T, nr>, nc> trans(const std::array<std::array<T, nc>, nr>& a1)
{
    std::array<std::array<T, nr>, nc> ans;

    for (auto i=0; i<nr;i++)
        for (auto j=0; j<nc;j++)
            ans[j][i]=a1[i][j];

    return ans;
}

template <typename T>
inline T determinant_3a(const std::array<std::array<T, 3>, 3> &m)
{
  const T a = m[0][0];
  const T b = m[0][1];
  const T c = m[0][2];
  const T d = m[1][0];
  const T e = m[1][1];
  const T f = m[1][2];
  const T g = m[2][0];
  const T h = m[2][1];
  const T k = m[2][2];
  const T determinant = (a * ((e * k) - (f * h))) - (b * ((k * d) - (f * g))) + (c * ((d * h) - (e * g)));
  return determinant;
}

template <typename T>
inline void InvMat3a(const std::array<std::array<T, 3>, 3> &m, std::array<std::array<T, 3>, 3> &Invm)
{
  T invdet;
  invdet = 1 / (determinant_3a(m));
  Invm[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  Invm[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  Invm[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  Invm[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  Invm[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  Invm[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  Invm[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  Invm[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  Invm[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}

template <typename T>
inline void Mat2VoigtStrain(const std::array<std::array<T, 3>, 3> &mat, std::array<std::array<T, 1>, 6> &voigt)
{
  voigt[0][0] = mat[0][0];
  voigt[1][0] = mat[1][1];
  voigt[2][0] = mat[2][2];
  voigt[3][0] = 2 * mat[1][2];
  voigt[4][0] = 2 * mat[0][2];
  voigt[5][0] = 2 * mat[0][1];
}

template <typename T>
inline void Mat2VoigtStress(const std::array<std::array<T, 3>, 3> &mat, std::array<std::array<T, 1>, 6> &voigt)
{
  voigt[0][0] = mat[0][0];
  voigt[1][0] = mat[1][1];
  voigt[2][0] = mat[2][2];
  voigt[3][0] = mat[1][2];
  voigt[4][0] = mat[0][2];
  voigt[5][0] = mat[0][1];
}
template <typename T>
inline void Voigt2MatStress(const std::array<std::array<T, 1>, 6> &voigt, std::array<std::array<T, 3>, 3>  &mat)
{
  mat[0][0] = voigt[0][0];
  mat[1][1] = voigt[1][0];
  mat[2][2] = voigt[2][0];
  mat[1][2] = voigt[3][0];
  mat[0][2] = voigt[4][0];
  mat[1][0] = mat[0][1];
}
template <typename T>
inline void Voigt2MatStrain(const std::array<std::array<T, 1>, 6> &voigt, std::array<std::array<T, 3>, 3>  &mat)
{
  mat[0][0] = voigt[0][0];
  mat[1][1] = voigt[1][0];
  mat[2][2] = voigt[2][0];
  mat[1][2] = 0.5*voigt[3][0];
  mat[0][2] = 0.5*voigt[4][0];
  mat[0][1] = 0.5*voigt[5][0];
  mat[2][1] = mat[1][2];
  mat[2][0] = mat[0][2];
  mat[1][0] = mat[0][1];
}

template <typename T, std::size_t nr, std::size_t nc>
inline std::array<std::array<T, nc>, nr> scaleA(T sc, const std::array<std::array<T, nc>, nr>& a1)
{
    std::array<std::array<T, nc>, nr> ans;

    for (auto i=0; i<nr;i++)
        for (auto j=0; j<nc;j++)
            ans[i][j]=sc*a1[i][j];

    return ans;
}
// const std::vector<std::vector<gxcData>> ShapeHexDxi{
//     {-0.125 * (1.0 - eta) * (1 + mu), -0.125 * (1.0 - eta) * (1.0 - mu), -0.125 * (1.0 + eta) * (1.0 - mu), -0.125 * (1.0 + eta) * (1.0 + mu), 0.125 * (1.0 - eta) * (1.0 + mu), 0.125 * (1.0 - eta) * (1.0 - mu), 0.125 * (1.0 + eta) * (1.0 - mu), 0.125 * (1.0 + eta) * (1.0 + mu)},
//     {-0.125 * (1 - xi) * (1 + mu), -0.125 * (1 - xi) * (1 - mu), 0.125 * (1 - xi) * (1 - mu), 0.125 * (1 - xi) * (1 + mu), -0.125 * (1 + xi) * (1 + mu), -0.125 * (1 + xi) * (1 - mu), 0.125 * (1 + xi) * (1 - mu), 0.125 * (1 + xi) * (1 + mu)},
//     {0.125 * (1 - xi) * (1 - eta), -0.125 * (1 - xi) * (1 - eta), -0.125 * (1 - xi) * (1 + eta), 0.125 * (1 - xi) * (1 + eta), 0.125 * (1 + xi) * (1 - eta), -0.125 * (1 + xi) * (1 - eta), -0.125 * (1 + xi) * (1 + eta), 0.125 * (1 + xi) * (1 + eta)}};
