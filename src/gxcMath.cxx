#include <gxcMath.h>


gxcData determinant_3(const std::vector<std::vector<gxcData>> &m)
{
  const gxcData a = m[0][0];
  const gxcData b = m[0][1];
  const gxcData c = m[0][2];
  const gxcData d = m[1][0];
  const gxcData e = m[1][1];
  const gxcData f = m[1][2];
  const gxcData g = m[2][0];
  const gxcData h = m[2][1];
  const gxcData k = m[2][2];
  const gxcData determinant = (a * ((e * k) - (f * h))) - (b * ((k * d) - (f * g))) + (c * ((d * h) - (e * g)));
  return determinant;
}

void InvMat3(const std::vector<std::vector<gxcData>> &m, std::vector<std::vector<gxcData>> &Invm)
{
  Invm.resize(3, std::vector<gxcData>(3, 0));
  gxcData invdet;
  invdet = 1 / (determinant_3(m));
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





gxcData LinearInterpolation(gxcData x,gxcData xa,gxcData xb,gxcData ya,gxcData yb)
{
    gxcData xl=xb-xa;
    gxcData yl=yb-ya;
    return ya+(x-xa)*yl/xl;
}

void DeviatoricStress(const std::vector<gxcData> Stress,double *DevStress, double &SphereStress)
{

    SphereStress=(Stress[0]+Stress[1]+Stress[2])/3;
    DevStress[0]=Stress[0]-SphereStress;
    DevStress[1]=Stress[1]-SphereStress;
    DevStress[2]=Stress[2]-SphereStress;
    DevStress[3]=Stress[3];
    DevStress[4]=Stress[4];
    DevStress[5]=Stress[5];		
}

gxcData StressNorm(const double *Stress)
{
  return std::sqrt(Stress[0]*Stress[0]+Stress[1]*Stress[1]+Stress[2]*Stress[2]+2*(Stress[3]*Stress[3]+Stress[4]*Stress[4]+Stress[5]*Stress[5]));
}

gxcData VecNorm3(const double *Vec)
{
  return std::sqrt(Vec[0]*Vec[0]+Vec[1]*Vec[1]+Vec[2]*Vec[2]);
}

std::array<std::array<double, 3>, 3> prodA33(const std::vector<std::vector<double>>& a1, const std::vector<std::vector<double>>& a2)
{
    std::array<std::array<double, 3>, 3> ans {};

    for (auto i=0; i<3;i++)
        for (auto j=0; j<3;j++)
            for (auto p=0; p<a2.size();p++)
                ans[i][j]+=a1[i][p]*a2[p][j];

    return ans;
}