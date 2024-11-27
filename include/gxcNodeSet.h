#pragma once
#include <vector>
#include "control.h"
#include "gxcMath.h"

class gxcEssentialNodeSet
{
public:
    // std::vector<int> nodeId;
    int nodeSetId{};
    std::vector<std::vector<double>> essentialBoundaryConditionSet;
    bool EBCSwitch[3] {};
    double currentEssentialBoundaryCondition[3] {};
    double currentVelocity[3] {};

    void updateCurrentBoundaryCondition(const control &simulationParameter);
    // void Predictor(std::vector<gxcNode> &Node, const double dlt);
    // void Corrector(std::vector<gxcNode> &Node, const double dlt);
};