#include "gxcNodeSet.h"


void gxcEssentialNodeSet::updateCurrentBoundaryCondition(const control &SimulationParameter)
{
    //find the current step
    for (unsigned j = 0; j < essentialBoundaryConditionSet.size(); j++)
    {
        if (SimulationParameter.timeCurrent <= essentialBoundaryConditionSet[j][0])
        {
            double LastStepVal{}, LastTime{};

            for (int k = 0; k < 3; k++)
            {
                if (j != 0)
                {
                    LastStepVal = essentialBoundaryConditionSet[j - 1][k + 1];
                    LastTime = essentialBoundaryConditionSet[j - 1][0];
                }
                currentEssentialBoundaryCondition[k] = LinearInterpolation(SimulationParameter.timeCurrent, LastTime, essentialBoundaryConditionSet[j][0], LastStepVal, essentialBoundaryConditionSet[j][k + 1]);
                currentVelocity[k] = (essentialBoundaryConditionSet[j][k + 1] - LastStepVal) / (essentialBoundaryConditionSet[j][0] - LastTime);
            }
        }
    }
    // EchoVarDebug(CurrentEssentialBoundaryCondition[0]);
}