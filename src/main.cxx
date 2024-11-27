#include "rkpmProlblem.h"

int main()
{
    grkpm problem;
    problem.readControl();
    problem.preprocess();
    problem.stateInitiate();
    problem.solve();
}