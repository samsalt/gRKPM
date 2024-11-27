#include "rkpmProlblem.h"

void grkpm::readControl()
{

    std::ifstream controlFile("control.dat");
    std::string line;
    int essentialNodeSetNum{};
    // bool HasFileName = 0;
    // std::string &ModelFileName;

    //     if (controlFile)
    //     {
    //         while (controlFile)
    //         {
    //             getline(controlFile, line);
    //             if (line == "# model file name")
    //                 getline(controlFile, ModelFileName);
    //             HasFileName = 1;
    //         }
    //         controlFile.close();
    //     }
    //     else
    //         EchoError("Cannot find control.dat file!");

    //     int sn = simulationParameter.FileName.length();
    //   char CharName[sn + 1];
    //   strcpy(CharName, SimulationParameter.FileName.c_str());

    // open and create input and output files
    int CPU_word_size, IO_word_size, error, InputExoid, OutputExoid;
    IO_word_size = sizeof(double);
    CPU_word_size = sizeof(double);

    float version;
    InputExoid = ex_open("etest.exo", EX_READ, &CPU_word_size, &IO_word_size, &version);
    char title[MAX_LINE_LENGTH + 1];
    int DimensionNum;
    error = ex_get_init(InputExoid, title, &DimensionNum, &simulationParameter.np, &nc, &simulationParameter.blockNum, &simulationParameter.nodeSetNum, &simulationParameter.sideSetNum);
    error = ex_close(InputExoid);
    simulationParameter.blockInfo.resize(simulationParameter.blockNum);

    if (controlFile)
    {
        while (controlFile)
        {
            getline(controlFile, line);
            if (line == "# time duration")
            {
                controlFile >> simulationParameter.timeEnd;
            }
            else if (line == "# time step size")
            {
                controlFile >> simulationParameter.dlt;
            }
            else if (line == "# time output period")
            {
                controlFile >> simulationParameter.timeOutputPeriod;
            }
            else if (line == "# LagrangianType")
            {
                // controlFile >> simulationParameter.LagrangianType;
            }
            else if (line == "# material property")
            {
                int blockId {};
                controlFile >> blockId;
                blockId--;
                
                controlFile >> simulationParameter.blockInfo[blockId].materialType;
                for (int i = 0; i < 20; i++)
                    controlFile >> simulationParameter.blockInfo[blockId].materialProperty[i];
            }
            else if (line == "# block set initial velocity")
            {
                int blockId {};
                controlFile >> blockId;
                blockId--;
                controlFile >> simulationParameter.blockInfo[blockId].velocityInitial[0];
                controlFile >> simulationParameter.blockInfo[blockId].velocityInitial[1];
                controlFile >> simulationParameter.blockInfo[blockId].velocityInitial[2];
            }
            else if (line == "# essential boundary condition")
            {
                gxcId nodeSetId;
                controlFile >> nodeSetId;
                gxcEssentialNodeSet ENTemp;
                essentialNodeSet.push_back(ENTemp);
                essentialNodeSet[essentialNodeSetNum].nodeSetId=nodeSetId-1;
                std::string DisplacementName = "Displacement" + std::to_string(nodeSetId) + ".dat";
                std::ifstream DisplacementFile(DisplacementName);
                std::string lineDisplacement;
                int DisplacementStep{};
                nodeSetId--;
                int swtemp{};
                controlFile >> swtemp;
                if (swtemp)
                    essentialNodeSet[essentialNodeSetNum].EBCSwitch[0] = true;
                controlFile >> swtemp;
                if (swtemp)
                    essentialNodeSet[essentialNodeSetNum].EBCSwitch[1] = true;
                controlFile >> swtemp;
                if (swtemp)
                    essentialNodeSet[essentialNodeSetNum].EBCSwitch[2] = true;

                if (DisplacementFile)
                {
                    std::getline(DisplacementFile, lineDisplacement);
                    DisplacementFile >> DisplacementStep;
                    essentialNodeSet[essentialNodeSetNum].essentialBoundaryConditionSet.resize(DisplacementStep, std::vector<gxcData>(4, 0));
                    std::getline(DisplacementFile, lineDisplacement);
                    std::getline(DisplacementFile, lineDisplacement);
                    for (int i = 0; i < DisplacementStep; i++)
                    {
                        DisplacementFile >> essentialNodeSet[essentialNodeSetNum].essentialBoundaryConditionSet[i][0];
                        DisplacementFile >> essentialNodeSet[essentialNodeSetNum].essentialBoundaryConditionSet[i][1];
                        DisplacementFile >> essentialNodeSet[essentialNodeSetNum].essentialBoundaryConditionSet[i][2];
                        DisplacementFile >> essentialNodeSet[essentialNodeSetNum].essentialBoundaryConditionSet[i][3];
                    }
                    if (essentialNodeSet[essentialNodeSetNum].essentialBoundaryConditionSet[DisplacementStep-1][0]<simulationParameter.timeEnd) EchoError("The time duration in "+DisplacementName+" cannot be smaller than the total simulation time");
                    DisplacementFile.close();
                    essentialNodeSetNum++;
                }
                else
                {
                    EchoError("Cannot find file" + DisplacementName+" for Essential boundary condition");
                }
            }
            // else if (line == "# natural boundary condition")
            // {
            //     gxcId SidesetId;
            //     std::vector<gxcId> LoadingSet;
            //     getline(controlFile, line);
            //     std::stringstream stream(line);
            //     while (1)
            //     {
            //         stream >> SidesetId;
            //         if (!stream)
            //             break;
            //         LoadingSet.push_back(SidesetId-1);
            //     }
            //     LoadingFaceSet.resize(LoadingSet.size());
            //     for (unsigned j = 0; j < LoadingSet.size(); j++)
            //     {
            //         gxcId BlockId, CellId, SetId;
            //         SetId=LoadingSet[j];
            //         BlockId = simulationParameter.BlockCellId[FaceTemp[SetId][0][1]][0];
            //         CellId = simulationParameter.BlockCellId[FaceTemp[SetId][0][1]][1];
            //         LoadingFaceSet[j].VertexNum = Block[BlockId].Cell[CellId]->SideNum;
            //         if (LoadingFaceSet[j].VertexNum == 6)
            //         {
            //             LoadingFaceSet[j].VertexNum = 4;
            //         }
            //         else if (LoadingFaceSet[j].VertexNum == 4)
            //         {
            //             LoadingFaceSet[j].VertexNum = 3;
            //         }
            //         gxcId FaceNodeId[6][4] = {{0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 4, 0}, {0, 1, 2, 3}, {4, 7, 6, 5}};
            //         LoadingFaceSet[SetId].Face.resize(FaceTemp[SetId].size());
            //         for (unsigned jj = 0; jj < FaceTemp[SetId].size(); jj++)
            //         {
            //             BlockId = simulationParameter.BlockCellId[FaceTemp[SetId][jj][1]][0];
            //             CellId = simulationParameter.BlockCellId[FaceTemp[SetId][jj][1]][1];

            //             LoadingFaceSet[SetId].Face[jj].FaceId = FaceTemp[SetId][jj][0];
            //             for (unsigned k = 0; k < 4; k++)
            //                 LoadingFaceSet[SetId].Face[jj].NodeList[k] = Block[BlockId].Cell[CellId]->VertexId[FaceNodeId[FaceTemp[SetId][jj][0]][k]];
            //         }
            //         SidesetId = LoadingSet[j];
            //         std::string LoadingName = "Loading" + std::to_string(SidesetId) + ".dat";
            //         std::ifstream LoadingFile(LoadingName);
            //         std::string lineLoading;
            //         int LoadingStep{};
            //         SidesetId--;
            //         if (LoadingFile)
            //         {
            //             std::getline(LoadingFile, lineLoading);
            //             std::cout << lineLoading;
            //             LoadingFile >> LoadingStep;
            //             LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet.resize(LoadingStep, std::vector<gxcData>(4, 0));
            //             LoadingFaceSet[SidesetId].CurrentNaturalBoundaryConditionSet.resize(3, 0);
            //             std::getline(LoadingFile, lineLoading);
            //             std::getline(LoadingFile, lineLoading);
            //             std::cout << lineLoading;
            //             for (int i = 0; i < LoadingStep; i++)
            //             {
            //                 LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][0];
            //                 LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][1];
            //                 LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][2];
            //                 LoadingFile >> LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[i][3];
            //             }
            //             if (LoadingFaceSet[SidesetId].NaturalBoundaryConditionSet[LoadingStep-1][0]<simulationParameter.TimeEnd) EchoError("Time duration in "+LoadingName+"cannot be smallter than the total simulation time");
            //             LoadingFile.close();
            //         }
            //         else
            //         {
            //             EchoError("Cannot find file" + LoadingName+ " for natural boundary condition");
            //         }
            //     }
            // }
            // else if (line == "# thread number")
            // {
            //     controlFile >> simulationParameter.ThreadNumber;
            // }
            else if (line == "# rk information")
            {
                // simulationParameter.MeshfreeSwitch = true;
                // controlFile >> simulationParameter.MeshfreeInfo.Deg;
                // controlFile >> simulationParameter.MeshfreeInfo.NormalWin;
            }
        }
        // EchoStr("Finish reading control file");
    }
    else
        EchoError("Cannot find control.dat file!");

    controlFile.close();
    // CheckControl();
}