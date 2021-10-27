#include "rkpmProlblem.h"

void grkpm::hexVolumePosition(const std::array<std::array<double, 3>, 8> &xyzel, double &volume, std::array<double, 3> &position, double &win)
{
  std::array<std::array<double, 3>, 3> at{};
  std::array<std::array<double, 3>, 3> bt{};
  std::array<std::array<double, 3>, 3> ct{};

  for (int i = 0; i < 3; i++)
  {
    at[0][i] = xyzel[6][i] - xyzel[0][i];
    bt[0][i] = xyzel[6][i] - xyzel[0][i];
    ct[0][i] = xyzel[6][i] - xyzel[0][i];
  }
  for (int i = 0; i < 3; i++)
  {
    at[1][i] = xyzel[1][i] - xyzel[0][i];
    bt[1][i] = xyzel[4][i] - xyzel[0][i];
    ct[1][i] = xyzel[3][i] - xyzel[0][i];
  }
  for (int i = 0; i < 3; i++)
  {
    at[2][i] = xyzel[2][i] - xyzel[5][i];
    bt[2][i] = xyzel[5][i] - xyzel[7][i];
    ct[2][i] = xyzel[7][i] - xyzel[2][i];
  }

  volume = (determinant_3a(at) + determinant_3a(bt) + determinant_3a(ct)) / 6;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 8; j++)
    {
      position[i] += xyzel[j][i];
    }
    position[i] = position[i] * 0.125;
  }
  
  win=0;
  double  DiamaterTemp{};
    for (int i = 0; i < 8; i++)
    {
        for (int j = i + 1; j < 8; j++)
        {
            for (gxcId kk = 0; kk < 3; kk++)
            {
                DiamaterTemp = std::abs(xyzel[i][kk] - xyzel[j][kk]);
                if (DiamaterTemp > win)
                    win = DiamaterTemp;
            }
        }
    }
    win = win*normalWin;
}
void grkpm::preprocess()
{
  int *ids;

  // open and create input and output files
  int CPU_word_size, IO_word_size, error, InputExoid, OutputExoid;
  IO_word_size = sizeof(double);
  CPU_word_size = sizeof(double);

  float version;
  InputExoid = ex_open("etest.exo", EX_READ, &CPU_word_size, &IO_word_size, &version);
  OutputExoid = ex_create("gxcOut.exo", EX_CLOBBER, &CPU_word_size, &IO_word_size);

  // read and write model information
  char title[MAX_LINE_LENGTH + 1];
  int DimensionNum;
  error = ex_get_init(InputExoid, title, &DimensionNum, &SimulationParameter.np, &nc, &SimulationParameter.blockNum, &SimulationParameter.nodeSetNum, &SimulationParameter.sideSetNum);
  error = ex_put_init(OutputExoid, "doublebase", 3, SimulationParameter.np, nc, SimulationParameter.blockNum, 0, 0);


  hostPosition = (cellPosition*)malloc(nc*sizeof(cellPosition));
  hostForce = (cellForce*)malloc(nc*sizeof(cellForce));
  hostDsp = (cellDsp*)malloc(nc*sizeof(cellDsp));

  EchoVar("Number of Blocks", SimulationParameter.blockNum);
  EchoVar("Number of Points", SimulationParameter.np);
  EchoVar("Number of Cells", nc);

  // read and write node coordinates
  double *xcoo, *ycoo, *zcoo;
  xcoo = (double *)calloc(SimulationParameter.np, sizeof(double));
  ycoo = (double *)calloc(SimulationParameter.np, sizeof(double));
  zcoo = (double *)calloc(SimulationParameter.np, sizeof(double));
  error = ex_get_coord(InputExoid, xcoo, ycoo, zcoo);
  error = ex_put_coord(OutputExoid, xcoo, ycoo, zcoo);

  char **coord_names;
  coord_names = (char **)malloc(3);
  coord_names[0] = (char *)malloc(2);
  coord_names[1] = (char *)malloc(2);
  coord_names[2] = (char *)malloc(2);
  std::string tempStr{"X"};
  strcpy(coord_names[0], tempStr.c_str());
  tempStr = "Y";
  strcpy(coord_names[1], tempStr.c_str());
  tempStr = "Z";
  strcpy(coord_names[2], tempStr.c_str());
  error = ex_put_coord_names(OutputExoid, coord_names);
  int cellIdGlobal {};

  // read all element blocks information and new cells in each block
  for (int BlockId = 0; BlockId < SimulationParameter.blockNum; BlockId++)
  {
    int localCellNum{}, vertexPerCellNum{}, AttributeNum{};
    char elem_type[MAX_STR_LENGTH + 1];

    error = ex_get_elem_block(InputExoid, BlockId + 1, elem_type, &localCellNum, &vertexPerCellNum, &AttributeNum);
    error = ex_put_elem_block(OutputExoid, BlockId + 1, elem_type, localCellNum, vertexPerCellNum, AttributeNum);

    if (localCellNum>SimulationParameter.cellNumMax) SimulationParameter.cellNumMax=localCellNum;


    int *connect;
    connect = (int *)calloc(vertexPerCellNum * localCellNum, sizeof(int));
    error = ex_get_elem_conn(InputExoid, BlockId + 1, connect);

    for (int i = 0; i < localCellNum; i++)
    {
      std::array<double, 3> localPosition{};
      std::array<std::array<double, 3>, 8> xyzel;
      for (int j = 0; j < 8; j++)
      {
        xyzel[j][0] = xcoo[connect[i * 8 + j] - 1];
        xyzel[j][1] = ycoo[connect[i * 8 + j] - 1];
        xyzel[j][2] = zcoo[connect[i * 8 + j] - 1];
      }
      double localVolume{}, localWin {};

      hexVolumePosition(xyzel,localVolume,localPosition,localWin);

      hostForce[cellIdGlobal].volume = localVolume;
      hostPosition[cellIdGlobal].win=localWin;
      if (SimulationParameter.winMax<localWin) SimulationParameter.winMax=localWin;
      for (int j = 0; j < 3; j++)
      {
        hostPosition[cellIdGlobal].coo[j] = localPosition[j];
        if (localPosition[j]>modelBound.max[j]) modelBound.max[j]=localPosition[j];
        if (localPosition[j]<modelBound.min[j]) modelBound.min[j]=localPosition[j];
      }


      cellIdGlobal++;
    }
    // error = ex_put_elem_conn(OutputExoid, BlockId + 1, connect); !!!!!!!!!!!!!!!!!!!!!!!! needs to be taken care
    free(connect);
  }

  free(xcoo);
  free(ycoo);
  free(zcoo);
};