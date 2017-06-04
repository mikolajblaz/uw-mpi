/*
  Autor: Mikołaj Błaż
  Nr indeksu: 346862
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "collisions-common.h"

MPI_Datatype MPI_STAR;

#define NEIGHBOURS_CNT 9

static inline int gridIdToRank(int myGridIdX, int myGridIdY, int gridSizeX) {
  return myGridIdY * gridSizeX + myGridIdX;
}

void gatherAllStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars, int numAllStars) {
  if (allStars->n != numAllStars) {
    freeStars(allStars);
    allStars->n = numAllStars;
    allStars->stars = malloc(numAllStars * sizeof(star_t));
  }
  int n = myStars->n;
  int * countInData = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(countInData);
  int * dispIn = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(dispIn);

  MPI_Allgather(&n, 1, MPI_INT,
                countInData, 1, MPI_INT, MPI_COMM_WORLD);

  calculateDisplacements(countInData, dispIn, numProcesses);

  MPI_Allgatherv(myStars->stars, n, MPI_STAR,
                 allStars->stars, countInData, dispIn, MPI_STAR, MPI_COMM_WORLD);

  free(countInData);
  free(dispIn);
}

static int initializeNeighbours(int * neighRanks, int * myGridId, int * gridSize) {
  int xId = myGridId[0];
  int yId = myGridId[1];
  int maxX = gridSize[0] - 1;
  int maxY = gridSize[1] - 1;
  int gridWidth = gridSize[0];

  int neighCnt = 0;

  neighRanks[neighCnt++] = gridIdToRank(xId, yId, gridWidth); // send to myself
  if (xId > 0) {
    neighRanks[neighCnt++] = gridIdToRank(xId - 1, yId, gridWidth);
    if (yId > 0) {
      neighRanks[neighCnt++] = gridIdToRank(xId - 1, yId - 1, gridWidth);
    }
    if (yId < maxY) {
      neighRanks[neighCnt++] = gridIdToRank(xId - 1, yId + 1, gridWidth);
    }
  }
  if (xId < maxX) {
    neighRanks[neighCnt++] = gridIdToRank(xId + 1, yId, gridWidth);
    if (yId > 0) {
      neighRanks[neighCnt++] = gridIdToRank(xId + 1, yId - 1, gridWidth);
    }
    if (yId < maxY) {
      neighRanks[neighCnt++] = gridIdToRank(xId + 1, yId + 1, gridWidth);
    }
  }
  if (yId > 0) {
    neighRanks[neighCnt++] = gridIdToRank(xId, yId - 1, gridWidth);
  }
  if (yId < maxY) {
    neighRanks[neighCnt++] = gridIdToRank(xId, yId + 1, gridWidth);
  }
  return neighCnt;
}

void prepareOutCount(int numProcesses, int starsCnt, int * countOutData, int * myGridId, int * gridSize) {
  static int neighCnt = -1;
  static int neighRanks[NEIGHBOURS_CNT];
  if (neighCnt == -1) {
    neighCnt = initializeNeighbours(neighRanks, myGridId, gridSize);
  }

  memset(countOutData, 0, numProcesses * sizeof(int));   // initialize to 0
  // send data to neigbours (and myself)
  for (int neigh = 0; neigh < neighCnt; neigh++) {
    countOutData[neighRanks[neigh]] = starsCnt;
  }
}


void gatherStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars, int * myGridId, int * gridSize) {
  int sumInData;

  int * countOutData = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(countOutData);
  int * countInData = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(countInData);
  int * dispOut = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(dispOut);
  int * dispIn = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(dispIn);

  // HEAD

  prepareOutCount(numProcesses, myStars->n, countOutData, myGridId, gridSize);

  MPI_Alltoall(countOutData, 1, MPI_INT,
               countInData, 1, MPI_INT, MPI_COMM_WORLD);
  // TODO in place?

  memset(dispOut, 0, numProcesses * sizeof(int));   // all processes receive the same data from array myStars->stars
  calculateDisplacements(countInData, dispIn, numProcesses);

  sumInData = 0;
  for (int i = 0; i < numProcesses; i++) {
    sumInData += countInData[i];
  }
  // printf("sumInData(%d): %d\n", myRank, sumInData);
  star_t * starsTmp = malloc(sumInData * sizeof(star_t));
  FAIL_IF_NULL(starsTmp);

  // printStars(myStars, myRank, -777);

  // printArrayS("myStarsOld", myRank, myStars.stars, myStars.n);

  MPI_Alltoallv(myStars->stars, countOutData, dispOut, MPI_STAR,
                starsTmp, countInData, dispIn, MPI_STAR, MPI_COMM_WORLD);

  //printArrayS("myStarsNew", myRank, starsTmp, sumInData);

  free(countInData);
  free(countOutData);
  free(dispIn);
  free(dispOut);

  freeStars(allStars);
  allStars->n = sumInData;
  allStars->stars = starsTmp;
}



int main(int argc, char * argv[]) {
  // stars
  int numStars[2];
  nstars_info_t myStars[2];
  nstars_info_t allStars[2];
  float initVelocities[2][2]; // [galaxy][x, y]
  float masses[2];

  // coordinates
  float minPosition[2];
  float maxPosition[2];
  float worldSize[2];
  float blockSize[2];
  int gridSize[2]; // {width, height}
  int myGridId[2];
  int numProcesses; // numProcesses = gridSize[0] * gridSize[1]
  int myRank;

  // computation
  int iter;
  int iterNum;
  float timeStep;
  float maxSimulationTime;

  // other
  int ret;
  bool verbose;
  char * filenameGal[2];

  /************************* INITIALIZATION ********************************/

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  initializeMpiStarType(&MPI_STAR);

  ret = parseArguments(argc, argv, gridSize, filenameGal, &timeStep, &maxSimulationTime, &verbose);
  if (ret == 0 && numProcesses != gridSize[0] * gridSize[1]) {
    fprintf(stderr, "ERROR: Number of active processes must be equal to (%d * %d)!\n", gridSize[0], gridSize[1]);
    ret = 7;
  }
  if (ret != 0) {
    MPI_Finalize();
    return ret;
  }

  // Here arguments are correct

  readInput(numProcesses, myRank, filenameGal, numStars, initVelocities, masses, myStars);
  // now everyone has filled and initialized: numStars, initVelocities, masses, myStars
  // process 0 has all stars

  computeWorldSize(numProcesses, myRank, gridSize, minPosition, maxPosition, worldSize, blockSize, myGridId, myStars);
  allStars[0] = initStars(numStars[0], 0);
  allStars[1] = initStars(numStars[1], 1);

  // now all processes know all parameters
  // we are ready for computation


  /************************* COMPUTATION ********************************/

  iterNum = (int) (maxSimulationTime / timeStep);
  for (iter = 0; iter < iterNum; iter++){
    // TODO: asynch
    exchangeStars(numProcesses, myRank, &myStars[0], minPosition, blockSize, gridSize[0]);
    exchangeStars(numProcesses, myRank, &myStars[1], minPosition, blockSize, gridSize[0]);  // gridSize[0]!

    if (verbose) {
      gatherAllStars(numProcesses, &myStars[0], &allStars[0], numStars[0]);
      gatherAllStars(numProcesses, &myStars[1], &allStars[1], numStars[1]);
      outputPositions(myRank, allStars, 0, iter);
      outputPositions(myRank, allStars, 1, iter);
    }

    // TODO: MPI_Bcast instead of gather stars
    gatherStars(numProcesses, &myStars[0], &allStars[0], myGridId, gridSize);
    gatherStars(numProcesses, &myStars[1], &allStars[1], myGridId, gridSize);

    computeCoordinates(&myStars[0], allStars, timeStep, masses, (iter == 0 ? initVelocities[0] : NULL), minPosition, maxPosition, worldSize);
    computeCoordinates(&myStars[1], allStars, timeStep, masses, (iter == 0 ? initVelocities[1] : NULL), minPosition, maxPosition, worldSize);
  }

  gatherAllStars(numProcesses, &myStars[0], &allStars[0], numStars[0]);
  gatherAllStars(numProcesses, &myStars[1], &allStars[1], numStars[1]);
  outputFinalPositions(myRank, allStars);

  freeStars(&myStars[0]);
  freeStars(&myStars[1]);
  freeStars(&allStars[0]);
  freeStars(&allStars[1]);

  MPI_Finalize();
  return 0;
}
