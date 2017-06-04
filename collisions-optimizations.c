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
// TODO remove

#include <math.h>
#include <stdbool.h>
#include "collisions-common.h"
#include "collisions-optimizations.h"

#ifdef OPTIMIZATION_1
/********************************** Collisions-1 *********************************************/
void gatherStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars, int * _dummy1, int * _dummy2) {
  // TODO
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

#endif



#ifdef OPTIMIZATION_2
/********************************** Collisions-2 *********************************************/

// TODO
#define PRINT_MSG_TAG 543
#define MPI_BACK_MESSAGE_TAG 1
#define MPI_FRONT_MESSAGE_TAG 2

#define NEIGHBOURS_CNT 8

static int initializeNeighbours(int * neighRanks, int * myGridId, int * gridSize) {
  int xId = myGridId[0];
  int yId = myGridId[1];
  int maxX = gridSize[0] - 1;
  int maxY = gridSize[1] - 1;
  int gridWidth = gridSize[0];

  int neighCnt = 0;

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

  printf("initializeNeighbours on proc (%d, %d): neighCount: %d!\n", myGridId[0], myGridId[1], neighCnt);
  return neighCnt;
}

void prepareOutCount(int numProcesses, int starsCnt, int * countOutData, int * myGridId, int * gridSize) {
  static int neighCnt = -1;
  static int neighRanks[NEIGHBOURS_CNT];
  if (neighCnt == -1) {
    neighCnt = initializeNeighbours(neighRanks, myGridId, gridSize);
  }

  memset(countOutData, 0, numProcesses * sizeof(int));   // initialize to 0
  // send data to neigbours
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

#endif



#ifdef OPTIMIZATION_3
/********************************** Collisions-3 *********************************************/

void gatherStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars, int * myGridId, int * gridSize) {}
#endif
