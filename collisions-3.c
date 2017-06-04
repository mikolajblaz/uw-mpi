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

static inline int myStarsCount(int numProcesses, int myRank, int numAllStars) {
  int numBase = numAllStars / numProcesses;
  int numInc = numAllStars % numProcesses;
  return (myRank >= numInc ? numBase : (numBase + 1));
}

void gatherAllStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars, int * countInData, int * dispIn) {
  int n = myStars->n;
  MPI_Allgatherv(myStars->stars, n, MPI_STAR,
                 allStars->stars, countInData, dispIn, MPI_STAR, MPI_COMM_WORLD);
}


void calculateProcessorCounts(int * processorCounts, int numProcesses, int myRank, int numAllStars) {
  for (int i = 0; i < numProcesses; i++) {
    processorCounts[i] = myStarsCount(numProcesses, i, numAllStars);
  }
}

void distributeStars(nstars_info_t * myStars, int numProcesses, int myRank, int * processorCounts, int * processorDisps) {
  int myStarsCnt = processorCounts[myRank];
  int * countOutData = processorCounts;
  int * dispOut = processorDisps;

  star_t * starsTmp = malloc(myStarsCnt * sizeof(star_t));
  FAIL_IF_NULL(starsTmp);

  MPI_Scatterv(myStars->stars, countOutData, dispOut, MPI_STAR,
               starsTmp, myStarsCnt, MPI_STAR, 0, MPI_COMM_WORLD);

  freeStars(myStars);
  myStars->n = myStarsCnt;
  myStars->stars = starsTmp;
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

  // who owns how many stars
  int * (processorCounts[2]);
  int * (processorDisps[2]);
  for (int dim = 0; dim < 2; dim++) {
    processorCounts[dim] = malloc(numProcesses * sizeof(int));
    FAIL_IF_NULL(processorCounts[dim]);
    processorDisps[dim] = malloc(numProcesses * sizeof(int));
    FAIL_IF_NULL(processorDisps[dim]);

    calculateProcessorCounts(processorCounts[dim], numProcesses, myRank, numStars[dim]);
    calculateDisplacements(processorCounts[dim], processorDisps[dim], numProcesses);
    distributeStars(&myStars[dim], numProcesses, myRank, processorCounts[dim], processorDisps[dim]);
  }

  // now all stars are distributed, continue computation like in 1 and 2, except from 'exchangeStars' step

  iterNum = (int) (maxSimulationTime / timeStep);
  for (iter = 0; iter < iterNum; iter++){
    // here we don't exchange stars between processors, just gather them

    gatherAllStars(numProcesses, &myStars[0], &allStars[0], processorCounts[0], processorDisps[0]);
    gatherAllStars(numProcesses, &myStars[1], &allStars[1], processorCounts[1], processorDisps[1]);

    if (verbose) {
      outputPositions(myRank, allStars, 0, iter);
      outputPositions(myRank, allStars, 1, iter);
    }

    computeCoordinates(&myStars[0], allStars, timeStep, masses, (iter == 0 ? initVelocities[0] : NULL), minPosition, maxPosition, worldSize);
    computeCoordinates(&myStars[1], allStars, timeStep, masses, (iter == 0 ? initVelocities[1] : NULL), minPosition, maxPosition, worldSize);
  }

  gatherAllStars(numProcesses, &myStars[0], &allStars[0], processorCounts[0], processorDisps[0]);
  gatherAllStars(numProcesses, &myStars[1], &allStars[1], processorCounts[1], processorDisps[1]);
  outputFinalPositions(myRank, allStars);

  for (int dim = 0; dim < 2; dim++) {
    freeStars(&myStars[dim]);
    freeStars(&allStars[dim]);
    free(processorCounts[dim]);
    free(processorDisps[dim]);
  }

  MPI_Finalize();
  return 0;
}
