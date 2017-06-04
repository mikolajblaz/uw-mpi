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

MPI_Datatype MPI_STAR;

void gatherAllStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars) {
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

    // TODO: MPI_Bcast instead of gather stars
    gatherAllStars(numProcesses, &myStars[0], &allStars[0]);
    gatherAllStars(numProcesses, &myStars[1], &allStars[1]);

    if (verbose) {
      outputPositions(myRank, allStars, 0, iter);
      outputPositions(myRank, allStars, 1, iter);
    }

    computeCoordinates(&myStars[0], allStars, timeStep, masses, (iter == 0 ? initVelocities[0] : NULL), minPosition, maxPosition, worldSize);
    computeCoordinates(&myStars[1], allStars, timeStep, masses, (iter == 0 ? initVelocities[1] : NULL), minPosition, maxPosition, worldSize);
  }

  gatherAllStars(numProcesses, &myStars[0], &allStars[0]);
  gatherAllStars(numProcesses, &myStars[1], &allStars[1]);
  // TODO: gather all stars!
  outputFinalPositions(myRank, allStars);

  freeStars(&myStars[0]);
  freeStars(&myStars[1]);
  freeStars(&allStars[0]);
  freeStars(&allStars[1]);

  MPI_Finalize();
  return 0;
}
