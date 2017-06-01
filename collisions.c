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
#include "collisions-helpers.h"


/* MPI signatures

int MPI_Alltoallv(const void *sendbuf, const int sendcounts[],
    const int sdispls[], MPI_Datatype sendtype,
    void *recvbuf, const int recvcounts[],
    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)

int MPI_Allgatherv(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
    const int displs[], MPI_Datatype recvtype, MPI_Comm comm)

*/


void distributeConfiguration(int myRank, int * numStars, float (*initVelocities)[2], float * mass) {
  float config[8];
  if (myRank == 0) {
    config[0] = (float) numStars[0];
    config[1] = (float) numStars[1];
    config[2] = initVelocities[0][0];
    config[3] = initVelocities[0][1];
    config[4] = initVelocities[1][0];
    config[5] = initVelocities[1][1];
    config[6] = mass[0];
    config[7] = mass[1];
  }

  MPI_Bcast(config, 8, MPI_FLOAT, 0, MPI_COMM_WORLD);

  if (myRank != 0) {
    numStars[0] = round(config[0]);
    numStars[1] = round(config[1]);
    initVelocities[0][0] = config[2];
    initVelocities[0][1] = config[3];
    initVelocities[1][0] = config[4];
    initVelocities[1][1] = config[5];
    mass[0] = config[6];
    mass[1] = config[7];
  }
}

/*
 * Process 0:
 * - reads all parameters,
 * - broadcasts all configuration numbers,
 * - reads all stars into variable 'myStars'.
 */
void readInput(int numProcesses, int myRank, char ** filenameGal,
                                     int * numStars, float (*initVelocities)[2], float * mass, nstars_info_t * myStars) {
  FILE * fp[2];
  int galaxy;

  if (myRank != 0) {
    distributeConfiguration(myRank, numStars, initVelocities, mass);
    myStars[0] = initStars(0, 0, true, true);
    myStars[0] = initStars(0, 1, true, true);

  } else { // process 0 reads from input
    for (galaxy = 0; galaxy < 2; galaxy++) {
      fp[galaxy] = fopen(filenameGal[galaxy], "r");
      if (fp[galaxy] == NULL) {
          fprintf(stderr, "ERROR: file \"%s\" could not be opened!\n", filenameGal[galaxy]);
          exit(1);
      }
      fscanf(fp[galaxy], "%d", &numStars[galaxy]);
      fscanf(fp[galaxy], "%f", &initVelocities[galaxy][0]);
      fscanf(fp[galaxy], "%f", &initVelocities[galaxy][1]);
      fscanf(fp[galaxy], "%f", &mass[galaxy]);
    }

    distributeConfiguration(myRank, numStars, initVelocities, mass);
    // TODO: asynchronous

    // read stars into myStars
    for (galaxy = 0; galaxy < 2; galaxy++) {
      myStars[galaxy] = initStars(numStars[galaxy], galaxy, true, true);
      for (int i = 0; i < numStars[galaxy]; i++) {
        for (int dim = 0; dim < 2; dim++) {
          fscanf(fp[galaxy], "%f", &myStars[galaxy].starsPositions[dim][i]);
        }
      }
      fclose(fp[galaxy]);
    }
  }
}

void printStars(nstars_info_t stars) {
  int n = stars.n;
  printf("STARS:\n");
  for (int i = 0; i < n; i++) {
    printf(" %f %f\n", stars.starsPositions[0][i], stars.starsPositions[1][i]);
  }
}

void checkConfig(int myRank, int * numStars, float (*initVelocities)[2], float * masses, nstars_info_t * myStars) {
  if (myRank > 1)
    return;

  #ifdef DEBUG
  for (int gal = 0; gal < 2; gal++) {
    printf("P[%d]: numStars[%d], v[%f, %f], m[%f]\n",
      myRank, numStars[gal], initVelocities[gal][0], initVelocities[gal][1], masses[gal]);
  }
  #endif
  printStars(myStars[0]);
  printStars(myStars[1]);
}

void computeWorldSize(int numProcesses, int myRank, int * gridSize, float * minPosition, float * maxPosition,
                      float * worldSize, float * blockSize, int * myGridId, nstars_info_t * myStars) {
  // calculate min and max position
  int gal, dim;
  float min, max, diffHalf;
  float minMax[4];

  // only process 0 computes min and max
  if (myRank == 0) {
    for (dim = 0; dim < 2; dim++) {
      // TODO: czy istnieje >= 1 gwiazda
      min = myStars[0].starsPositions[dim][0];
      max = myStars[0].starsPositions[dim][0];

      // find min and max in both galaxies
      for (gal = 0; gal < 2; gal++) {
        countMinMax(myStars[gal].starsPositions[dim], myStars[gal].n, &min, &max);
      }
      #ifdef DEBUG
      printf("MinMax: %f %f \n", min, max);
      #endif
      diffHalf = (max - min) / 2;
      minMax[dim] = min - diffHalf;
      minMax[2 + dim] = max + diffHalf;
      // these are expanded universe boundaries
    }
  }

  // all processes receive min and max
  MPI_Bcast(minMax, 4, MPI_FLOAT, 0, MPI_COMM_WORLD);

  for (dim = 0; dim < 2; dim++) {
    min = minMax[dim];
    max = minMax[2 + dim];
    minPosition[dim] = min;
    maxPosition[dim] = max;
    worldSize[dim] = max - min;
    blockSize[dim] = (max - min) / gridSize[dim];
  }
  rankToGridId(myRank, myGridId, gridSize[0]);
  // 0 1 2 --> (0,0) (1,0) (2,0)
  // 3 4 5 --> (0,1) (1,1) (2,1)

  #ifdef DEBUG
  // TODO remove?
  assert(myRank == gridIdToRank(myGridId[0], myGridId[1], gridSize[0]));

  for (int dim = 0; dim < 2; dim++) {
    printf("$WORLD: P[%d][%d, %d]: dim: %d world from [%f] to [%f], size: [%f], blockSize: [%f]\n",
      myRank, myGridId[0], myGridId[1], dim, minPosition[dim],  maxPosition[dim], worldSize[dim], blockSize[dim]);
  }
  #endif
}

void exchangeCountData(int numProcesses, int * countOutData, int * countInData, int * sum) {
  MPI_Alltoall(countOutData, numProcesses, MPI_INT,
               countInData, numProcesses, MPI_INT, MPI_COMM_WORLD);
  // TODO in place?
  int lsum = 0;
  for (int i = 0; i < numProcesses; i++) {
    lsum += countInData[i];
  }
  *sum = lsum;
}

nstars_info_t exchangeStars(int numProcesses, int myRank, nstars_info_t myStars, float * minPosition, float * blockSize, int gridSizeX) {
  nstars_info_t myNewStars;
  int sumInData;

  int * countOutData = malloc(numProcesses * sizeof(int));
  int * countInData = malloc(numProcesses * sizeof(int));

  sortStars(numProcesses, &myStars, countOutData, minPosition, blockSize, gridSizeX);
  exchangeCountData(numProcesses, countOutData, countInData, &sumInData);

  myNewStars = initStars(sumInData, myStars.galaxy, myStars.withAccelerations, myStars.withAllInfo);

  // TODO!
  MPI_Alltoallv(myStars.starsPositions[0], countOutData, countOutData, MPI_FLOAT,
                myNewStars.starsPositions[0], countInData, countInData, MPI_FLOAT, MPI_COMM_WORLD);

  freeStars(myStars);
  free(countInData);
  free(countOutData);

  return myNewStars;
}

void gatherStars(int numProcesses, const nstars_info_t myStars, nstars_info_t allStars) {
  // TODO
  int n = myStars.n;
  int * countInData = malloc(numProcesses * sizeof(int));

  MPI_Allgather(&n, 1, MPI_INT,
                countInData, 1, MPI_INT, MPI_COMM_WORLD);

  // TODO!
  MPI_Allgatherv(myStars.starsPositions[0], n, MPI_FLOAT,
                 allStars.starsPositions[0], countInData, countInData, MPI_FLOAT, MPI_COMM_WORLD);
}

void computeNewPositions(nstars_info_t stars, float dt) {
  int n = stars.n;
  for (int i = 0; i < n; i++) {
    stars.starsPositions[0][i] += stars.starsVelocities[0][i] * dt + stars.starsAccelerations[0][i] * dt * dt / 2;
  }
}

void computeNewAccelerationsAndVelocities(nstars_info_t myStars, nstars_info_t allStars, float dt, const float otherMass) {
  int n = myStars.n;
  int allN = allStars.n;
  int i,j;
  float x, y;
  float distX, distY;
  float dist3;
  float newA[2];
  const float GM = G * otherMass;

  for (i = 0; i < n; i++) {
    newA[0] = 0;
    newA[1] = 0;
    x = myStars.starsPositions[0][i];
    y = myStars.starsPositions[1][i];

    for (j = 0; j < allN; j++) {
      distX = allStars.starsPositions[0][j] - x;
      distY = allStars.starsPositions[1][j] - y;
      dist3 = distX * distX + distY * distY;
      dist3 = dist3 * sqrt(dist3);  // r^3

      newA[0] += distX / dist3;
      newA[1] += distY / dist3;
    }
    newA[0] *= GM;
    newA[1] *= GM;

    myStars.starsVelocities[0][i] += (myStars.starsAccelerations[0][i] + newA[0]) * dt / 2;
    myStars.starsVelocities[1][i] += (myStars.starsAccelerations[1][i] + newA[1]) * dt / 2;
    myStars.starsAccelerations[0][i] = newA[0];
    myStars.starsAccelerations[1][i] = newA[1];
    // TODO what if F > FLOAT_MAX / 2
    // TODO what to do with me?
  }
  // TODO subtract force from me
  // TODO
}

void outputPositions(int numProcesses, int myRank, nstars_info_t myStars, int galaxy, int iter) {
  // TODO
}

void outputFinalPositions(int numProcesses, int myRank, nstars_info_t myStars, int galaxy) {
  // TODO
}


int main(int argc, char * argv[]) {
  // stars
  int galaxy;
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

  galaxy = 0;

  readInput(numProcesses, myRank, filenameGal, numStars, initVelocities, masses, myStars);
  // now everyone has filled and initialized: numStars, initVelocities, masses, myStars
  // process 0 has all stars

  #ifdef DEBUG
  checkConfig(myRank, numStars, initVelocities, masses, myStars);
  // TODO remove
  #endif

  computeWorldSize(numProcesses, myRank, gridSize, minPosition, maxPosition, worldSize, blockSize, myGridId, myStars);
  allStars[galaxy] = initStars(numStars[galaxy], galaxy, false, false);


  // now all processes know all parameters
  // we are ready for computation


  /************************* COMPUTATION ********************************/

  // HEAD

  // TODO: iteration 0
  // allStars are ready in process 0
  // need to compute only accelerations, not velocities
  // if verbose: outputPositions

  iterNum = (int) (maxSimulationTime / timeStep);
  for (iter = 0; iter < iterNum; iter++) {
    myStars[galaxy] = exchangeStars(numProcesses, myRank, myStars[galaxy], minPosition, blockSize, gridSize[0]);

    computeNewPositions(myStars[galaxy], timeStep);

    if (iter < iterNum - 1) {  // in last iteration further computation is not necessary
      gatherStars(numProcesses, myStars[galaxy], allStars[galaxy]);
      // TODO other galaxy
      computeNewAccelerationsAndVelocities(myStars[galaxy], allStars[galaxy], timeStep, masses[galaxy]);
      if (verbose) {
        outputPositions(numProcesses, myRank, myStars[galaxy], galaxy, iter);
      }
    }
  }

  outputFinalPositions(numProcesses, myRank, myStars[galaxy], galaxy);

  freeStars(myStars[galaxy]);
  freeStars(allStars[galaxy]);
  // TODO free memory


  MPI_Finalize();
  return 0;
}
