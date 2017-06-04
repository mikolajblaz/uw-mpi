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

MPI_Datatype MPI_STAR;

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
  int ret;

  if (myRank != 0) {
    distributeConfiguration(myRank, numStars, initVelocities, mass);
    myStars[0] = initStars(0, 0);
    myStars[1] = initStars(0, 1);

  } else { // process 0 reads from input
    for (galaxy = 0; galaxy < 2; galaxy++) {
      fp[galaxy] = fopen(filenameGal[galaxy], "r");
      if (fp[galaxy] == NULL) {
          fprintf(stderr, "ERROR: file \"%s\" could not be opened!\n", filenameGal[galaxy]);
          MPI_Finalize();
          exit(1);
      }
      ret = fscanf(fp[galaxy], "%d %f %f %f", &numStars[galaxy], &initVelocities[galaxy][0], &initVelocities[galaxy][1], &mass[galaxy]);
      if (ret != 4) {
        fprintf(stderr, "ERROR while reading configuration parameters from input!\n");
        MPI_Finalize();
        exit(1);
      }
    }

    distributeConfiguration(myRank, numStars, initVelocities, mass);
    // TODO: asynchronous

    // read stars into myStars
    for (galaxy = 0; galaxy < 2; galaxy++) {
      myStars[galaxy] = initStars(numStars[galaxy], galaxy);
      for (int i = 0; i < numStars[galaxy]; i++) {
        ret = fscanf(fp[galaxy], "%f %f", &myStars[galaxy].stars[i].position[0], &myStars[galaxy].stars[i].position[1]);
        if (ret != 2) {
          fprintf(stderr, "ERROR while reading stars positions (line %d)!\n", i);
          MPI_Finalize();
          exit(1);
        }
      }
      fclose(fp[galaxy]);
    }
  }
}

void printStars(nstars_info_t * stars, int myRank, int iter) {
  int n = stars->n;
  printf("STARS[%d] on proc [%d] in iter [%d]: [", stars->galaxy, myRank, iter);
  for (int i = 0; i < n; i++) {
    printf("(%f,%f) ", stars->stars[i].position[0], stars->stars[i].position[1]);
  }
  printf("]\n");
}

void checkConfig(int myRank, int * numStars, float (*initVelocities)[2], float * masses, nstars_info_t * myStars) {
  if (myRank > 1)
    return;

  for (int gal = 0; gal < 2; gal++) {
    printf("P[%d]: numStars[%d], v[%f, %f], m[%f]\n",
      myRank, numStars[gal], initVelocities[gal][0], initVelocities[gal][1], masses[gal]);
  }

  printStars(&myStars[0], myRank, 0);
  printStars(&myStars[1], myRank, 0);
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
      min = myStars[0].stars[0].position[dim];
      max = myStars[0].stars[0].position[dim];

      // find min and max in both galaxies
      for (gal = 0; gal < 2; gal++) {
        countMinMax(myStars[gal].stars, dim, myStars[gal].n, &min, &max);
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

  // if (myRank == 1) {
  //   for (int dim = 0; dim < 2; dim++) {
  //     printf("$WORLD: P[%d][%d, %d]: dim: %d world from [%f] to [%f], size: [%f], blockSize: [%f]\n",
  //       myRank, myGridId[0], myGridId[1], dim, minPosition[dim],  maxPosition[dim], worldSize[dim], blockSize[dim]);
  //   }
  // }
  #endif
}

void exchangeCountData(int numProcesses, int * countOutData, int * countInData, int * sum) {
  //printArray("countOutData", myRank, countOutData, numProcesses);
  MPI_Alltoall(countOutData, 1, MPI_INT,
               countInData, 1, MPI_INT, MPI_COMM_WORLD);
  // TODO in place?
  //printArray("countInData", myRank, countInData, numProcesses);
  int lsum = 0;
  for (int i = 0; i < numProcesses; i++) {
    lsum += countInData[i];
  }
  *sum = lsum;
}

void calculateDisplacements(int * counts, int * disp, const int size) {
  disp[0] = 0;
	for (int i = 1; i < size; i++){
		disp[i] = counts[i - 1] + disp[i - 1];
	}
}

void exchangeStars(int numProcesses, int myRank, nstars_info_t * myStars, float * minPosition, float * blockSize, int gridSizeX) {
  int sumInData;

  int * countOutData = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(countOutData);
  int * countInData = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(countInData);
  int * dispOut = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(dispOut);
  int * dispIn = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(dispIn);

  sortStars(numProcesses, myStars, countOutData, minPosition, blockSize, gridSizeX, myRank);

  exchangeCountData(numProcesses, countOutData, countInData, &sumInData);
  calculateDisplacements(countOutData, dispOut, numProcesses);
  calculateDisplacements(countInData, dispIn, numProcesses);

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

  freeStars(myStars);
  myStars->n = sumInData;
  myStars->stars = starsTmp;
}

void gatherStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars) {
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


void computeCoordinates(nstars_info_t * myStars, nstars_info_t * allStars, float dt, float * masses, float * initVelocities,
                        float * minPosition, float * maxPosition, float * worldSize) {
  int n = myStars->n;
  int allN;
  int i,j;
  int gal;
  int index;
  float x, y;
  float distX, distY;
  float dist3;
  float newA[2][2];
  const float GM[2] = {G * masses[0], G * masses[1]};
  bool immobilize;

  for (i = 0; i < n; i++) {
    x = myStars->stars[i].position[0];
    y = myStars->stars[i].position[1];
    index = myStars->stars[i].index;
    immobilize = false;

    for (gal = 0; gal < 2 && !immobilize; gal++) {
      allN = allStars[gal].n;
      newA[gal][0] = 0;
      newA[gal][1] = 0;
      for (j = 0; j < allN; j++) {
        if (index == allStars[gal].stars[j].index)  // not with myself
          continue;

        distX = allStars[gal].stars[j].position[0] - x;
        distY = allStars[gal].stars[j].position[1] - y;
        dist3 = distX * distX + distY * distY;
        // TODO: if dist != 0
        if (dist3 == 0) {     // division by zero - danger!
          immobilize = true;
          break;
        }
        dist3 = dist3 * sqrt(dist3);  // r^3

        newA[gal][0] += distX / dist3;
        newA[gal][1] += distY / dist3;
      }
      newA[gal][0] *= GM[gal];    // proper units
      newA[gal][1] *= GM[gal];
    }
    newA[0][0] += newA[1][0];     // sum from both galaxies
    newA[0][1] += newA[1][1];

    if (newA[0][0] > MAX_NUM || newA[0][1] > MAX_NUM) {
      immobilize = true;
    }

    if (!immobilize) {
      myStars->stars[i].acceleration[0] = newA[0][0];
      myStars->stars[i].acceleration[1] = newA[0][1];

      if (initVelocities != NULL) {   // first iteration only
        myStars->stars[i].velocity[0] = initVelocities[0];
        myStars->stars[i].velocity[1] = initVelocities[1];
      } else {
        myStars->stars[i].velocity[0] += (myStars->stars[i].acceleration[0] + newA[0][0]) * dt / 2;
        myStars->stars[i].velocity[1] += (myStars->stars[i].acceleration[1] + newA[0][1]) * dt / 2;
      }

      myStars->stars[i].position[0] += myStars->stars[i].velocity[0] * dt + newA[0][0] * dt * dt / 2;
      myStars->stars[i].position[0] = cyclePosition(myStars->stars[i].position[0], minPosition[0], maxPosition[0], worldSize[0]);
      myStars->stars[i].position[1] += myStars->stars[i].velocity[1] * dt + newA[0][1] * dt * dt / 2;
      myStars->stars[i].position[1] = cyclePosition(myStars->stars[i].position[1], minPosition[1], maxPosition[1], worldSize[1]);
    }
  }
}

void outputPositions(int myRank, nstars_info_t * myStars, int galaxy, int iter) {
  if (myRank == 0) {
    char filename[FILENAME_LENGTH];
    snprintf(filename, FILENAME_LENGTH, "res%d_%d.txt", galaxy + 1, iter);
    writeStarsToFile(myStars[galaxy], filename);
  }
}

void outputFinalPositions(int myRank, nstars_info_t * myStars) {
  if (myRank == 0) {
    writeStarsToFile(myStars[0], "res1.txt");
    writeStarsToFile(myStars[1], "res2.txt");
  }
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

  galaxy = 0;

  readInput(numProcesses, myRank, filenameGal, numStars, initVelocities, masses, myStars);
  // now everyone has filled and initialized: numStars, initVelocities, masses, myStars
  // process 0 has all stars

  #ifdef DEBUG
  // checkConfig(myRank, numStars, initVelocities, masses, myStars);
  // TODO remove
  #endif


  computeWorldSize(numProcesses, myRank, gridSize, minPosition, maxPosition, worldSize, blockSize, myGridId, myStars);
  allStars[0] = initStars(numStars[0], 0);
  allStars[1] = initStars(numStars[1], 1);


  // now all processes know all parameters
  // we are ready for computation


  /************************* COMPUTATION ********************************/

  iter = 0;
  printStars(&myStars[0], myRank, -99);
  printStars(&myStars[1], myRank, -99);
  exchangeStars(numProcesses, myRank, &myStars[0], minPosition, blockSize, gridSize[0]);
  exchangeStars(numProcesses, myRank, &myStars[1], minPosition, blockSize, gridSize[0]);
  printStars(&myStars[0], myRank, -3);
  printStars(&myStars[1], myRank, -3);
  // TODO: MPI_Bcast instead of gather stars
  gatherStars(numProcesses, &myStars[0], &allStars[0]);
  gatherStars(numProcesses, &myStars[1], &allStars[1]);
  printStars(&allStars[0], myRank, -20);
  printStars(&allStars[1], myRank, -20);
  computeCoordinates(&myStars[0], allStars, timeStep, masses, initVelocities[0], minPosition, maxPosition, worldSize);
  computeCoordinates(&myStars[1], allStars, timeStep, masses, initVelocities[1], minPosition, maxPosition, worldSize);
  printStars(&myStars[0], myRank, -1);
  printStars(&myStars[1], myRank, -1);

  iterNum = (int) (maxSimulationTime / timeStep);
  for (iter = 1; iter < iterNum; iter++) {
    if (verbose) {
      outputPositions(myRank, myStars, 0, iter - 1);
      outputPositions(myRank, myStars, 1, iter - 1);
      printStars(&myStars[0], myRank, iter - 1);
      printStars(&myStars[1], myRank, iter - 1);
    }

    exchangeStars(numProcesses, myRank, &myStars[0], minPosition, blockSize, gridSize[0]);
    exchangeStars(numProcesses, myRank, &myStars[1], minPosition, blockSize, gridSize[0]);

    gatherStars(numProcesses, &myStars[0], &allStars[0]);
    gatherStars(numProcesses, &myStars[1], &allStars[1]);

    computeCoordinates(&myStars[0], allStars, timeStep, masses, NULL, minPosition, maxPosition, worldSize);
    computeCoordinates(&myStars[1], allStars, timeStep, masses, NULL, minPosition, maxPosition, worldSize);
  }

  outputFinalPositions(myRank, myStars);

  freeStars(&myStars[0]);
  freeStars(&myStars[1]);
  freeStars(&allStars[0]);
  freeStars(&allStars[1]);
  // TODO free memory


  MPI_Finalize();
  return 0;
}
