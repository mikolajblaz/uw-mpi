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

#define OPT_VERBOSE "-v"
#define OPT_HOR "--hor"
#define OPT_VER "--ver"
#define OPT_GAL1 "--gal1"
#define OPT_GAL2 "--gal2"
#define OPT_DELTA "--delta"
#define OPT_TOTAL "--total"
#define OPT_MAX 8


/************************* HELPERS ********************************/

// compute min and max, with min and max already initialized!
void countMinMax(star_t * stars, const int dim, const int size, float * min, float * max) {
  float lmin, lmax, curr;
  lmin = *min;
  lmax = *max;

  for (int i = 0; i < size; i++) {
    curr = stars[i].position[dim];
    if (curr < lmin)
      lmin = curr;
    if (curr > lmax)
      lmax = curr;
  }

  *min = lmin;
  *max = lmax;
}


void printArray(char * name, int myRank, int * A, int size) {
  printf("%s (proc: %d): [ ", name, myRank);
  for (int i = 0; i < size; i++) {
    printf("%d ", A[i]);
  }
  printf("]\n");
}
void printArrayS(char * name, int myRank, star_t * A, int size) {
  printf("%s (proc: %d): [ ", name, myRank);
  for (int i = 0; i < size; i++) {
    printf("(%f,%f) ", A[i].position[0], A[i].position[1]);
  }
  printf("]\n");
}

// perform counting sort to sort stars
void sortStars(int numProcesses, nstars_info_t * stars, int * countOutData, float * minPosition, float * blockSize, int gridSizeX, int myRank) {
  memset(countOutData, 0, numProcesses * sizeof(int));

  const int n = stars->n;
  if (n == 0)
    return;

  int i;
  int who;
  int * whoOwns = malloc(n * sizeof(int));
  FAIL_IF_NULL(whoOwns);
  int * countPrefixSum = malloc(numProcesses * sizeof(int));
  FAIL_IF_NULL(countPrefixSum);

  for (i = 0; i < n; i++) {
    who = whoOwnsStarInRank(
      stars->stars[i].position,
      minPosition,
      blockSize,
      gridSizeX
    );
    #ifdef DEBUG
    assert(0 <= who && who < numProcesses);
    // printf("Who: %d\n", who);
    #endif
    whoOwns[i] = who;
    countOutData[who] += 1;
  }

  //printArray("countOutData", myRank, countOutData, numProcesses);

  // TODO: send could (should!) be here
  // count[i] += count[i-1]

  countPrefixSum[0] = countOutData[0];
  for (i = 1; i < numProcesses; i++) {
    countPrefixSum[i] = countPrefixSum[i - 1] + countOutData[i];
  }
  // output[countPrefixSum[whoOwns[i]] - 1] = whoOwns[i];

  //printArray("countPrefixSum", myRank, countPrefixSum, numProcesses);

  star_t * starsTmp = malloc(n * sizeof(star_t));
  FAIL_IF_NULL(starsTmp);

  int dest;
  for (i = 0; i < n; i++) {
    dest = --countPrefixSum[whoOwns[i]];
    starsTmp[dest] = stars->stars[i];
  }

  free(countPrefixSum);
  free(whoOwns);

  freeStars(stars);
  stars->stars = starsTmp;
}

void sortAllStarsForPrinting(nstars_info_t * stars) {
  int n = stars->n;
  star_t * starsTmp = malloc(n * sizeof(star_t));
  for (int i = 0; i < n; i++) {
    starsTmp[stars->stars[i].index] = stars->stars[i];
  }
  freeStars(stars);
  stars->stars = starsTmp;
}

void writeStarsToFile(nstars_info_t * stars, char * filename) {
  sortAllStarsForPrinting(stars);
  int n = stars->n;
  FILE * fp = fopen(filename, "w");
  if (fp == NULL) {
      fprintf(stderr, "ERROR: file \"%s\" could not be created!\n", filename);
      MPI_Finalize();
      exit(1);
  }

  #ifdef DEBUG
  for (int i = 0; i < n; i++) {
    fprintf(fp, "%d %.1f %.1f\n", stars->stars[i].index, stars->stars[i].position[0], stars->stars[i].position[1]);
  }
  #endif
  for (int i = 0; i < n; i++) {
    fprintf(fp, "%.1f %.1f\n", stars->stars[i].position[0], stars->stars[i].position[1]);
  }

  fclose(fp);
}


/************************* INLINES ********************************/
// inline functions declared in collisions-helpers.h

void rankToGridId(int myRank, int * myGridId, int gridSizeX) {
  // 3 4 5 --> (0,1) (1,1) (2,1)
  // 0 1 2 --> (0,0) (1,0) (2,0)
  myGridId[0] = myRank % gridSizeX;
  myGridId[1] = myRank / gridSizeX;
}

int gridIdToRank(int myGridIdX, int myGridIdY, int gridSizeX) {
  return myGridIdY * gridSizeX + myGridIdX;
}

// make world a torus (in terms of positions) along 1 (some) dimension
float cyclePosition(float pos, float minPos, float maxPos, float worldSize) {
  // works for (pos < minPos) too!
  return (pos > maxPos || pos < minPos) ? (pos - floor((pos - minPos) / worldSize) * worldSize) : pos;
}

// who owns star inside the world along 1 (some) dimension in grid ids
int whoOwnsStarInGridId(float pos, float minPos, float blockSize) {
  return (pos - minPos) / blockSize;
}

// who (rank) owns star
int whoOwnsStarInRank(float * position, float * minPosition, float * blockSize, int gridSizeX) {
  return gridIdToRank(
    whoOwnsStarInGridId(position[0], minPosition[0], blockSize[0]),
    whoOwnsStarInGridId(position[1], minPosition[1], blockSize[1]),
    gridSizeX
  );
}


/************************* API ********************************/

void printUsage(char * progName) {
  fprintf(stderr, "Usage:\n"
    "    %s [-v] --hor <H> --ver <V> --gal1 <G1> --gal2 <G2> --delta <D> --total <T>\n"
    "    (arguments order is not relevant)\n",
    progName);
}

int parseArguments(int argc, char * argv[], int * gridSize, char ** filenameGal,
                   float * timeStep, float * maxSimulationTime, bool * verbose) {
  int ret = 0;
  * verbose = false;
  bool argumentsAreSet[6];
  size_t optind;

  if (argc < 13) {
    fprintf(stderr, "ERROR: Too few arguments!\n");
    ret = 1;
  } else if (argc > 14) {
    fprintf(stderr, "ERROR: Too many arguments!\n");
    ret = 2;
  } else {
    for (optind = 1; optind < argc && argv[optind][0] == '-'; optind++) {
      if (strncmp(argv[optind], OPT_VERBOSE, OPT_MAX) == 0) {
        * verbose = true;
      } else if (optind + 1 >= argc) {
        fprintf(stderr, "ERROR: No argument for option '%s'!\n", argv[optind]);
        ret = 3;
      } else {  // another argument is waiting - OK
        if (strncmp(argv[optind], OPT_HOR, OPT_MAX) == 0) {
          argumentsAreSet[0] = true;
          gridSize[0] = atof(argv[++optind]);
        } else if (strncmp(argv[optind], OPT_VER, OPT_MAX) == 0) {
          argumentsAreSet[1] = true;
          gridSize[1] = atof(argv[++optind]);
        } else if (strncmp(argv[optind], OPT_GAL1, OPT_MAX) == 0) {
          argumentsAreSet[2] = true;
          filenameGal[0] = argv[++optind];
        } else if (strncmp(argv[optind], OPT_GAL2, OPT_MAX) == 0) {
          argumentsAreSet[3] = true;
          filenameGal[1] = argv[++optind];
        } else if (strncmp(argv[optind], OPT_DELTA, OPT_MAX) == 0) {
          argumentsAreSet[4] = true;
          *timeStep = atof(argv[++optind]);
        } else if (strncmp(argv[optind], OPT_TOTAL, OPT_MAX) == 0) {
          argumentsAreSet[5] = true;
          *maxSimulationTime = atof(argv[++optind]);
        } else {
          fprintf(stderr, "ERROR: Unexpected option '%s'!\n", argv[optind]);
          ret = 4;
        }
      }
    }
  }

  if (ret != 0) {
    printUsage(argv[0]);
    return ret;
  }

  // check correctness
  for (int i = 0; i < 6; i++) {
    if (argumentsAreSet[i] == false) {
      fprintf(stderr, "ERROR: some required arguments not passed!\n");
      printUsage(argv[0]);
      return 5;
    }
  }

  if (gridSize[0] <= 0 || gridSize[1] <= 0) {
    fprintf(stderr, "ERROR: grid sizes must be positive!\n");
    printUsage(argv[0]);
    return 6;
  }

  #ifdef DEBUG
  // TODO remove
  printf("#%d #%d #%s #%s #%f #%f\n", gridSize[0], gridSize[1], filenameGal[0], filenameGal[1], *timeStep, *maxSimulationTime);
  #endif

  return 0;
}

void initializeMpiStarType(MPI_Datatype * datatype) {
  MPI_Datatype type[2] = { MPI_FLOAT, MPI_INT };
  int blocklen[2] = { 6, 1 };
  MPI_Aint disp[2];
  star_t starExample;
  disp[0] = (void*) &starExample.position[0] - (void*) &starExample;
  disp[1] = (void*) &starExample.index - (void*) &starExample;
  #ifdef DEBUG
  printf("DISPLACEMENT: %ld %ld \n", disp[0], disp[1]);
  #endif

  MPI_Type_create_struct(2, blocklen, disp, type, datatype);
  MPI_Type_commit(datatype);
}

nstars_info_t initStars(int n, int galaxy) {
  nstars_info_t ret;
  ret.n = n;
  ret.galaxy = galaxy;
  if (n == 0) {
    ret.stars = NULL;
  } else {
    ret.stars = malloc(n * sizeof(star_t));
    FAIL_IF_NULL(ret.stars);
  }
  return ret;
}

void freeStars(nstars_info_t * stars) {
  free(stars->stars);
}

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
        myStars[galaxy].stars[i].index = i;
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

  MPI_Alltoall(countOutData, 1, MPI_INT,
               countInData, 1, MPI_INT, MPI_COMM_WORLD);
  // TODO in place?

  calculateDisplacements(countOutData, dispOut, numProcesses);
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

  freeStars(myStars);
  myStars->n = sumInData;
  myStars->stars = starsTmp;
}

void computeCoordinates(nstars_info_t * myStars, nstars_info_t * allStars, const float dt, float * masses, float * initVelocities,
                        float * minPosition, float * maxPosition, float * worldSize) {
  int n = myStars->n;
  int allN;
  int i,j;
  int gal;
  const int myGal = myStars->galaxy;
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
        // TODO: optimize
        if (myGal == gal && index == allStars[gal].stars[j].index)  // not with myself
          continue;

        distX = allStars[gal].stars[j].position[0] - x; // minus sign before F is here
        distY = allStars[gal].stars[j].position[1] - y;
        dist3 = distX * distX + distY * distY;
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

    if (newA[0][0] > MAX_NUM || newA[0][1] > MAX_NUM) { // force too big
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

      // these are positions for the next iteration (t + dt)
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
    writeStarsToFile(&myStars[galaxy], filename);
  }
}

void outputFinalPositions(int myRank, nstars_info_t * myStars) {
  if (myRank == 0) {
    writeStarsToFile(&myStars[0], "res1.txt");
    writeStarsToFile(&myStars[1], "res2.txt");
  }
}
