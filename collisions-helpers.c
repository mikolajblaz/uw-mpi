#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <assert.h>

#include <math.h>
#include <stdbool.h>
#include "collisions-helpers.h"

#define OPT_VERBOSE "-v"
#define OPT_HOR "--hor"
#define OPT_VER "--ver"
#define OPT_GAL1 "--gal1"
#define OPT_GAL2 "--gal2"
#define OPT_DELTA "--delta"
#define OPT_TOTAL "--total"
#define OPT_MAX 8
// etc.

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

void quicksort(int *A, int len) {
  if (len < 2) return;

  int pivot = A[len / 2];

  int i, j;
  for (i = 0, j = len - 1; ; i++, j--)
  {
    while (A[i] < pivot) i++;
    while (A[j] > pivot) j--;

    if (i >= j) break;

    int temp = A[i];
    A[i]     = A[j];
    A[j]     = temp;
  }

  quicksort(A, i);
  quicksort(A + i, len - i);
}

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

/******************** STARS SPECIFIC FUNCTIONS ******************************/

void initializeMpiStarType(MPI_Datatype * datatype) {
  MPI_Datatype type[2] = { MPI_FLOAT, MPI_INT };
  int blocklen[2] = { 8, 1 };
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
  // HEAD
  //memset(countOutData, 0, numProcesses * sizeof(int));
  // TODO memeset

  for (int i = 0; i < numProcesses; i++) {
    countOutData[i] = 0;
  }

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

void writeStarsToFile(nstars_info_t stars, char * filename) {
  int n = stars.n;
  FILE * fp = fopen(filename, "w");
  if (fp == NULL) {
      fprintf(stderr, "ERROR: file \"%s\" could not be created!\n", filename);
      MPI_Finalize();
      exit(1);
  }

  for (int i = 0; i < n; i++) {
    fprintf(fp, "%.1f %.1f\n", stars.stars[i].position[0], stars.stars[i].position[1]);
  }

  fclose(fp);
}




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
