#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <assert.h>

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

  // TODO remove
  printf("#%d #%d #%s #%s #%f #%f\n", gridSize[0], gridSize[1], filenameGal[0], filenameGal[1], *timeStep, *maxSimulationTime);

  return 0;
}

nstars_info_t initStars(int n, int galaxy, bool onlyPositions) {
  nstars_info_t ret;
  ret.n = n;
  ret.galaxy = galaxy;
  ret.starsPositions[0] = malloc(n * sizeof(float));
  ret.starsPositions[1] = malloc(n * sizeof(float));
  ret.onlyPositions = onlyPositions;
  if (!onlyPositions) {
    for (int dim = 0; dim < 2; dim++) {
      ret.starsVelocities[dim] = malloc(n * sizeof(float));
      ret.starsAccelerations[dim] = malloc(n * sizeof(float));
    }
    ret.indices = malloc(n * sizeof(float));
  }
  return ret;
}

void freeStars(nstars_info_t stars) {
  free(stars.starsPositions[0]);
  free(stars.starsPositions[1]);
  if (!stars.onlyPositions) {
    for (int galaxy = 0; galaxy < 2; galaxy++) {
      free(stars.starsVelocities[galaxy]);
      free(stars.starsAccelerations[galaxy]);
    }
    free(stars.indices);
  }
}

void sortStars(int numProcesses, nstars_info_t * stars, int * countOutData) {
  // TODO
}
