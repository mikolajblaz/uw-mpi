#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <assert.h>

#include <stdbool.h>
#include "collisions-helpers.h"

#define OPTION_VERBOSE "--verbose"
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

void printUsage(char * filename) {
  // TODO
}

int parseArguments(int argc, char * argv[], char ** filenameGal1, char ** filenameGal2) {
  // TODO
  int ret = 0;
  bool verbose;
  if (argc < 2) {
      fprintf(stderr, "ERROR: Too few arguments!\n");
      ret = 1;
  }
  else {
    int argIdx = 1;
    if (argc == 3) {
        if (strncmp(argv[argIdx], OPTION_VERBOSE, strlen(OPTION_VERBOSE)) != 0) {
            fprintf(stderr, "ERROR: Unexpected option '%s'!\n", argv[argIdx]);
            ret = 3;
        }
        verbose = true;
        ++argIdx;
    }
    // numPointsPerDimension = atoi(argv[argIdx]);
  }
  if (ret != 0) {
    printUsage(argv[0]);
    MPI_Finalize();
  }
  return ret;
}

nstars_info_t initStars(int n, int galaxy, bool onlyPositions) {
  nstars_info_t ret;
  ret.n = n;
  ret.galaxy = galaxy;
  ret.starsPositions[0] = malloc(n * sizeof(float));
  ret.starsPositions[1] = malloc(n * sizeof(float));
  ret.onlyPositions = onlyPositions;
  if (!onlyPositions) {
    for (int galaxy = 0; galaxy < 2; galaxy++) {
      ret.starsVelocities[galaxy] = malloc(n * sizeof(float));
      ret.starsAccelerations[galaxy] = malloc(n * sizeof(float));
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

// TODO remove
// void initializeMpiStarType(MPI_Datatype * datatype) {
//   MPI_Datatype type[3] = { MPI_INT, MPI_BOOL, MPI_FLOAT };
//   int blocklen[3] = { 1, 1, 4 };
//   MPI_Aint disp[3];
//   star_t starExample;
//   disp[0] = &starExample.index - &starExample;
//   disp[1] = &starExample.isInFirstGalaxy - &starExample;
//   disp[2] = &starExample.positionX - &starExample;
//
//   MPI_Type_create_struct(3, blocklen, disp, type, datatype);
//   MPI_Type_commit(datatype);
// }
