#include "collisions_helpers.h"

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


int parseArguments(int argc, char const * argv[], char * filenameGal1, char * filenameGal2) {
  // TODO
  int ret = 0;
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
        verbose = 1;
        ++argIdx;
    }
    numPointsPerDimension = atoi(argv[argIdx]);
  }
  if (ret != 0) {
    printUsage(argv[0]);
    MPI_Finalize();
  }
  return ret;
}

void initializeMpiStarType(MPI_Datatype * datatype) {
  MPI_Datatype type[3] = { MPI_INT, MPI_BOOL, MPI_FLOAT };
  int blocklen[3] = { 1, 1, 4 };
  MPI_Aint disp[3];
  star_t starExample;
  disp[0] = &starExample.index - &starExample;
  disp[1] = &starExample.isInFirstGalaxy - &starExample;
  disp[2] = &starExample.positionX - &starExample;

  MPI_Type_create_struct(3, blocklen, disp, type, datatype);
  MPI_Type_commit(datatype);
}
