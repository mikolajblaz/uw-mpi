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

#define OPTION_VERBOSE "--verbose"

#define PRINT_MSG_TAG 543

static int const MPI_BACK_MESSAGE_TAG = 1;
static int const MPI_FRONT_MESSAGE_TAG = 2;

typedef struct star_s
{
    int index;
    bool isInFirstGalaxy;
    float positionX;
    float positionY;
    float velocityX;
    float velocityY;
} star_t;




int parseArguments(int argc, char const * argv[], char * filenameGal1, char * filenameGal2) {
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

void readAndDistributeInput(int numProcesses, int myRank, const char * filenameGal) {
  if (myRank == 0) {
    // read
  }

  if (myRank == 0) {
    // send
  } else {
    // receive
  }

}


int main(int argc, char const * argv[]) {
  // variables describing whole system
  int numProcesses;
  int worldW;          // x dimension
  int worldH;          // y dimension

  // variables describing me
  int myRank;
  int myX;
  int myY;
  int ret;

  char * filenameGal1;
  char * filenameGal2;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  ret = parseArguments(argc, argv);
  if (ret != 0) {
    return ret;
  }

  readAndDistributeInput(numProcesses, filenameGal1);
  readAndDistributeInput(numProcesses, filenameGal2);


  return 0;
}
