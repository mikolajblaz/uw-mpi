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
#include <stdbool.h>
#include "collisions-helpers.h"


MPI_Datatype MPI_Star;


/* MPI signatures

int MPI_Alltoallv(const void *sendbuf, const int sendcounts[],
    const int sdispls[], MPI_Datatype sendtype,
    void *recvbuf, const int recvcounts[],
    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)

int MPI_Allgatherv(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
    const int displs[], MPI_Datatype recvtype, MPI_Comm comm)

*/

/*
 * Process 0:
 * - reads all parameters,
 * - broadcasts all configuration numbers,
 * - reads all stars into variable 'myStars'.
 */
nstars_info_t readAndDistributeInput(int numProcesses, int myRank, const char * filenameGal, int * numStars, int galaxy) {
  // TODO
  nstars_info_t myStars;
  if (myRank == 0) {
    // read numStars, init velocities, masses
    // MPI_Bcast
    myStars = initStars((*numStars), galaxy, false);
    // read stars into myStars
  } else {
    myStars = initStars(0, galaxy, false);
  }
  return myStars;
}

void exchangeCountData(int numProcesses, int myRank, int * countOutData, int * countInData) {
  // TODO
  // MPI_AllToAll(countInData, countOutData)
}

nstars_info_t exchangeStars(int numProcesses, int myRank, nstars_info_t myStars) {
  // TODO
  nstars_info_t myNewStars;

  int * countOutData = malloc(numProcesses * sizeof(int));
  int * countInData = malloc(numProcesses * sizeof(int));

  sortStars(numProcesses, &myStars, countOutData);
  exchangeCountData(numProcesses, myRank, countOutData, countInData);

  // MPI_AllToAll(countInData, countOutData)

  free(countInData);
  free(countOutData);

  // TODO
  return myNewStars;
}

void gatherStars(int numProcesses, int myRank, const nstars_info_t myStars, nstars_info_t allStars) {
  // TODO
}

void computeNewPositions(nstars_info_t myStars) {
  // TODO
}

void computeNewAccelerationsAndVelocities(nstars_info_t myStars, nstars_info_t allStars) {
  // TODO
}

void outputPositions(int numProcesses, int myRank, nstars_info_t myStars, int galaxy, int iter) {
  // TODO
}

void outputFinalPositions(int numProcesses, int myRank, nstars_info_t myStars, int galaxy) {
  // TODO
}


int main(int argc, char * argv[]) {
  // variables describing whole system
  int numProcesses;
  int worldW;          // x dimension
  int worldH;          // y dimension
  int numStars[2];
  int galaxy;
  nstars_info_t allStars[2];

  // variables describing me
  int myRank;
  int myX;
  int myY;
  nstars_info_t myStars[2];
  nstars_info_t myNewStars[2];

  // computation
  int iter = 0;
  int iterNum;
  float timeStep;
  float maxSimulationTime;

  // other
  int ret;
  int verbose;
  char * filenameGal[2];


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  ret = parseArguments(argc, argv, &filenameGal[0], &filenameGal[1]);
  if (ret != 0) {
    MPI_Finalize();
    return ret;
  }

  galaxy = 0;

  // TODO remove
  // initializeMpiStarType(&MPI_Star);

  myStars[galaxy] = readAndDistributeInput(numProcesses, myRank, filenameGal[galaxy], &numStars[galaxy], galaxy);
  // now process 0 has all stars

  allStars[galaxy] = initStars(numStars[galaxy], galaxy, true);

  // TODO: iteration 0
  // allStars are ready in process 0
  // need to compute only accelerations, not velocities
  // if verbose: outputPositions

  iterNum = (int) (maxSimulationTime / timeStep);
  for (int iter = 0; iter < iterNum; iter++) {
    myNewStars[galaxy] = exchangeStars(numProcesses, myRank, myStars[galaxy]);
    freeStars(myStars[galaxy]);
    myStars[galaxy] = myNewStars[galaxy];

    computeNewPositions(myStars[galaxy]);

    if (iter < iterNum - 1) {  // in last iteration further computation is not necessary
      gatherStars(numProcesses, myRank, myStars[galaxy], allStars[galaxy]);
      computeNewAccelerationsAndVelocities(myStars[galaxy], allStars[galaxy]);
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
