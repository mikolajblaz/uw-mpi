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


void readAndDistributeInput(int numProcesses, int myRank, const char * filenameGal1, const char * filenameGal2) {
  // TODO
  if (myRank == 0) {
    // read
  }

  if (myRank == 0) {
    // send

  } else {
    myStars = initStars(0);
    // receive
  }

}

int * exchangeCountData(int numProcesses, int myRank, int * countOutData) {
  // TODO
  int * countInData = malloc(numProcesses * sizeof(int));
  // MPI_AllToAll(countInData, countOutData)
  return countInData;
}

void exchangeStars(int numProcesses, int myRank, nstars_info_t myStars) {
  // TODO
  int * countOutData;
  nstars_info_t myNewStars;

  countOutData = sortStars(&myStars);
  exchangeCountData(numProcesses, myRank, countOutData);

  // MPI_AllToAll(countInData, countOutData)

  // TODO
}

void gatherStars(int numProcesses, int myRank, const nstars_info_t myStars, nstars_info_t * allStars) {
  // TODO
}

nstars_info_t initStars(int n, bool onlyPositions) {
  nstars_info_t ret;
  ret.n = n;
  ret.starsPositions[0] = malloc(n * sizeof(float));
  ret.starsPositions[1] = malloc(n * sizeof(float));
  ret.onlyPositions = onlyPositions;
  if (!onlyPositions) {

  }
  return ret;
}

void freeStars(nstars_info_t stars) {
  free(stars.starsPositions[0]);
  free(stars.starsPositions[1]);
}


int main(int argc, char const * argv[]) {
  // variables describing whole system
  int numProcesses;
  int worldW;          // x dimension
  int worldH;          // y dimension
  int numStars;

  int iter = 0;
  int iterNum;
  float timeStep;
  float maxSimulationTime;

  // variables describing me
  int myRank;
  int myX;
  int myY;
  nstars_info_t myStars;
  nstars_info_t myNewStars;
  nstars_info_t allStars;

  // other
  int ret;
  char * filenameGal1;
  char * filenameGal2;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  ret = parseArguments(argc, argv);
  if (ret != 0) {
    MPI_Finalize();
    return ret;
  }

  initializeMpiStarType(&MPI_Star);
  readAndDistributeInput(numProcesses, myRank, filenameGal1, filenameGal2);

  allStars = initStars(numStars);

  iterNum = (int) (maxSimulationTime / timeStep);
  for (int iter = 0; iter < iterNum; iter++) {
    // TODO
    gatherStars(numProcesses, myRank, myStars, &allStars);
    computeForcesAndMove();

    if (iter < iterNum - 1) {  // not in last iteration
      myNewStars = exchangeStars(myStars);
    }
    freeStars(myStars);
    myStars = myNewStars;
  }

  freeStars(allStars);
  // TODO free memory

  MPI_Finalize();
  return 0;
}
