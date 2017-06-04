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
#include "collisions-optimizations.h"

#ifdef OPTIMIZATION_1
/********************************** Collisions-1 *********************************************/
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

#endif

#ifdef OPTIMIZATION_2
/********************************** Collisions-1 *********************************************/
void gatherStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars) {}

#endif

#ifdef OPTIMIZATION_3
/********************************** Collisions-1 *********************************************/

void gatherStars(int numProcesses, nstars_info_t * myStars, nstars_info_t * allStars) {}
#endif
