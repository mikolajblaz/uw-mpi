/*
  Autor: Mikołaj Błaż
  Nr indeksu: 346862
*/

#ifndef _COLLISIONS_HELPERS_
#define _COLLISIONS_HELPERS_

#include <float.h>

#define G (155893.597)
#define FILENAME_LENGTH 20
#define MAX_NUM (FLT_MAX / 2)

#define FAIL_IF_NULL(ptr) if ((ptr) == NULL) \
  { fprintf(stderr, "ERROR [%s, line %d]: couldn't allocate memory!\n", __FILE__, __LINE__); MPI_Finalize(); exit(1); }

// all information about stars velocities, acceleration and indices
typedef struct star_s {
    float position[2];
    float velocity[2];
    float acceleration[2];
    int index;    // TODO: float?
} star_t;

typedef struct nstars_info_s {
  int n;
  int galaxy;
  star_t * stars;
} nstars_info_t;

extern MPI_Datatype MPI_STAR;

// inline functions
extern inline void rankToGridId(int myRank, int * myGridId, int gridSizeX);
extern inline int gridIdToRank(int myGridIdX, int myGridIdY, int gridSizeX);
// make world a torus (in terms of positions) along 1 (some) dimension
extern inline float cyclePosition(float pos, float minPos, float maxPos, float worldSize);
// who owns star inside the world along 1 (some) dimension in grid ids
extern inline int whoOwnsStarInGridId(float pos, float minPos, float blockSize);
// who (rank) owns star
extern inline int whoOwnsStarInRank(float * position, float * minPosition, float * blockSize, int gridSizeX);


/************************* API ********************************/
void printUsage(char * progName);
int parseArguments(int argc, char * argv[], int * gridSize, char ** filenameGal,
                   float * timeStep, float * maxSimulationTime, bool * verbose);
void initializeMpiStarType(MPI_Datatype * datatype);
nstars_info_t initStars(int n, int galaxy);
void freeStars(nstars_info_t * stars);
void sortStars(int numProcesses, nstars_info_t * stars, int * countOutData, float * minPosition, float * blockSize, int gridSizeX, int myRank);

void distributeConfiguration(int myRank, int * numStars, float (*initVelocities)[2], float * mass);
void readInput(int numProcesses, int myRank, char ** filenameGal,
                                     int * numStars, float (*initVelocities)[2], float * mass, nstars_info_t * myStars);
void printStars(nstars_info_t * stars, int myRank, int iter);
void checkConfig(int myRank, int * numStars, float (*initVelocities)[2], float * masses, nstars_info_t * myStars);
void computeWorldSize(int numProcesses, int myRank, int * gridSize, float * minPosition, float * maxPosition,
                      float * worldSize, float * blockSize, int * myGridId, nstars_info_t * myStars);
void calculateDisplacements(int * counts, int * disp, const int size);
void exchangeStars(int numProcesses, int myRank, nstars_info_t * myStars, float * minPosition, float * blockSize, int gridSizeX);
void computeCoordinates(nstars_info_t * myStars, nstars_info_t * allStars, const float dt, float * masses, float * initVelocities,
                        float * minPosition, float * maxPosition, float * worldSize);
void outputPositions(int myRank, nstars_info_t * myStars, int galaxy, int iter);
void outputFinalPositions(int myRank, nstars_info_t * myStars);

#endif /* _COLLISIONS_HELPERS_ */
