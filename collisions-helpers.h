#ifndef _COLLISIONS_HELPERS_
#define _COLLISIONS_HELPERS_

#define G 155893.597
#define FILENAME_LENGTH 20

// TODO
#define PRINT_MSG_TAG 543
#define MPI_BACK_MESSAGE_TAG 1
#define MPI_FRONT_MESSAGE_TAG 2

#define FAIL_IF_NULL(ptr) if ((ptr) == NULL) \
  { fprintf(stderr, "ERROR [%s, line %d]: couldn't allocate memory!\n", __FILE__, __LINE__); MPI_Finalize(); exit(1); }

// all information about stars velocities, acceleration and indices
typedef struct nstars_info_s {
    int n;
    int galaxy;
    bool withAccelerations; // true if space for accelerations is allocated
    bool withAllInfo;       // true if space for velocities and indices is allocated
    float * starsPositions[2];
    float * starsVelocities[2];
    float * starsAccelerations[2];
    int * indices;
} nstars_info_t;


// helper functions
void countMinMax(float * A, int size, float * min, float * max);
void quicksort(int *A, int len);
int parseArguments(int argc, char * argv[], int * gridSize, char ** filenameGal,
                   float * timeStep, float * maxSimulationTime, bool * verbose);
void writeStarsToFile(nstars_info_t stars, char * filename);

// stars related functions
nstars_info_t initStars(int n, int galaxy, bool withAccelerations, bool withAllInfo);
void freeStars(nstars_info_t stars);
void sortStars(int numProcesses, nstars_info_t * stars, int * countOutData, float * minPosition, float * blockSize, int gridSizeX);



// inline functions
inline void rankToGridId(int myRank, int * myGridId, int gridSizeX) {
  // 0 1 2 --> (0,0) (1,0) (2,0)
  // 3 4 5 --> (0,1) (1,1) (2,1)
  myGridId[0] = myRank % gridSizeX;
  myGridId[1] = myRank / gridSizeX;
}

inline int gridIdToRank(int myGridIdX, int myGridIdY, int gridSizeX) {
  return myGridIdY * gridSizeX + myGridIdX;
}

// make world a torus (in terms of positions) along 1 (some) dimension
inline float cyclePosition(float pos, float minPos, float maxPos, float worldSize) {
  // works for (pos < minPos) too!
  return (pos > maxPos || pos < minPos) ? (pos - floor((pos - minPos) / worldSize) * worldSize) : pos;
}

// who owns star inside the world along 1 (some) dimension in grid ids
inline int whoOwnsStarInGridId(float pos, float minPos, float blockSize) {
  return (pos - minPos) / blockSize;
}

// who (rank) owns star
inline int whoOwnsStarInRank(float positionX, float positionY, float minPositionX, float minPositionY,
                             float blockSizeX, float blockSizeY, int gridSizeX) {
  return gridIdToRank(
    whoOwnsStarInGridId(positionX, minPositionX, blockSizeX),
    whoOwnsStarInGridId(positionY, minPositionY, blockSizeY),
    gridSizeX
  );
}


#endif /* _COLLISIONS_HELPERS_ */
