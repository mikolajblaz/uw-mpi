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
extern inline void rankToGridId(int myRank, int * myGridId, int gridSizeX);

extern inline int gridIdToRank(int myGridIdX, int myGridIdY, int gridSizeX);

// make world a torus (in terms of positions) along 1 (some) dimension
extern inline float cyclePosition(float pos, float minPos, float maxPos, float worldSize);

// who owns star inside the world along 1 (some) dimension in grid ids
extern inline int whoOwnsStarInGridId(float pos, float minPos, float blockSize);

// who (rank) owns star
extern inline int whoOwnsStarInRank(float positionX, float positionY, float minPositionX, float minPositionY,
                             float blockSizeX, float blockSizeY, int gridSizeX);


#endif /* _COLLISIONS_HELPERS_ */
