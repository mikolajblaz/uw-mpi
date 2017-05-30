#ifndef _COLLISIONS_HELPERS_
#define _COLLISIONS_HELPERS_

// TODO
#define PRINT_MSG_TAG 543
#define MPI_BACK_MESSAGE_TAG 1
#define MPI_FRONT_MESSAGE_TAG 2

// extra information about stars velocities, acceleration and indices
typedef struct nstars_info_s {
    int n;
    int galaxy;
    bool withAccelerations;
    bool withAllInfo;
    float * starsPositions[2];
    float * starsVelocities[2];
    float * starsAccelerations[2];
    int * indices;
} nstars_info_t;


// helper functions
void countMinMax(float * A, int size, float * min, float * max);
void quicksort(int *A, int len);
void sortStars(int numProcesses, nstars_info_t * stars, int * countOutData);
int parseArguments(int argc, char * argv[], int * gridSize, char ** filenameGal,
                   float * timeStep, float * maxSimulationTime, bool * verbose);
nstars_info_t initStars(int n, int galaxy, bool withAccelerations, bool withAllInfo);
void freeStars(nstars_info_t stars);
void initializeMpiStarType(MPI_Datatype * datatype);

#endif /* _COLLISIONS_HELPERS_ */
