#ifndef _COLLISIONS_HELPERS_
#define _COLLISIONS_HELPERS_

// TODO
#define PRINT_MSG_TAG 543
#define MPI_BACK_MESSAGE_TAG 1
#define MPI_FRONT_MESSAGE_TAG 2

// // TODO remove
// typedef struct star_extra_info_s {
//     int index;
//     int galaxy;
//     float velocity[2];
//     float acceleration[2];
// } star_extra_info_t;
//
// // information about stars positions
// typedef struct nstars_pos_s {
//     int n;
//     int galaxy;
//     float * starsPositions[2];
// } nstars_pos_t;

// extra information about stars velocities, acceleration and indices
typedef struct nstars_info_s {
    int n;
    int galaxy;
    bool onlyPositions;
    int * indices;
    float * starsPositions[2];
    float * starsVelocities[2];
    float * starsAccelerations[2];
} nstars_info_t;


// helper functions
void quicksort(int *A, int len);
int parseArguments(int argc, char const * argv[], char * filenameGal1, char * filenameGal2);
void initializeMpiStarType(MPI_Datatype * datatype);

#endif /* _COLLISIONS_HELPERS_ */
