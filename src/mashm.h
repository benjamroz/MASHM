#include "mpi.h"

#include "intraNodeComm.h"

typedef struct {
  MPI_Comm comm;
  int size;
  int rank;
  mashmBool isMasterProc;
  intraNodeComm intraComm;
  int numSharedMemNodes;
  int sharedMemIndex;
  int isInit;
} Mashm;

int mashmInit(Mashm* in_mashm, MPI_Comm in_comm);
MPI_Comm mashmGetComm(const Mashm in_mashm);
int mashmGetSize(const Mashm in_mashm);
int mashmGetRank(const Mashm in_mashm);
