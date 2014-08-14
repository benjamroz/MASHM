#include "intraNodeComm.h"

int intraNodeInit(intraNodeComm* intraComm, MPI_Comm in_comm) {
  int ierr;

  /* Store the parent commincator */
  intraComm->parentComm = in_comm;

  /* Create the shared memory subcommunicator groups */
  ierr = MPI_Comm_split_type(intraComm->parentComm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &(intraComm->comm));

  /* Determine the number of nodal comm groups */

  ierr = MPI_Comm_size(intraComm->comm, &(intraComm->size));
  ierr = MPI_Comm_rank(intraComm->comm, &(intraComm->rank));
  ierr = MPI_Comm_rank(intraComm->parentComm, &(intraComm->parentRank));

  /* Subcommunicator rank == 0 is the master process */
  if (intraComm->rank == 0) {
    intraComm->isMasterProc = true;
  }
  else {
    intraComm->isMasterProc = false;
  }

  intraComm->parentRanksOnNode = (int *) malloc(sizeof(int)*intraComm->size);

  ierr = MPI_Allgather(&(intraComm->parentRank), 1, MPI_INT, 
                       intraComm->parentRanksOnNode, 1, MPI_INT, intraComm->comm);

  return 0;
}

MPI_Comm intraNodeGetComm(const intraNodeComm intraComm) {
  return intraComm.comm;
}

int intraNodeGetSize(const intraNodeComm intraComm) {
  return intraComm.size;
}

int intraNodeGetRank(const intraNodeComm intraComm) {
  return intraComm.rank;
}

int intraNodeDetermineGlobalInfo(intraNodeComm* intraComm);
int intraNodeDetermineNodalInfo(intraNodeComm* intraComm);

void intraNodePrintInfo(const intraNodeComm intraComm) {
  if (intraComm.isMasterProc) {
    printf("  Shared memory node has size %d\n", intraComm.size);
  }
}
