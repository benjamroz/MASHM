#include "MashmIntraNodeComm.h"

int MashmIntraNodeCommInit(MashmIntraNodeComm* intraComm, MPI_Comm in_comm) {
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

MPI_Comm MashmIntraNodeCommGetComm(const MashmIntraNodeComm intraComm) {
  return intraComm.comm;
}

int MashmIntraNodeCommGetSize(const MashmIntraNodeComm intraComm) {
  return intraComm.size;
}

int MashmIntraNodeCommGetRank(const MashmIntraNodeComm intraComm) {
  return intraComm.rank;
}

void MashmIntraNodeCommPrintInfo(const MashmIntraNodeComm intraComm) {
  if (intraComm.isMasterProc) {
    printf("  Shared memory node has size %d\n", intraComm.size);
  }
}

int MashmIntraNodeCommGetSharedRank(const MashmIntraNodeComm intraComm, int pairRank) {
  int i;
  for (i = 0; i < intraComm.size; i++) {
    if (pairRank == intraComm.parentRanksOnNode[i]) {
      return i;
    }
  }
  return -1;
}

int MashmIntraNodeCommDestroy(MashmIntraNodeComm* intraComm) {
  int ierr;

  /* Free the shared memory subcommunicator groups */
  ierr = MPI_Comm_free(&(intraComm->comm));

  return 0;
}
