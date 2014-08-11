#include "mashm.h"

int mashmInit(Mashm* in_mashm, MPI_Comm in_comm) {
  int ierr;
  int numSharedMemNodes, sharedMemNodeRank;

  /* Temporary subcommunicator to determine the number of shared memory nodes */
  MPI_Comm rankComm; 

  /* Set the communicator and get the size and rank */
  in_mashm->comm = in_comm;
  ierr = MPI_Comm_size(mashmGetComm(*in_mashm), &(in_mashm->size));
  ierr = MPI_Comm_rank(mashmGetComm(*in_mashm), &(in_mashm->rank));

  if (in_mashm->rank == 0) {
    in_mashm->isMasterProc = true;
  }
  else {
    in_mashm->isMasterProc = false;
  }

  /* Initialize the intra-node subcommunicator */
  init(&(in_mashm->intraComm),in_mashm->comm);

  /* Now calculate the number of shared memory indices */
  ierr = MPI_Comm_split(in_mashm->comm, in_mashm->intraComm.rank, in_mashm->rank, &rankComm);

  /* Only the nodal root is participates */
  if (in_mashm->intraComm.rank == 0) {
    ierr = MPI_Comm_size(rankComm, &numSharedMemNodes);
    ierr = MPI_Comm_rank(rankComm, &sharedMemNodeRank);
    /* The number of shared memory nodes */
    in_mashm->numSharedMemNodes = numSharedMemNodes;
    /* The index of each shared memory node */
    in_mashm->sharedMemIndex = sharedMemNodeRank;
    if (in_mashm->sharedMemIndex == 0) {
      printf("Number of shared memory nodes %d\n", in_mashm->numSharedMemNodes);
    }
  }

  /* Destroy this comm */
  ierr = MPI_Comm_free(&rankComm);

  /* Broadcast (to the shared sub comm) the number of shared memory nodes */
  ierr = MPI_Bcast(&(in_mashm->numSharedMemNodes), 1, MPI_INT, 0, in_mashm->comm);
  /* Broadcast (to the shared sub comm) the index of each shared memory nodes */
  ierr = MPI_Bcast(&(in_mashm->sharedMemIndex), 1, MPI_INT, 0, in_mashm->comm);

  in_mashm->isInit = true;

  return 0;
}

MPI_Comm mashmGetComm(const Mashm in_mashm) {
  return in_mashm.comm;
}

int mashmGetSize(const Mashm in_mashm) {
  return in_mashm.size;
}

int mashmGetRank(const Mashm in_mashm) {
  return in_mashm.rank;
}

void mashmPrintInfo(const Mashm in_mashm) {
  int iNode;

  if (in_mashm.isMasterProc) {
    printf("Number of shared memory nodes %d\n", in_mashm.numSharedMemNodes);

  }

  for (iNode = 0; iNode < in_mashm.numSharedMemNodes; iNode++) {
    if (in_mashm.isMasterProc) {
      printf("  Node %d\n", iNode);
    }
    if (in_mashm.sharedMemIndex == iNode) {
      intraNodeCommPrintInfo(in_mashm.intraComm);
    }

  }


}
