#include "mashm.h"

#include "mashmBool.h"
#include "intraNodeComm.h"
#include "mashmCommCycle.h"

struct MashmPrivate {
  MPI_Comm comm;
  int size;
  int rank;
  mashmBool isMasterProc;
  intraNodeComm intraComm;
  int numSharedMemNodes;
  int sharedMemIndex;
  int isInit;
  MashmCommCollection commCollection;
} ;

int mashmInit(Mashm* in_mashm, MPI_Comm in_comm) {
  int ierr;
  int numSharedMemNodes, sharedMemNodeRank;

  /* Temporary subcommunicator to determine the number of shared memory nodes */
  MPI_Comm rankComm; 

  in_mashm->p = (_p_mashm*) malloc(sizeof(_p_mashm));

  /* Set the communicator and get the size and rank */
  in_mashm->p->comm = in_comm;
  ierr = MPI_Comm_size(mashmGetComm(*in_mashm), &(in_mashm->p->size));
  ierr = MPI_Comm_rank(mashmGetComm(*in_mashm), &(in_mashm->p->rank));

  if (in_mashm->p->rank == 0) {
    in_mashm->p->isMasterProc = true;
  }
  else {
    in_mashm->p->isMasterProc = false;
  }

  /* Initialize the intra-node subcommunicator */
  init(&(in_mashm->p->intraComm),in_mashm->p->comm);

  /* Now calculate the number of shared memory indices */
  ierr = MPI_Comm_split(in_mashm->p->comm, in_mashm->p->intraComm.rank, in_mashm->p->rank, &rankComm);

  /* Only the nodal root is participates */
  if (in_mashm->p->intraComm.rank == 0) {
    ierr = MPI_Comm_size(rankComm, &numSharedMemNodes);
    ierr = MPI_Comm_rank(rankComm, &sharedMemNodeRank);
    /* The number of shared memory nodes */
    in_mashm->p->numSharedMemNodes = numSharedMemNodes;
    /* The index of each shared memory node */
    in_mashm->p->sharedMemIndex = sharedMemNodeRank;
    if (in_mashm->p->sharedMemIndex == 0) {
      printf("Number of shared memory nodes %d\n", in_mashm->p->numSharedMemNodes);
    }
  }

  /* Destroy this comm */
  ierr = MPI_Comm_free(&rankComm);

  /* Broadcast (to the shared sub comm) the number of shared memory nodes */
  ierr = MPI_Bcast(&(in_mashm->p->numSharedMemNodes), 1, MPI_INT, 0, in_mashm->p->comm);
  /* Broadcast (to the shared sub comm) the index of each shared memory nodes */
  ierr = MPI_Bcast(&(in_mashm->p->sharedMemIndex), 1, MPI_INT, 0, in_mashm->p->comm);

  /* Initialize the MashmCommCollection */
  MashmCommCollectionInit(&(in_mashm->p->commCollection));

  in_mashm->p->isInit = true;

  return 0;
}

MPI_Comm mashmGetComm(const Mashm in_mashm) {
  return in_mashm.p->comm;
}

int mashmGetSize(const Mashm in_mashm) {
  return in_mashm.p->size;
}

int mashmGetRank(const Mashm in_mashm) {
  return in_mashm.p->rank;
}

void mashmPrintInfo(const Mashm in_mashm) {
  int iNode;

  if (in_mashm.p->isMasterProc) {
    printf("Number of shared memory nodes %d\n", in_mashm.p->numSharedMemNodes);
  }

  for (iNode = 0; iNode < in_mashm.p->numSharedMemNodes; iNode++) {
    if (in_mashm.p->isMasterProc) {
      printf("  Node %d\n", iNode);
    }
    if (in_mashm.p->sharedMemIndex == iNode) {
      intraNodeCommPrintInfo(in_mashm.p->intraComm);
    }
  }
}

void mashmAddSymComm(Mashm in_mashm, int pairRank, int msgSize) {
  MashmCommCollectionAddComm(&(in_mashm.p->commCollection), pairRank, msgSize, msgSize);
}

/**
 * @brief Print out all of the communications
 */
void mashmPrintCommCollection(const Mashm in_mashm) {
  int i;
  for (i = 0; i < in_mashm.p->size; i++) {
    if (i == in_mashm.p->rank) {
      printf("Rank %d has communication:\n", in_mashm.p->rank);
      MashmCommCollectionPrint(in_mashm.p->commCollection);
    }
  }
}

/**
 * @brief Call when finished adding all messages.
 *
 * @param in_mashm Set precalculation of modified messaging
 */
void mashmCommFinish(Mashm in_mashm) {
}


/* @brief Begin nodal communication
 *
 * @param in_mash
 *
 * Begin Isend/Irecv communication. The waitalls for these are called with mashmInterNodeCommEnd
 */
void mashmInterNodeCommBegin(Mashm myMashm);

/* @brief Finish nodal communication
 *
 * @param in_mash
 *
 * Wait for all internode communication to be completed. Here, we call the MPI_Waitall corresponding to the MPI_Irecv/MPI_Isend calls in mashmInterNodeCommBegin.
 */
void mashmInterNodeCommEnd(Mashm myMashm);

/* @brief Perform intranode communication
 *
 * @param in_mash
 *
 * Perform intranode communication. Depending upon the method used this will call different algorithms.
 */
void mashmIntraNodeExchange(Mashm myMashm);

