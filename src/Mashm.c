#include "stddef.h"

#include "mpi.h"

#include "MashmIntraNodeComm.h"
#include "MashmCommCycle.h"
#include "MashmPrivate.h"
#include "Mashm.h"


int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

/* Need a method to convert a Fortran MPI_Comm (integer) to  
 *   a C MPI_Comm type */
void MashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm) {
  MPI_Comm c_comm;
  c_comm = MPI_Comm_f2c(f_comm);
  MashmInit(in_mashm, c_comm);
}

void MashmInit(Mashm* in_mashm, MPI_Comm in_comm) {

  /* Allocate the MashmPrivate structure */
  in_mashm->p = (struct MashmPrivate*) malloc(sizeof(struct MashmPrivate));

  /* Initialize the MashmPrivate structure */
  p_MashmInit(in_mashm->p, in_comm);

}

MPI_Comm MashmGetComm(const Mashm in_mashm) {
  return in_mashm.p->comm;
}

int MashmGetSize(const Mashm in_mashm) {
  return in_mashm.p->size;
}

int MashmGetRank(const Mashm in_mashm) {
  return in_mashm.p->rank;
}

void MashmPrintInfo(const Mashm in_mashm) {
  int iNode;

  if (in_mashm.p->isMasterProc) {
    printf("Number of shared memory nodes %d\n", in_mashm.p->numSharedMemNodes);
  }

  for (iNode = 0; iNode < in_mashm.p->numSharedMemNodes; iNode++) {
    if (in_mashm.p->isMasterProc) {
      printf("  Node %d\n", iNode);
    }
    if (in_mashm.p->sharedMemIndex == iNode) {
      MashmIntraNodeCommPrintInfo(in_mashm.p->intraComm);
    }
  }
}

void MashmSetNumComms(Mashm in_mashm, int numComms) {
  MashmCommCollectionSetSize(&(in_mashm.p->commCollection), numComms);
}

void MashmSetComm(Mashm in_mashm, int commIndex, int pairRank, int msgSize) {
  MashmCommCollectionSetComm(&(in_mashm.p->commCollection), commIndex, pairRank, msgSize, msgSize);
}

/**
 * @brief Print out all of the communications
 */
void MashmPrintCommCollection(const Mashm in_mashm) {
  int i;
  for (i = 0; i < in_mashm.p->size; i++) {
    if (i == in_mashm.p->rank) {
      printf("Rank %d has communication:\n", in_mashm.p->rank);
      MashmCommCollectionPrint(in_mashm.p->commCollection);
    }
  }
}

/**
 * @brief Set the communication type
 *
 * Default (if this is not called) is MASHM_COMM_STANDARD
 */
void MashmSetCommMethod(Mashm in_mashm, MashmCommType commType) {
  in_mashm.p->commType = commType;
}

/**
 * @brief Get the communication type
 */
MashmCommType MashmGetCommMethod(Mashm in_mashm) {
  return in_mashm.p->commType;
}

/**
 * @brief Call when finished adding all messages.
 * 
 * Allocate buffer data for communication methods. 
 * If intranode then just regular data
 * If intranode shared then 
 *     regular data for internode comm
 *     shared data for intranode comm
 * If minimal nodal then
 *     block of shared data for each nodal message
 *       or
 *     block of shared data for all messages
 *      
 * TODO: This function is kind of a mess - it needs to be simplified
 *
 * @param in_mashm Set precalculation of modified messaging
 */
void MashmCommFinish(Mashm in_mashm) {
  p_MashmFinish(in_mashm.p);
}

MashmBool MashmIsMsgOnNode(Mashm in_mashm, int msgIndex) {
  return in_mashm.p->onNodeMessage[msgIndex];

}

double* MashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive) {
  /* Call the private routine */
  return p_MashmGetBufferPointer(in_mashm.p, msgIndex, sendReceive);
}

void MashmRetireBufferPointer(Mashm in_mashm, double** bufPtr) {
  /* Call the private routine */
  return p_MashmRetireBufferPointer(in_mashm.p, bufPtr);
}

double* MashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive) {
  int iRank;
  for (iRank = 0; iRank < in_mashm.p->intraComm.size; iRank++) {
    if (destRank == in_mashm.p->intraComm.parentRanksOnNode[iRank]) {
      if (sendReceive == MASHM_SEND) {
        return in_mashm.p->sendBufferPointers[iRank];
      }
      else {
        return in_mashm.p->recvBufferPointers[iRank];
      }
    }
  }
  return NULL;
}

void MashmDestroy(Mashm* in_mashm) {

  /* Destroy the MashmPrivate data */
  p_Destroy(in_mashm->p);

  /* Free the pointer to the MashmPrivate data */
  free(in_mashm->p);

}

/* @brief Begin nodal communication
 *
 * @param in_mash
 *
 * Begin Isend/Irecv communication. The waitalls for these are called with MashmInterNodeCommEnd
 */
void MashmInterNodeCommBegin(Mashm in_mashm) {
  in_mashm.p->p_interNodeCommBegin(in_mashm.p);
}

void MashmIntraNodeCommBegin(Mashm in_mashm) {
  in_mashm.p->p_intraNodeCommBegin(in_mashm.p);
}

void MashmIntraNodeCommEnd(Mashm in_mashm) {
  in_mashm.p->p_intraNodeCommEnd(in_mashm.p);
}

void MashmInterNodeCommEnd(Mashm in_mashm) {
  in_mashm.p->p_interNodeCommEnd(in_mashm.p);
}

MashmBool MashmIsIntraNodeRank(Mashm in_mashm, int pairRank) {
  return p_MashmIsIntraNodeRank(in_mashm.p, pairRank);
}

