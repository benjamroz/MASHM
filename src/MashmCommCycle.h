#ifndef MASHM_COMM_H
#define MASHM_COMM_H

#include "stdio.h"
#include "stdlib.h"

/* Structure representing one point to point communication */
typedef struct {
  int pairRank;
  int sendSize;
  int recvSize;
} MashmComm;

void MashmCommInit(MashmComm* in_comm) {
  in_comm->pairRank = -1;
  in_comm->sendSize = -1;
  in_comm->recvSize = -1;
}

void MashmCommSetComm(MashmComm* in_comm, int pairRank, int sendSize, int recvSize) {
  in_comm->pairRank = pairRank;
  in_comm->sendSize = sendSize;
  in_comm->recvSize = recvSize;
}

int MashmCommGetPairRank(const MashmComm in_comm) {
  return in_comm.pairRank;
}

int MashmCommGetSendSize(const MashmComm in_comm) {
  return in_comm.sendSize;
}

int MashmCommGetRecvSize(const MashmComm in_comm) {
  return in_comm.recvSize;
}

void MashmCommCopy(const MashmComm in_comm, MashmComm* copy_comm) {
  copy_comm->pairRank = MashmCommGetPairRank(in_comm);
  copy_comm->sendSize = MashmCommGetSendSize(in_comm);
  copy_comm->recvSize = MashmCommGetRecvSize(in_comm);
}

void MashmCommPrint(const MashmComm in_comm) {
  printf("Comm dest %d, send size %d, recv size %d\n", in_comm.pairRank, in_comm.sendSize, in_comm.recvSize);
}

/**
 * @brief Basically a struct to keep information relevant to one message 
 */
typedef struct {
  MashmComm* commArray;
  int commArraySize;
  MashmBool isInit;
  MashmBool isAllocated;
} MashmCommCollection;

/** 
 * @brief Do nothing routine for now
 */
void MashmCommCollectionInit(MashmCommCollection* commCollection) {

  commCollection->commArraySize = 0;
  commCollection->commArray = NULL;
  commCollection->isInit = true;
  commCollection->isAllocated = false;
}

/**
 * Set the size of the comm collection
 */
void MashmCommCollectionSetSize(MashmCommCollection* commCollection, int numComms) {

  /* Allocate at least two */
  if (commCollection->isAllocated) {
    free(commCollection->commArray);
    commCollection->isAllocated = false;
  }
  commCollection->commArraySize = numComms;
  commCollection->commArray = (MashmComm*) malloc(commCollection->commArraySize * sizeof(MashmComm));
  commCollection->isAllocated = true;
}

void MashmCommCollectionSetComm(MashmCommCollection* commCollection, int commIndex, int pairRank, int sendSize, int recvSize) {

  if (commIndex >= commCollection->commArraySize) {
    /* Comm Collection array is at maximum size
     * Need to extend the memory */
    printf("Error: Index %d outside of bounds of MashmCommCollection size %d\n", commIndex, commCollection->commArraySize);
  }
  (commCollection->commArray[commIndex]).pairRank = pairRank; 
  (commCollection->commArray[commIndex]).sendSize = sendSize;
  (commCollection->commArray[commIndex]).recvSize = recvSize;
}

void MashmCommCollectionPrint(const MashmCommCollection commCollection) {
  int i;
  if (!commCollection.isInit) {
    printf("MashmCommCollection not initialized\n");
    return;
  }
  for (i = 0; i < commCollection.commArraySize; i++) {
    MashmCommPrint(commCollection.commArray[i]); 
  }
}

void MashmCommCollectionReset(MashmCommCollection* commCollection) {

  commCollection->commArraySize = 0;
  free(commCollection->commArray);

  commCollection->isInit = true;
  commCollection->isAllocated = false;
}

void MashmCommCollectionDestroy(MashmCommCollection* commCollection) {

  free(commCollection->commArray);
  commCollection->isInit = false;
  commCollection->isAllocated = false;
  commCollection->commArraySize = 0;

}


#endif 
