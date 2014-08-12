#ifndef MASHM_COMM_H
#define MASHM_COMM_H

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

typedef struct {
  MashmComm* commArray;
  int commArraySize;
  int commArrayReserveSize;
  mashmBool isInit;
} MashmCommCollection;

void MashmCommCollectionInit(MashmCommCollection* commCollection) {

  /* Allocate at least two */
  commCollection->commArrayReserveSize = 2;
  commCollection->commArray = (MashmComm*) malloc(commCollection->commArrayReserveSize * sizeof(MashmComm));
  commCollection->commArraySize = 0;

  commCollection->isInit = true;
}

void MashmCommCollectionExtend(MashmCommCollection* commCollection) {
  MashmComm* tempCollection;
  int i;
  int commArrayNewReserveSize = 2*commCollection->commArrayReserveSize;

  tempCollection = (MashmComm*) malloc(commArrayNewReserveSize*sizeof(MashmComm));

  for (i = 0; i < commCollection->commArraySize; i++) {
    MashmCommCopy(commCollection->commArray[i],&(tempCollection[i]));
  }
  
  free(commCollection->commArray);
  commCollection->commArray = tempCollection;
  commCollection->commArrayReserveSize = commArrayNewReserveSize;

}

void MashmCommCollectionAddComm(MashmCommCollection* commCollection, int pairRank, int sendSize, int recvSize) {
  int commIndex;

  if (commCollection->commArraySize == commCollection->commArrayReserveSize) {
    /* Comm Collection array is at maximum size
     * Need to extend the memory */
    MashmCommCollectionExtend(commCollection);
  }
  commIndex = commCollection->commArraySize;
  (commCollection->commArray[commIndex]).pairRank = pairRank; 
  (commCollection->commArray[commIndex]).sendSize = sendSize;
  (commCollection->commArray[commIndex]).recvSize = recvSize;
  commCollection->commArraySize = commCollection->commArraySize + 1;
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

#endif 
