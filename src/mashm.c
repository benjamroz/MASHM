#include "stddef.h"

#include "mashm.h"

#include "mashmBool.h"
#include "intraNodeComm.h"
#include "mashmCommCycle.h"

#include "mpi.h"

struct MashmPrivate {
  /* Root communicator */
  MPI_Comm comm;
  int size;
  int rank;
  MashmBool isMasterProc;
  /* Intranodal communicator struct */
  intraNodeComm intraComm;
  int numSharedMemNodes;
  int sharedMemIndex;
  int isInit;
  MashmCommCollection commCollection;
  MashmCommType commType;
  /* Pointers to the underlying data. Used by user to fill buffer. */
  double** sendBufferPointers; 
  double** recvBufferPointers; 
  /* Actual allocated memory. Hidden from user, pointed to by bufferPointers. */
  double* p_regularSendBuffer; 
  double* p_regularRecvBuffer; 
  MashmBool p_regularBufferIsAlloc;
  /* Actual allocated shared MPI memory. Hidden from user, pointed to by bufferPointers. */
  double* p_sharedSendBuffer;
  double* p_sharedRecvBuffer;
  MashmBool p_sharedIsAlloc;

  int numOrigMessages;
  int numMpiMsgs;
  int bufferSize;
  int sharedBufferSize;
  int numConnectedNodes;

  /* MPI data */
  MPI_Request* recvRequests;
  MPI_Request* sendRequests;
  MPI_Status* recvStatuses;
  MPI_Status* sendStatuses;

  /* MPI data for intranode comm 
   * Should this go in intraComm?
   * */
  MPI_Request* intraRecvRequests;
  MPI_Request* intraSendRequests;
  MPI_Status* intraRecvStatuses;
  MPI_Status* intraSendStatuses;

};

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
  intraNodeInit(&(in_mashm->p->intraComm),in_mashm->p->comm);

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

  in_mashm->p->commType = MASHM_COMM_STANDARD;

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
      intraNodePrintInfo(in_mashm.p->intraComm);
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
 * @brief Set the communication type
 *
 * Default (if this is not called) is MASHM_COMM_STANDARD
 */
void mashmSetCommMethod(Mashm in_mashm, MashmCommType commType) {
  in_mashm.p->commType = commType;
}

/**
 * @brief Get the communication type
 */
MashmCommType mashmGetCommMethod(Mashm in_mashm) {
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
 * TODO: We should have an option to declare which communication strategy we are targeting.
 *
 * @param in_mashm Set precalculation of modified messaging
 */
void mashmCommFinish(Mashm in_mashm) {
  int regularBufferOffset, sharedBufferOffset;
  int iMsg, msgDestRank;
  
  /* Allocate data for the number of messages */
  in_mashm.p->numOrigMessages = in_mashm.p->commCollection.commArraySize;
  in_mashm.p->sendBufferPointers = (double**) malloc(sizeof(double*)*in_mashm.p->numOrigMessages);
  in_mashm.p->recvBufferPointers = (double**) malloc(sizeof(double*)*in_mashm.p->numOrigMessages);

  /* Determine the number of intranode and internode messages and calculate the sizes */
  mashmCalcNumMpiMsgs(in_mashm);
  mashmCalcMsgBufferSize(in_mashm);

  /* Allocate MPI buffer */
  in_mashm.p->p_regularSendBuffer = (double*) malloc(sizeof(double)*in_mashm.p->bufferSize);
  in_mashm.p->p_regularRecvBuffer = (double*) malloc(sizeof(double)*in_mashm.p->bufferSize);
  /* Allocate MPI shared memory */
  in_mashm.p->p_sharedSendBuffer = (double*) malloc(sizeof(double)*in_mashm.p->sharedBufferSize);
  in_mashm.p->p_sharedRecvBuffer = (double*) malloc(sizeof(double)*in_mashm.p->sharedBufferSize);

  /* Note that this depends upon which communication method we choose */
  switch (in_mashm.p->commType) {
    case MASHM_COMM_STANDARD:
      mashmSetupInterNodeComm(in_mashm);
      break;

    case MASHM_COMM_INTRA_MSG:
      mashmSetupIntraNodeComm(in_mashm);
      break;

    case MASHM_COMM_INTRA_SHARED:
      printf("Communication method MASHM_INTRA_SHARED not yet implemented\n");
      break;

    case MASHM_COMM_MIN_AGG:
      printf("Communication method MASHM_MIN_AGG not yet implemented\n");
      break;

  }

  /* Set pointers to the buffer and shared buffer */
  regularBufferOffset = 0;
  sharedBufferOffset = 0;
  for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
    msgDestRank = (in_mashm.p->commCollection.commArray[iMsg]).pairRank;
    if ( (in_mashm.p->commType == MASHM_COMM_INTRA_SHARED ||
          in_mashm.p->commType == MASHM_COMM_MIN_AGG) 
        && mashmIsIntraNodeRank(in_mashm, msgDestRank)) {
      in_mashm.p->sendBufferPointers[iMsg] = &(in_mashm.p->p_sharedSendBuffer[sharedBufferOffset]);
      in_mashm.p->recvBufferPointers[iMsg] = &(in_mashm.p->p_sharedRecvBuffer[sharedBufferOffset]);
      sharedBufferOffset = sharedBufferOffset + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
    }
    else {
      in_mashm.p->sendBufferPointers[iMsg] = &(in_mashm.p->p_regularSendBuffer[regularBufferOffset]);
      in_mashm.p->recvBufferPointers[iMsg] = &(in_mashm.p->p_regularRecvBuffer[regularBufferOffset]);
      regularBufferOffset = regularBufferOffset + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
    }
  }
}

double* mashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive) {
  if (sendReceive == MASHM_SEND) {
    return in_mashm.p->sendBufferPointers[msgIndex];
  }
  else {
    return in_mashm.p->recvBufferPointers[msgIndex];
  }
}

double* mashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive) {
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

void mashmSetupInterNodeComm(Mashm in_mashm) {
  int numMsgs = in_mashm.p->numMpiMsgs;
  /* Allocate the MPI data */
  in_mashm.p->recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numMsgs);
  in_mashm.p->sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numMsgs);
  in_mashm.p->recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numMsgs);
  in_mashm.p->sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numMsgs);
}

void mashmSetupIntraNodeComm(Mashm in_mashm) {
  int numMsgs = in_mashm.p->numMpiMsgs;

  int numIntraNodeMsgs = mashmNumIntraNodeMsgs(in_mashm);
  int numInterNodeMsgs = numMsgs - numIntraNodeMsgs;

  /* Allocate the internode MPI data */
  in_mashm.p->recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numInterNodeMsgs);
  in_mashm.p->sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numInterNodeMsgs);
  in_mashm.p->recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numInterNodeMsgs);
  in_mashm.p->sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numInterNodeMsgs);

  /* Allocate the intranode MPI data */
  in_mashm.p->intraRecvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  in_mashm.p->intraSendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  in_mashm.p->intraRecvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  in_mashm.p->intraSendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  
}


/* @brief Begin nodal communication
 *
 * @param in_mash
 *
 * Begin Isend/Irecv communication. The waitalls for these are called with mashmInterNodeCommEnd
 */
void mashmInterNodeCommBegin(Mashm in_mashm) {
  int tmp;
  switch (in_mashm.p->commType) {
    case (MASHM_COMM_STANDARD):
      mashmStandardCommBegin(in_mashm); 
      break;
    case (MASHM_COMM_INTRA_MSG):
      //mashmIntraMsgsCommBegin(in_mashm); 
      tmp = 1;
      break;
    case (MASHM_COMM_INTRA_SHARED):
      tmp = 1;
      break;
    case (MASHM_COMM_MIN_AGG):
      tmp = 1;
      break;
  }
}

void mashmIntraNodeCommBegin(Mashm in_mashm) {
}

void mashmIntraNodeCommEnd(Mashm in_mashm) {
}



void mashmInterNodeCommEnd(Mashm in_mashm) {
  int tmp;
  switch (in_mashm.p->commType) {
    case (MASHM_COMM_STANDARD):
      mashmStandardCommEnd(in_mashm); 
      break;
    case (MASHM_COMM_INTRA_MSG):
     // mashmIntraMsgsCommEnd(in_mashm); 
      tmp = 1;
      break;
    case (MASHM_COMM_INTRA_SHARED):
      tmp = 1;
      break;
    case (MASHM_COMM_MIN_AGG):
      tmp = 1;
      break;
  }
}

void mashmStandardCommBegin(Mashm in_mashm) {
  int ierr;
  int iMsg;

  /* First do the Irecvs */
  for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
    ierr = MPI_Irecv(mashmGetBufferPointer(in_mashm, iMsg, MASHM_RECEIVE), 
                     in_mashm.p->commCollection.commArray[iMsg].recvSize, 
                     MPI_DOUBLE,
                     in_mashm.p->commCollection.commArray[iMsg].pairRank, 
                     1, in_mashm.p->comm, &(in_mashm.p->recvRequests[iMsg]));

  }
  /* Next do the Isends */
  for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
    ierr = MPI_Isend(mashmGetBufferPointer(in_mashm, iMsg, MASHM_SEND), 
                     in_mashm.p->commCollection.commArray[iMsg].sendSize, 
                     MPI_DOUBLE,
                     in_mashm.p->commCollection.commArray[iMsg].pairRank, 
                     1, in_mashm.p->comm, &(in_mashm.p->sendRequests[iMsg]));

  }

}


void mashmIntraMsgsCommBegin(Mashm in_mashm) {
}

/* @brief Finish nodal communication
 *
 * @param in_mash
 *
 * Wait for all internode communication to be completed. Here, we call the MPI_Waitall corresponding to the MPI_Irecv/MPI_Isend calls in mashmInterNodeCommBegin.
 */
void mashmStandardCommEnd(Mashm in_mashm) {
  int ierr;
 
  ierr = MPI_Waitall(in_mashm.p->commCollection.commArraySize, in_mashm.p->recvRequests, 
                     in_mashm.p->recvStatuses);
  ierr = MPI_Waitall(in_mashm.p->commCollection.commArraySize, in_mashm.p->sendRequests, 
                     in_mashm.p->sendStatuses);
}

/* @brief Perform intranode communication
 *
 * @param in_mash
 *
 * Perform intranode communication. Depending upon the method used this will call different algorithms.
 */
void mashmIntraNodeExchange(Mashm myMashm);

/* @brief determine how many MPI messages that will be sent
 */
void mashmCalcNumMpiMsgs(Mashm in_mashm) {
  int numMpiMsgs;
  if (in_mashm.p->commType == MASHM_COMM_STANDARD ||
      in_mashm.p->commType == MASHM_COMM_INTRA_MSG) {
    numMpiMsgs = in_mashm.p->commCollection.commArraySize;
  }
  else if (in_mashm.p->commType == MASHM_COMM_INTRA_SHARED) {
    /* Determine how many communications are intranode and subtract for the total number */
    numMpiMsgs = in_mashm.p->commCollection.commArraySize - mashmNumIntraNodeMsgs(in_mashm);
  }
  else if (in_mashm.p->commType == MASHM_COMM_MIN_AGG) {
    mashmCalcNumConnectedNodes(in_mashm);
    numMpiMsgs = in_mashm.p->numConnectedNodes;
  }
  in_mashm.p->numMpiMsgs = numMpiMsgs;
}

int mashmNumIntraNodeMsgs(Mashm in_mashm) {
  int iMsg;
  int iRank;
  int numIntraNodeMsgs;
  int msgDestRank;

  numIntraNodeMsgs = 0;
  /* Cycle through messages
   * Compare rank with ranks of nodal communicator 
   */
  for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
    msgDestRank = (in_mashm.p->commCollection.commArray[iMsg]).pairRank;
    if (mashmIsIntraNodeRank(in_mashm, msgDestRank)) {
      numIntraNodeMsgs = numIntraNodeMsgs + 1;
    }
  }
  return numIntraNodeMsgs;
}

MashmBool mashmIsIntraNodeRank(Mashm in_mashm, int pairRank) {
  int iRank;
  for (iRank = 0; iRank < in_mashm.p->intraComm.size; iRank++) {
    if (pairRank == in_mashm.p->intraComm.parentRanksOnNode[iRank]) {
      return true;
    }
  }
  return false;
}

void mashmCalcMsgBufferSize(Mashm in_mashm) {
  int iMsg;
  int iRank;
  int bufferSize;
  int sharedBufferSize;
  int tmpSize;
  int msgDestRank;

  bufferSize = 0;
  sharedBufferSize = 0;
  /* Cycle through messages
   * Compare rank with ranks of nodal communicator 
   */
  for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
    msgDestRank = (in_mashm.p->commCollection.commArray[iMsg]).pairRank;
    if (in_mashm.p->commType == MASHM_COMM_INTRA_SHARED && 
        mashmIsIntraNodeRank(in_mashm, msgDestRank)) {
      sharedBufferSize = sharedBufferSize + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
    }
    else {
      bufferSize = bufferSize + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
    }
  }
  /* If MIN_AGG then all memory is shared */
  if (MASHM_COMM_MIN_AGG) {
    tmpSize = bufferSize;
    bufferSize = sharedBufferSize;
    sharedBufferSize = tmpSize;
  }

  in_mashm.p->bufferSize = bufferSize;
  in_mashm.p->sharedBufferSize = sharedBufferSize;
}

void mashmCalcNumConnectedNodes(Mashm in_mashm) {
  in_mashm.p->numConnectedNodes = -1;
}
