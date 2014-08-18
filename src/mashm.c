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
  MashmBool buffersInit;

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

  void (* p_interNodeCommBegin)(struct MashmPrivate*);
  void (* p_interNodeCommEnd)(struct MashmPrivate*);
  void (* p_intraNodeCommBegin)(struct MashmPrivate*);
  void (* p_intraNodeCommEnd)(struct MashmPrivate*);
};

/* Need a method to convert a Fortran MPI_Comm (integer) to  
 *   a C MPI_Comm type */
void mashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm) {
  MPI_Comm c_comm;
  c_comm = MPI_Comm_f2c(f_comm);
  mashmInit(in_mashm, c_comm);
}

void mashmInit(Mashm* in_mashm, MPI_Comm in_comm) {
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
  ierr = MPI_Bcast(&(in_mashm->p->numSharedMemNodes), 1, MPI_INT, 0, in_mashm->p->intraComm.comm);
  /* Broadcast (to the shared sub comm) the index of each shared memory nodes */
  ierr = MPI_Bcast(&(in_mashm->p->sharedMemIndex), 1, MPI_INT, 0, in_mashm->p->intraComm.comm);

  /* Initialize the MashmCommCollection */
  MashmCommCollectionInit(&(in_mashm->p->commCollection));

  in_mashm->p->commType = MASHM_COMM_STANDARD;

  in_mashm->p->isInit = true;
  in_mashm->p->buffersInit = false;

  //return 0;
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
  int i;
  
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
  for (i = 0; i < in_mashm.p->bufferSize; i++) {
    in_mashm.p->p_regularSendBuffer[i] = 666;
    in_mashm.p->p_regularRecvBuffer[i] = 667;
  }

  /* Allocate MPI shared memory */
  in_mashm.p->p_sharedSendBuffer = (double*) malloc(sizeof(double)*in_mashm.p->sharedBufferSize);
  in_mashm.p->p_sharedRecvBuffer = (double*) malloc(sizeof(double)*in_mashm.p->sharedBufferSize);
  for (i = 0; i < in_mashm.p->sharedBufferSize; i++) {
    in_mashm.p->p_sharedSendBuffer[i] = 668;
    in_mashm.p->p_sharedRecvBuffer[i] = 669;
  }
  in_mashm.p->buffersInit = true;

  /* Note that this depends upon which communication method we choose */
  switch (in_mashm.p->commType) {
    case MASHM_COMM_STANDARD:
      mashmSetupInterNodeComm(in_mashm);
      in_mashm.p->p_interNodeCommBegin = p_mashmStandardCommBegin;
      in_mashm.p->p_interNodeCommEnd = p_mashmStandardCommEnd;

      in_mashm.p->p_intraNodeCommBegin = p_nullFunction;
      in_mashm.p->p_intraNodeCommEnd = p_nullFunction;

      break;

    case MASHM_COMM_INTRA_MSG:
      mashmSetupIntraNodeComm(in_mashm);

      printf("Communication method MASHM_INTRA_MSGS not yet implemented\n");

      in_mashm.p->p_interNodeCommBegin = p_nullFunction;
      in_mashm.p->p_interNodeCommEnd = p_nullFunction;
      in_mashm.p->p_intraNodeCommBegin = p_nullFunction;
      in_mashm.p->p_intraNodeCommEnd = p_nullFunction;

      break;

    case MASHM_COMM_INTRA_SHARED:
      printf("Communication method MASHM_INTRA_SHARED not yet implemented\n");

      in_mashm.p->p_interNodeCommBegin = p_nullFunction;
      in_mashm.p->p_interNodeCommEnd = p_nullFunction;
      in_mashm.p->p_intraNodeCommBegin = p_nullFunction;
      in_mashm.p->p_intraNodeCommEnd = p_nullFunction;

      break;

    case MASHM_COMM_MIN_AGG:
      printf("Communication method MASHM_MIN_AGG not yet implemented\n");

      in_mashm.p->p_interNodeCommBegin = p_nullFunction;
      in_mashm.p->p_interNodeCommEnd = p_nullFunction;
      in_mashm.p->p_intraNodeCommBegin = p_nullFunction;
      in_mashm.p->p_intraNodeCommEnd = p_nullFunction;

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
  /* Call the private routine */
  return p_mashmGetBufferPointer(in_mashm.p, msgIndex, sendReceive);
}


double* p_mashmGetBufferPointer(struct MashmPrivate* p_mashm, int msgIndex, MashmSendReceive sendReceive) {
  if (sendReceive == MASHM_SEND) {
    return p_mashm->sendBufferPointers[msgIndex];
  }
  else {
    return p_mashm->recvBufferPointers[msgIndex];
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
  in_mashm.p->p_interNodeCommBegin(in_mashm.p);
}

void mashmIntraNodeCommBegin(Mashm in_mashm) {
  in_mashm.p->p_intraNodeCommBegin(in_mashm.p);
}

void mashmIntraNodeCommEnd(Mashm in_mashm) {
  in_mashm.p->p_intraNodeCommEnd(in_mashm.p);
}

void mashmInterNodeCommEnd(Mashm in_mashm) {
  in_mashm.p->p_interNodeCommEnd(in_mashm.p);
}

void p_mashmStandardCommBegin(struct MashmPrivate* p_mashm) {
  int ierr;
  int iMsg;
  int numMsgs = p_mashm->commCollection.commArraySize;
  int i;
  double* buf;

  /* First do the Irecvs */
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    buf = p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE);
    for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
      if (buf[i] != 667 && buf[i] != 669) {
        printf("Error with recv buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
      }
    }

    ierr = MPI_Irecv(p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE), 
                     p_mashm->commCollection.commArray[iMsg].recvSize, 
                     MPI_DOUBLE,
                     p_mashm->commCollection.commArray[iMsg].pairRank, 
                     1, p_mashm->comm, &(p_mashm->recvRequests[iMsg]));
    if (ierr != 0) {
      printf("Error in receiving message %d\n", iMsg);
    }

  }
  /* Next do the Isends */
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    buf = p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND);
    for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
      if (buf[i] != 666 && buf[i] != 668) {
        printf("Error with send buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
      }
    }

    ierr = MPI_Isend(p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND), 
                     p_mashm->commCollection.commArray[iMsg].sendSize, 
                     MPI_DOUBLE,
                     p_mashm->commCollection.commArray[iMsg].pairRank, 
                     1, p_mashm->comm, &(p_mashm->sendRequests[iMsg]));
    if (ierr != 0) {
      printf("Error in sending message %d\n", iMsg);
    }

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
void p_mashmStandardCommEnd(struct MashmPrivate* p_mashm) {
  int ierr;
 
  ierr = MPI_Waitall(p_mashm->commCollection.commArraySize, p_mashm->recvRequests, 
                     p_mashm->recvStatuses);
  ierr = MPI_Waitall(p_mashm->commCollection.commArraySize, p_mashm->sendRequests, 
                     p_mashm->sendStatuses);
}

/* @brief Blank function for function pointers
 *
 * This will be called for communication methods that do no nothing in certain steps
 */
void p_nullFunction(struct MashmPrivate* p_mashm) {
}

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
  if (in_mashm.p->commType == MASHM_COMM_MIN_AGG) {
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


void mashmDestroy(Mashm* in_mashm) {

  /* Destroy the MashmCommCollection 
   * TODO: this should be destroyed in the finish method */
  MashmCommCollectionDestroy(&(in_mashm->p->commCollection));

  /* Destroy the intra-node subcommunicator */
  intraNodeDestroy(&(in_mashm->p->intraComm));
 
  if (in_mashm->p->buffersInit) {
    free(in_mashm->p->sendBufferPointers);
    free(in_mashm->p->recvBufferPointers);
    free(in_mashm->p->p_regularSendBuffer);
    free(in_mashm->p->p_regularRecvBuffer);
    free(in_mashm->p->p_sharedSendBuffer);
    free(in_mashm->p->p_sharedRecvBuffer);
  }


  /* Destroy the MashmPrivate data */
  free(in_mashm->p);

  //return 0;
}


void mashmInitTmp(MPI_Fint f_comm) {
  int ierr;
  int tmpSize;
  MPI_Comm c_comm;
  c_comm = MPI_Comm_f2c(f_comm);
  ierr = MPI_Comm_size(c_comm, &tmpSize);
  printf("mashmInitTmp size = %d\n", tmpSize);
}

