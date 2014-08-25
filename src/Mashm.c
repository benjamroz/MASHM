#include "stddef.h"

#include "Mashm.h"

#include "MashmBool.h"
#include "intraNodeComm.h"
#include "MashmCommCycle.h"

#include "mpi.h"

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

typedef struct {
  int srcSharedMemRank;
  int destGlobalRank;
  int destNodeIndex;
  int msgSize;
  int srcNodeIndex;

} commTuple;

/* Sort on first element of comm tuple */
int commTupleCmpFunc (const void * a, const void * b) {
   int diff1;
   int diff2;
   int diff3;

   /* If (only) one of the destination nodes is the source node prefer that commTuple */
   if ((*(commTuple*)a).destNodeIndex != (*(commTuple*)b).destNodeIndex) {
     if ((*(commTuple*)a).destNodeIndex == (*(commTuple*)a).srcNodeIndex) {
       return -1;
     }
     else if ((*(commTuple*)b).destNodeIndex == (*(commTuple*)b).srcNodeIndex) {
       return 1;
     }
   }

   diff1 =  ( (*(commTuple*)a).destNodeIndex - (*(commTuple*)b).destNodeIndex );
   if (diff1 != 0) {
     return diff1;
   }
   else {
     diff2 = ( (*(commTuple*)a).srcSharedMemRank - (*(commTuple*)b).srcSharedMemRank );
     if (diff2 != 0) {
       return diff2;
     }
     else {
       diff3 = ( (*(commTuple*)a).destGlobalRank - (*(commTuple*)b).destGlobalRank );
       return diff3;
     }
   }
}


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
  int numNodalMsgs;
  int numCommNodes;

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

  /* Whether a message is on node of off node */
  MashmBool* onNodeMessage;

  int numInterNodePtrs;
  int numIntraNodePtrs;

  int numInterNodeMsgs;
  int numIntraNodeMsgs;

  MPI_Win sendSharedMemWindow;
  MPI_Win recvSharedMemWindow;

  double** sharedMemRecvBufferIndex;

  int* pairRanks;
  int* pairSharedRanks;

  int* sendIntraMsgOffsets;
  int* recvIntraMsgOffsets;

  int* msgOffsets;
  int* msgNodeIndices;
  int numNodalSubMsgs;
  int* nodalMsgSizes;

  int* nodalMsgOwner;
  MashmMinAggType minAggScheme;
  int numOwnedNodalMsgs;
  int nodalSharedBufferSize;
  int* nodalOffsets; /* Needed to offset nodal messages when one MPI process owns more than one nodal message */
  int* uniqueNodeIndices;
  int* nodalRecvRank;

  double* p_sendNodalSharedBuffer;
  double* p_recvNodalSharedBuffer;
  MPI_Win sendNodalSharedMemWindow;
  MPI_Win recvNodalSharedMemWindow;
  double** sendNodalSharedBufferIndex;
  double** recvNodalSharedBufferIndex;

};

/* Need a method to convert a Fortran MPI_Comm (integer) to  
 *   a C MPI_Comm type */
void MashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm) {
  MPI_Comm c_comm;
  c_comm = MPI_Comm_f2c(f_comm);
  MashmInit(in_mashm, c_comm);
}

void MashmInit(Mashm* in_mashm, MPI_Comm in_comm) {
  int ierr;
  int numSharedMemNodes, sharedMemNodeRank;

  /* Temporary subcommunicator to determine the number of shared memory nodes */
  MPI_Comm rankComm; 

  in_mashm->p = (_p_mashm*) malloc(sizeof(_p_mashm));

  /* Set the communicator and get the size and rank */
  in_mashm->p->comm = in_comm;

  ierr = MPI_Comm_size(MashmGetComm(*in_mashm), &(in_mashm->p->size));
  ierr = MPI_Comm_rank(MashmGetComm(*in_mashm), &(in_mashm->p->rank));

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

  in_mashm->p->minAggScheme = MASHM_MIN_AGG_ROUND_ROBIN;
  //return 0;
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
      intraNodePrintInfo(in_mashm.p->intraComm);
    }
  }
}

void MashmSetNumComms(Mashm in_mashm, int numComms) {
  MashmCommCollectionSetSize(&(in_mashm.p->commCollection), numComms);
}

void MashmSetComm(Mashm in_mashm, int commIndex, int pairRank, int msgSize) {
  MashmCommCollectionAddComm(&(in_mashm.p->commCollection), commIndex, pairRank, msgSize, msgSize);
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
  int regularBufferOffset, sharedBufferOffset;
  int iMsg, msgDestRank;
  int i;
  int numMessages;

  /* Allocate data for the number of messages */
  numMessages = in_mashm.p->commCollection.commArraySize;

  in_mashm.p->numOrigMessages = numMessages;

  in_mashm.p->sendBufferPointers = (double**) malloc(sizeof(double*)*in_mashm.p->numOrigMessages);
  in_mashm.p->recvBufferPointers = (double**) malloc(sizeof(double*)*in_mashm.p->numOrigMessages);
  in_mashm.p->onNodeMessage = (MashmBool*) malloc(sizeof(MashmBool)*in_mashm.p->numOrigMessages);
  in_mashm.p->pairRanks = (int*) malloc(sizeof(int)*numMessages);
  in_mashm.p->pairSharedRanks = (int*) malloc(sizeof(int)*numMessages);

  /* */
  for (i = 0; i < in_mashm.p->numOrigMessages; i++) {
    in_mashm.p->pairRanks[i] = (in_mashm.p->commCollection.commArray[i]).pairRank;
    in_mashm.p->pairSharedRanks[i] = intraNodeGetSharedRank(in_mashm.p->intraComm, in_mashm.p->pairRanks[i]);
  }

  /* Determine the number of intranode and internode messages and calculate the sizes */
  /* TODO: these two methods should be condensed into one -
   *       they are doing similar things 
   */
  MashmCalcNumMpiMsgs(in_mashm);
  MashmCalcMsgBufferSize(in_mashm);
  printf("Buf size %d, shared buf size %d\n", in_mashm.p->bufferSize,in_mashm.p->sharedBufferSize);
  printf("Rank %d, intra rank %d sends %d MPI messages\n", in_mashm.p->rank, in_mashm.p->intraComm.rank, in_mashm.p->numInterNodeMsgs);
  
  /* Allocate MPI (non-shared) buffer */
  in_mashm.p->p_regularSendBuffer = (double*) malloc(sizeof(double)*in_mashm.p->bufferSize);
  in_mashm.p->p_regularRecvBuffer = (double*) malloc(sizeof(double)*in_mashm.p->bufferSize);
  for (i = 0; i < in_mashm.p->bufferSize; i++) {
    in_mashm.p->p_regularSendBuffer[i] = 666;
    in_mashm.p->p_regularRecvBuffer[i] = 667;
  }

  /* Allocate MPI shared memory */
  if (in_mashm.p->commType == MASHM_COMM_INTRA_SHARED ||
      in_mashm.p->commType == MASHM_COMM_MIN_AGG) {
    p_mashmAllocateSharedMemory(in_mashm.p, in_mashm.p->sharedBufferSize);
  }

  for (i = 0; i < in_mashm.p->sharedBufferSize; i++) {
    in_mashm.p->p_sharedSendBuffer[i] = 668;
    in_mashm.p->p_sharedRecvBuffer[i] = 669;
  }
  in_mashm.p->buffersInit = true;

  /* Set pointers to the buffer and shared buffer */
  if (in_mashm.p->commType == MASHM_COMM_MIN_AGG) {
    p_mashmCalculateNodalMsgSchedule(in_mashm.p);
    p_mashmSetupAggType(in_mashm.p);
    p_mashmAllocateSharedMemoryMinAgg(in_mashm.p);
    p_mashmCalcMsgIndicesMinAgg(in_mashm.p);

  }
  else {
    regularBufferOffset = 0;
    sharedBufferOffset = 0;
    for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
      msgDestRank = (in_mashm.p->commCollection.commArray[iMsg]).pairRank;
      if ( (in_mashm.p->commType == MASHM_COMM_INTRA_SHARED ||
            in_mashm.p->commType == MASHM_COMM_MIN_AGG) 
          && MashmIsIntraNodeRank(in_mashm, msgDestRank)) {
        /* Set the shared memory pointers */
        in_mashm.p->sendBufferPointers[iMsg] = &(in_mashm.p->p_sharedSendBuffer[sharedBufferOffset]);
        in_mashm.p->recvBufferPointers[iMsg] = &(in_mashm.p->p_sharedRecvBuffer[sharedBufferOffset]);
        sharedBufferOffset = sharedBufferOffset + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
        in_mashm.p->onNodeMessage[iMsg] = true;
      }
      else {
        /* Set the pointers to regular memory */
        in_mashm.p->sendBufferPointers[iMsg] = &(in_mashm.p->p_regularSendBuffer[regularBufferOffset]);
        in_mashm.p->recvBufferPointers[iMsg] = &(in_mashm.p->p_regularRecvBuffer[regularBufferOffset]);
        regularBufferOffset = regularBufferOffset + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
        if (in_mashm.p->commType == MASHM_COMM_INTRA_MSG &&
            MashmIsIntraNodeRank(in_mashm, msgDestRank)) {
          /* On node MPI_Isend/MPI_Irecv */
          in_mashm.p->onNodeMessage[iMsg] = true;
        }
        else {
          in_mashm.p->onNodeMessage[iMsg] = false;
        }
      }
    }
  }

  /* All communication types need to setup the data needed for MPI_Irecv/MPI_Isend */
  p_mashmSetupInterNodeComm(in_mashm.p);

  /* Perform some initialization for each comm method and set pointers to proper routines */
  switch (in_mashm.p->commType) {
    case MASHM_COMM_STANDARD:
      in_mashm.p->p_interNodeCommBegin = p_mashmStandardCommBegin;
      in_mashm.p->p_interNodeCommEnd = p_mashmStandardCommEnd;

      in_mashm.p->p_intraNodeCommBegin = p_nullFunction;
      in_mashm.p->p_intraNodeCommEnd = p_nullFunction;

      break;

    case MASHM_COMM_INTRA_MSG:
      p_mashmSetupIntraMsgComm(in_mashm.p);

      in_mashm.p->p_interNodeCommBegin = p_mashmStandardCommBegin;
      in_mashm.p->p_interNodeCommEnd = p_mashmStandardCommEnd;

      in_mashm.p->p_intraNodeCommBegin = p_mashmIntraMsgsCommBegin;
      in_mashm.p->p_intraNodeCommEnd = p_mashmIntraMsgsCommEnd;

      break;

    case MASHM_COMM_INTRA_SHARED:
      /* Calculate shared memory indices etc. */
      p_mashmSetupIntraSharedComm(in_mashm.p);

      in_mashm.p->p_interNodeCommBegin = p_mashmStandardCommBegin;
      in_mashm.p->p_interNodeCommEnd = p_mashmStandardCommEnd;

      in_mashm.p->p_intraNodeCommBegin = p_mashmIntraSharedCommBegin;
      in_mashm.p->p_intraNodeCommEnd = p_mashmIntraSharedCommEnd;

      break;

    case MASHM_COMM_MIN_AGG:
      /* Calculate */
      p_mashmSetupIntraSharedComm(in_mashm.p);

      in_mashm.p->p_interNodeCommBegin = p_mashmMinAggCommBegin;
      in_mashm.p->p_interNodeCommEnd = p_mashmStandardCommEnd;

      in_mashm.p->p_intraNodeCommBegin = p_mashmIntraSharedCommBegin;
      in_mashm.p->p_intraNodeCommEnd = p_mashmIntraSharedCommEnd;

      break;

  }

}

MashmBool MashmIsMsgOnNode(Mashm in_mashm, int msgIndex) {
  return in_mashm.p->onNodeMessage[msgIndex];

}

double* MashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive) {
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

/**
 * Should be private
 */
void p_mashmSetupInterNodeComm(struct MashmPrivate* p_mashm) {
  int numMsgs;
  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
    numMsgs = p_mashm->numOwnedNodalMsgs;
  }
  else {
    numMsgs = p_mashm->numInterNodeMsgs;
  }
  /* Allocate the MPI data */
  p_mashm->recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numMsgs);
  p_mashm->sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numMsgs);
  p_mashm->recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numMsgs);
  p_mashm->sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numMsgs);
}

/**
 * Should be private
 */
void p_mashmSetupIntraMsgComm(struct MashmPrivate* p_mashm) {
  int numIntraNodeMsgs = p_mashm->numIntraNodeMsgs;

  /* Allocate the intranode MPI data */
  p_mashm->intraRecvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  p_mashm->intraSendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  p_mashm->intraRecvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  p_mashm->intraSendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  
}


void p_mashmSetupIntraSharedComm(struct MashmPrivate* p_mashm) {
  /* Need to figure out the indices and ... for the other ...*/
  int iMsg, msgCounter;
  int ierr;
  MPI_Aint otherWinSize;
  int otherDispUnit;

  p_mashm->sharedMemRecvBufferIndex = (double**) malloc(sizeof(double*)*p_mashm->numIntraNodeMsgs);

  /* Iterate through the intranodal message and get a pointer to the dest data */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
      /* TODO: Need shared rank here - NOT global pair Rank */
      ierr = MPI_Win_shared_query(p_mashm->sendSharedMemWindow, 
                                  p_mashm->pairSharedRanks[iMsg], &otherWinSize, 
                                  &otherDispUnit, (void*)&(p_mashm->sharedMemRecvBufferIndex[msgCounter]));
      if (ierr != 0) {
        printf("Error in MPI_Win_shared_query\n");
      }

    }
  }
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

void p_mashmStandardCommBegin(struct MashmPrivate* p_mashm) {
  int ierr;
  int iMsg;
  int numMsgs = p_mashm->commCollection.commArraySize;
  int i;
  double* buf;
  int msgCounter;

  /* First do the Irecvs */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (! p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
#if 0
      buf = p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE);
      for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
        if (buf[i] != 667 && buf[i] != 669) {
          printf("Error with recv buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
        }
      }
#endif

      ierr = MPI_Irecv(p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE), 
                       p_mashm->commCollection.commArray[iMsg].recvSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->recvRequests[msgCounter]));
      if (ierr != 0) {
        printf("Error in receiving message %d\n", iMsg);
      }
    }
  }
  /* Next do the Isends */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (! p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
#if 0
      buf = p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND);
      for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
        if (buf[i] != 666 && buf[i] != 668) {
          printf("Error with send buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
        }
      }
#endif
      ierr = MPI_Isend(p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND), 
                       p_mashm->commCollection.commArray[iMsg].sendSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->sendRequests[msgCounter]));
      if (ierr != 0) {
        printf("Error in sending message %d\n", iMsg);
      }
    }
  }
}

/* @brief Finish nodal communication
 *
 * @param in_mash
 *
 * Wait for all internode communication to be completed. Here, we call the MPI_Waitall corresponding to the MPI_Irecv/MPI_Isend calls in MashmInterNodeCommBegin.
 */
void p_mashmStandardCommEnd(struct MashmPrivate* p_mashm) {
  int ierr;

  int numMsgs;
  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
    numMsgs = p_mashm->numOwnedNodalMsgs;
  }
  else {
    numMsgs = p_mashm->numInterNodeMsgs;
  }
 
  ierr = MPI_Waitall(numMsgs, p_mashm->recvRequests, 
                     p_mashm->recvStatuses);
  ierr = MPI_Waitall(numMsgs, p_mashm->sendRequests, 
                     p_mashm->sendStatuses);
}

void p_mashmIntraMsgsCommBegin(struct MashmPrivate* p_mashm) {
  int ierr;
  int iMsg;
  int numMsgs = p_mashm->commCollection.commArraySize;
  int i;
  double* buf;
  int msgCounter;

  /* First do the Irecvs */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
#if 0
      buf = p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE);
      for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
        if (buf[i] != 667 && buf[i] != 669) {
          printf("Error with recv buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
        }
      }
#endif

      ierr = MPI_Irecv(p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE), 
                       p_mashm->commCollection.commArray[iMsg].recvSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->intraRecvRequests[msgCounter]));
      if (ierr != 0) {
        printf("Error in receiving message %d\n", iMsg);
      }
    }
  }

  /* Next do the Isends */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
#if 0
      buf = p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND);
      for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
        if (buf[i] != 666 && buf[i] != 668) {
          printf("Error with send buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
        }
      }
#endif

      ierr = MPI_Isend(p_mashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND), 
                       p_mashm->commCollection.commArray[iMsg].sendSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->intraSendRequests[msgCounter]));
      if (ierr != 0) {
        printf("Error in sending message %d\n", iMsg);
      }

    }
  }
}


void p_mashmIntraMsgsCommEnd(struct MashmPrivate* p_mashm) {
  int ierr;
 
  ierr = MPI_Waitall(p_mashm->numIntraNodeMsgs, p_mashm->intraRecvRequests, 
                     p_mashm->intraRecvStatuses);
  ierr = MPI_Waitall(p_mashm->numIntraNodeMsgs, p_mashm->intraSendRequests, 
                     p_mashm->intraSendStatuses);
}


void p_mashmIntraSharedCommBegin(struct MashmPrivate* p_mashm) {
  int iMsg, msgCounter;
  int ierr;
  int i;
  int mpiFenceMode;
  int sharedBufferOffset;

  mpiFenceMode = MPI_MODE_NOPUT;


  ierr = MPI_Win_fence(MPI_MODE_NOPUT,p_mashm->sendSharedMemWindow);

  /* Iterate through the intranodal message and get a pointer to the dest data */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
      sharedBufferOffset = p_mashm->recvIntraMsgOffsets[msgCounter];
      for (i = 0; i < (p_mashm->commCollection.commArray[iMsg]).sendSize; i++) {
        /* Copy data from the shared mem rank's send buffer to this recvBuffer */
        *(p_mashm->recvBufferPointers[iMsg]+i) =  *(p_mashm->sharedMemRecvBufferIndex[msgCounter]+sharedBufferOffset+i);

      }
    }
  }
}

void p_mashmIntraSharedCommEnd(struct MashmPrivate* p_mashm) {
  int ierr;
  ierr = MPI_Win_fence(MPI_MODE_NOSTORE,p_mashm->sendSharedMemWindow);
}

/* @brief Blank function for function pointers
 *
 * This will be called for communication methods that do no nothing in certain steps
 */
void p_nullFunction(struct MashmPrivate* p_mashm) {
}

/* @brief determine how many MPI messages that will be sent
 */
void MashmCalcNumMpiMsgs(Mashm in_mashm) {

  int tmpInt;
  int numNodalMessages;
  int floorNumMessagesPerRank;
  int modNumMessagesPerRank;

  /* This is usually zero */

  switch (in_mashm.p->commType) {
    case MASHM_COMM_STANDARD:
      in_mashm.p->numIntraNodeMsgs = 0;
      in_mashm.p->numInterNodeMsgs = in_mashm.p->commCollection.commArraySize;
      break;
    case MASHM_COMM_INTRA_MSG:
      in_mashm.p->numIntraNodeMsgs = MashmNumIntraNodeMsgs(in_mashm);
      in_mashm.p->numInterNodeMsgs = in_mashm.p->commCollection.commArraySize - in_mashm.p->numIntraNodeMsgs;
      break;
    case MASHM_COMM_INTRA_SHARED:
      in_mashm.p->numIntraNodeMsgs = MashmNumIntraNodeMsgs(in_mashm);
      in_mashm.p->numInterNodeMsgs = in_mashm.p->commCollection.commArraySize - in_mashm.p->numIntraNodeMsgs;
      break;
    case MASHM_COMM_MIN_AGG:
      /* Shared memory messages */
      in_mashm.p->numIntraNodeMsgs = MashmNumIntraNodeMsgs(in_mashm);
      /* Each node sends in_mashm.p->numSharedMemNodes MPI Messages 
       * Assign MPI messages in a round robin fashion */
      numNodalMessages = in_mashm.p->numSharedMemNodes;
      floorNumMessagesPerRank = numNodalMessages/in_mashm.p->intraComm.size;
      modNumMessagesPerRank = numNodalMessages % in_mashm.p->intraComm.size;

      /* All ranks have at least this many messages */
      in_mashm.p->numInterNodeMsgs = floorNumMessagesPerRank;

      /* Lower ranks may have one extra MPI message */
      if (in_mashm.p->intraComm.rank < modNumMessagesPerRank) {
        in_mashm.p->numInterNodeMsgs = floorNumMessagesPerRank + 1;
      }

      break;
  }
}

int MashmNumIntraNodeMsgs(Mashm in_mashm) {
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
    if (MashmIsIntraNodeRank(in_mashm, msgDestRank)) {
      numIntraNodeMsgs = numIntraNodeMsgs + 1;
    }
  }
  return numIntraNodeMsgs;
}

MashmBool MashmIsIntraNodeRank(Mashm in_mashm, int pairRank) {
  return p_MashmIsIntraNodeRank(in_mashm.p, pairRank);
}

void MashmCalcMsgBufferSize(Mashm in_mashm) {
  int iMsg;
  int iRank;
  int bufferSize;
  int sharedBufferSize;
  int tmpSize;
  int msgDestRank;
  int numIntraNodePtrs = 0;
  int numInterNodePtrs = 0;

  bufferSize = 0;
  sharedBufferSize = 0;

  /* Cycle through messages
   * Compare rank with ranks of nodal communicator 
   */
  for (iMsg = 0; iMsg < in_mashm.p->commCollection.commArraySize; iMsg++) {
    msgDestRank = (in_mashm.p->commCollection.commArray[iMsg]).pairRank;
    if ((in_mashm.p->commType == MASHM_COMM_INTRA_SHARED ||
         in_mashm.p->commType == MASHM_COMM_MIN_AGG) &&
        MashmIsIntraNodeRank(in_mashm, msgDestRank)) {
      sharedBufferSize = sharedBufferSize + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
      numIntraNodePtrs = numIntraNodePtrs + 1;

    }
    else {
      bufferSize = bufferSize + (in_mashm.p->commCollection.commArray[iMsg]).sendSize;
      numInterNodePtrs = numInterNodePtrs + 1;
    }
  }

  /* If MIN_AGG then all memory is shared */
#if 0
  if (in_mashm.p->commType == MASHM_COMM_MIN_AGG) {
    /* Swap the buffer sizes */
    tmpSize = bufferSize;
    bufferSize = sharedBufferSize;
    sharedBufferSize = tmpSize;

    /* Swap the number of pointers */
    tmpSize = numInterNodePtrs;
    numInterNodePtrs = numIntraNodePtrs;
    numIntraNodePtrs = tmpSize;
  }
#endif
  in_mashm.p->bufferSize = bufferSize;
  in_mashm.p->sharedBufferSize = sharedBufferSize;
  in_mashm.p->numInterNodePtrs = numInterNodePtrs;
  in_mashm.p->numIntraNodePtrs = numIntraNodePtrs;
}

void MashmDestroy(Mashm* in_mashm) {
  int ierr;

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
    free(in_mashm->p->onNodeMessage);
    free(in_mashm->p->pairRanks);
    free(in_mashm->p->pairSharedRanks);
  }

  /* Deallocate the shared memory window created with the
   *   call to MPI_Win_allocate_shared
   *   Note that this also frees the underlying shared memory */
  if (in_mashm->p->commType == MASHM_COMM_INTRA_SHARED ||
      in_mashm->p->commType == MASHM_COMM_MIN_AGG) {
    ierr = MPI_Win_free(&(in_mashm->p->sendSharedMemWindow));
    ierr = MPI_Win_free(&(in_mashm->p->recvSharedMemWindow));

  }

  /* Destroy the MashmPrivate data */
  free(in_mashm->p);

  //return 0;
}

/* 
 * Allocates contiguous shared memory data for the send and receive buffers
 */
void p_mashmAllocateSharedMemory(struct MashmPrivate* p_mashm, int bufferSize) {
  int ierr;
  int i;
  int numOrigMsgs;
  int counter, runningOffset;
  int numIntraNodeMsgs;

  MPI_Request *recvRequests, *sendRequests;
  MPI_Status *recvStatuses, *sendStatuses;

  ierr = MPI_Win_allocate_shared(sizeof(double)*bufferSize, sizeof(double),
                                 MPI_INFO_NULL, p_mashm->intraComm.comm,
                                 &(p_mashm->p_sharedSendBuffer),&(p_mashm->sendSharedMemWindow));
  if (ierr != 0) {
    printf("Error in MPI_Win_allocate_shared1\n");
  }

  ierr = MPI_Win_allocate_shared(sizeof(double)*bufferSize, sizeof(double),
                                 MPI_INFO_NULL, p_mashm->intraComm.comm,
                                 &(p_mashm->p_sharedRecvBuffer),&(p_mashm->recvSharedMemWindow));
  if (ierr != 0) {
    printf("Error in MPI_Win_allocate_shared2\n");
  }

  /* Initial synchronization of this memory */
  ierr = MPI_Win_fence(0,p_mashm->sendSharedMemWindow);
  ierr = MPI_Win_fence(0,p_mashm->recvSharedMemWindow);

  /** 
   * Next we need to communicate the offsets into this buffer of each 
   * submessage with the communication partner 
   */
  numIntraNodeMsgs = p_mashm->numIntraNodeMsgs;

  p_mashm->sendIntraMsgOffsets = (int*) malloc(sizeof(int)*numIntraNodeMsgs);
  p_mashm->recvIntraMsgOffsets = (int*) malloc(sizeof(int)*numIntraNodeMsgs);

  numOrigMsgs = p_mashm->numOrigMessages;

  counter = 0;
  runningOffset = 0;
  for (i = 0; i < numOrigMsgs; i++) {
    if (p_MashmIsIntraNodeRank(p_mashm, (p_mashm->commCollection.commArray[i]).pairRank)) {
      p_mashm->sendIntraMsgOffsets[counter] = runningOffset;
      counter += 1;
      runningOffset += (p_mashm->commCollection.commArray[i]).sendSize;
    }
  }

  /* Now exchange with neighbors */
  recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);

  /* Usual point to point communication */
  counter = 0;
  for (i = 0; i < numOrigMsgs; i++) {
    if (p_MashmIsIntraNodeRank(p_mashm, (p_mashm->commCollection.commArray[i]).pairRank)) {
      ierr = MPI_Irecv(&(p_mashm->recvIntraMsgOffsets[counter]), 1, MPI_INT, p_mashm->commCollection.commArray[i].pairRank, 0, p_mashm->comm, &recvRequests[counter]);
      counter += 1;
    }
  }
  counter = 0;
  for (i = 0; i < numOrigMsgs; i++) {
    if (p_MashmIsIntraNodeRank(p_mashm, (p_mashm->commCollection.commArray[i]).pairRank)) {
      ierr = MPI_Isend(&(p_mashm->sendIntraMsgOffsets[counter]), 1, MPI_INT, p_mashm->commCollection.commArray[i].pairRank, 0, p_mashm->comm, &sendRequests[counter]);
      counter += 1;
    }
  }

  ierr = MPI_Waitall(numIntraNodeMsgs,recvRequests,recvStatuses); 
  ierr = MPI_Waitall(numIntraNodeMsgs,sendRequests,sendStatuses); 
}

void p_mashmAllocateSharedMemoryMinAgg(struct MashmPrivate* p_mashm) {

  int ierr;
  int allocateSize;
  int iMsg, msgCounter;
  MPI_Aint otherWinSize;
  int otherDispUnit;

  allocateSize = p_mashm->nodalSharedBufferSize;
  /* Some MPI implementations don't like size zero arrays */
  if (allocateSize == 0) {
    allocateSize = 1;
  }

  ierr = MPI_Win_allocate_shared(sizeof(double)*allocateSize, sizeof(double),
                                 MPI_INFO_NULL, p_mashm->intraComm.comm,
                                 &(p_mashm->p_sendNodalSharedBuffer),&(p_mashm->sendNodalSharedMemWindow));
  if (ierr != 0) {
    printf("Error in MPI_Win_allocate_shared1\n");
  }

  ierr = MPI_Win_allocate_shared(sizeof(double)*allocateSize, sizeof(double),
                                 MPI_INFO_NULL, p_mashm->intraComm.comm,
                                 &(p_mashm->p_recvNodalSharedBuffer),&(p_mashm->recvNodalSharedMemWindow));
  if (ierr != 0) {
    printf("Error in MPI_Win_allocate_shared2\n");
  }

  /* Initial synchronization of this memory */
  ierr = MPI_Win_fence(0,p_mashm->sendNodalSharedMemWindow);
  ierr = MPI_Win_fence(0,p_mashm->recvNodalSharedMemWindow);

  /* Queries... */

  /* Need to figure out the indices and ... for the other ...*/
  p_mashm->sendNodalSharedBufferIndex = (double**) malloc(sizeof(double*)*p_mashm->numNodalMsgs);
  p_mashm->recvNodalSharedBufferIndex = (double**) malloc(sizeof(double*)*p_mashm->numNodalMsgs);
  printf("MINAGGG\n");
  /* Iterate through the intranodal message and get a pointer to the dest data */
  for (iMsg = 0; iMsg < p_mashm->numNodalMsgs; iMsg++) {
    /* TODO: Need shared rank here - NOT global pair Rank */
    ierr = MPI_Win_shared_query(p_mashm->sendNodalSharedMemWindow, 
                                p_mashm->nodalMsgOwner[iMsg], &otherWinSize, 
                                &otherDispUnit, (void*)&(p_mashm->sendNodalSharedBufferIndex[iMsg]));
    if (ierr != 0) {
      printf("Error in MPI_Win_shared_query\n");
    }
    ierr = MPI_Win_shared_query(p_mashm->recvNodalSharedMemWindow, 
                                p_mashm->nodalMsgOwner[iMsg], &otherWinSize, 
                                &otherDispUnit, (void*)&(p_mashm->recvNodalSharedBufferIndex[iMsg]));
    if (ierr != 0) {
      printf("Error in MPI_Win_shared_query\n");
    }
  }
}

void p_mashmCalculateNodalMsgSchedule(struct MashmPrivate* p_mashm) {
  int i;
  int ierr;
  int numOrigMsgs;
  int* msgNodeIndex;
  int* msgNodeOffset;

  MPI_Request* recvRequests;
  MPI_Request* sendRequests;
  MPI_Status* recvStatuses;
  MPI_Status* sendStatuses;

  int* srcSharedMemRanks;
  int* destGlobalRanks;
  int* msgSizes;
  int sumNumMsgs;

  int numNodes;
  int nodeCounter;
  int* allMsgNodeIndices;
  int* allSrcSharedMemRanks;
  int* allDestGlobalRanks;
  int* allMsgSizes;
  int iNodeIndex;

  int iRank;
  int* allNumMsgs;
  int* displs;

  commTuple* commArray;
  int* msgOffsets;
  int tmpMsgSize, j;
  int commArrayOffset;
  int iOffset;


  numOrigMsgs = p_mashm->commCollection.commArraySize;

  msgNodeIndex = (int*) malloc(sizeof(int)*numOrigMsgs);
  msgNodeOffset = (int*) malloc(sizeof(int)*numOrigMsgs);

  /* Need a communication cycle to get exchange node Indices with each pairRank */
  p_mashm->msgNodeIndices = (int*) malloc(sizeof(int)*numOrigMsgs);
  srcSharedMemRanks = (int*) malloc(sizeof(int)*numOrigMsgs);
  destGlobalRanks = (int*) malloc(sizeof(int)*numOrigMsgs);
  msgSizes = (int*) malloc(sizeof(int)*numOrigMsgs);

  for (i = 0; i < numOrigMsgs; i++) {
    srcSharedMemRanks[i] = p_mashm->intraComm.rank;
    destGlobalRanks[i] = p_mashm->commCollection.commArray[i].pairRank;
    msgSizes[i] = p_mashm->commCollection.commArray[i].sendSize;
  }

  // MPI_Isend, MPI_Irecv to retrieve nodeIndex
  recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numOrigMsgs);
  sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numOrigMsgs);
  recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numOrigMsgs);
  sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numOrigMsgs);

  /* Usual point to point communication */
  for (i = 0; i < numOrigMsgs; i++) {
    ierr = MPI_Irecv(&(p_mashm->msgNodeIndices[i]), 1, MPI_INT, p_mashm->commCollection.commArray[i].pairRank, 0, p_mashm->comm, &recvRequests[i]);
  }
  for (i = 0; i < numOrigMsgs; i++) {
    ierr = MPI_Isend(&(p_mashm->sharedMemIndex), 1, MPI_INT, p_mashm->commCollection.commArray[i].pairRank, 0, p_mashm->comm, &sendRequests[i]);
  }

  ierr = MPI_Waitall(numOrigMsgs,recvRequests,recvStatuses); 
  ierr = MPI_Waitall(numOrigMsgs,sendRequests,sendStatuses); 

  /* Reduce to shared memory root */
  if (p_mashm->intraComm.isMasterProc) {
    allNumMsgs = (int*) malloc(sizeof(int)*p_mashm->size);
  }
  else {
    /* Debugger throws error unless malloc'ed */
    allNumMsgs = (int*) malloc(sizeof(int)*0);
  }
  ierr = MPI_Gather(&numOrigMsgs, 1, MPI_INT, allNumMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);

  if (p_mashm->intraComm.isMasterProc) {
    sumNumMsgs = 0;
    for (i = 0; i < p_mashm->intraComm.size; i++) {
      sumNumMsgs += allNumMsgs[i];
    }
    allMsgNodeIndices = (int*) malloc(sizeof(int)*sumNumMsgs);
    allSrcSharedMemRanks = (int*) malloc(sizeof(int)*sumNumMsgs);
    allDestGlobalRanks = (int*) malloc(sizeof(int)*sumNumMsgs);
    allMsgSizes = (int*) malloc(sizeof(int)*sumNumMsgs);
    displs = (int*) malloc(sizeof(int)*sumNumMsgs);
    displs[0] = 0;
    for (i = 1; i < p_mashm->intraComm.size; i++) {
      displs[i] = displs[i-1] + allNumMsgs[i-1];
    }
  }
  else { 
    allMsgNodeIndices = (int*) malloc(sizeof(int)*0);
    displs = (int*) malloc(sizeof(int)*0);
    allSrcSharedMemRanks = NULL;
    allDestGlobalRanks = NULL;
    allMsgSizes = NULL;
  }
  ierr = MPI_Gatherv(p_mashm->msgNodeIndices, numOrigMsgs, MPI_INT,
                     allMsgNodeIndices, allNumMsgs, displs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Gatherv(srcSharedMemRanks, numOrigMsgs, MPI_INT,
                     allSrcSharedMemRanks, allNumMsgs, displs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Gatherv(destGlobalRanks, numOrigMsgs, MPI_INT,
                     allDestGlobalRanks, allNumMsgs, displs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Gatherv(msgSizes, numOrigMsgs, MPI_INT,
                     allMsgSizes, allNumMsgs, displs, MPI_INT, 0, p_mashm->intraComm.comm);

  
  /* sort list */
  if (p_mashm->intraComm.isMasterProc) {

    commArray = (commTuple*) malloc(sizeof(commTuple)*sumNumMsgs);

    /* Combine into a struct for the call to qsort */
    for (i = 0; i < sumNumMsgs; i++) {
      commArray[i].srcSharedMemRank = allSrcSharedMemRanks[i];
      commArray[i].destGlobalRank = allDestGlobalRanks[i];
      commArray[i].destNodeIndex = allMsgNodeIndices[i];
      commArray[i].msgSize = allMsgSizes[i];
      commArray[i].srcNodeIndex = p_mashm->sharedMemIndex;
    }


    qsort(commArray, sumNumMsgs, sizeof(commTuple), commTupleCmpFunc);

    /* Print the comm Schedule */
    for (i = 0; i < sumNumMsgs; i++) {
      printf("Rank %d: msg %i: srcRank %d, nodeIndex %d, destRank %d, size %d\n",
             p_mashm->rank, i, commArray[i].srcSharedMemRank, commArray[i].destNodeIndex,
             commArray[i].destGlobalRank, commArray[i].msgSize);
    }

    /* Now advance the commArray to the first non-self node */
    for (i = 0; i < sumNumMsgs; i++) {
      if (commArray[i].destNodeIndex != p_mashm->sharedMemIndex) {
        commArrayOffset = i;
        break;
      }
    }

    /* Count unique nodes to which we'll be messaging */
    nodeCounter = 0;
    for (i = commArrayOffset; i < sumNumMsgs; i++) {
      if (i == commArrayOffset) {
        nodeCounter += 1;
      }
      else if (commArray[i].destNodeIndex != commArray[i-1].destNodeIndex) {
        nodeCounter += 1;
      }
    }
    p_mashm->numNodalSubMsgs = sumNumMsgs - commArrayOffset;
    printf("Rank %d, commArrayOffset %d, numNodalSubMsgs %d\n", p_mashm->rank,
           commArrayOffset, p_mashm->numNodalSubMsgs);

    p_mashm->numNodalMsgs = nodeCounter;
    numNodes = p_mashm->numNodalMsgs;

    /* Count again including self node */
    nodeCounter = 0;
    for (i = 0; i < sumNumMsgs; i++) {
      if (i == 0) {
        nodeCounter += 1;
      }
      else if (commArray[i].destNodeIndex != commArray[i-1].destNodeIndex) {
        nodeCounter += 1;
      }
    }
    p_mashm->numCommNodes = nodeCounter;

    printf("Rank %d sends messages to num %d nodes\n", p_mashm->rank, numNodes);

    p_mashm->uniqueNodeIndices = (int*) malloc(sizeof(int)*(p_mashm->numCommNodes));

    nodeCounter = 0;
    p_mashm->uniqueNodeIndices[nodeCounter] = commArray[i].destNodeIndex;
    nodeCounter = nodeCounter + 1;
    for (i = 1; i < sumNumMsgs; i++) {
      if (commArray[i].destNodeIndex != commArray[i-1].destNodeIndex) {
        p_mashm->uniqueNodeIndices[nodeCounter] = commArray[i].destNodeIndex;
        nodeCounter = nodeCounter + 1;
      }
    }

    /* Calculate the offsets into each nodal message */
    /* Also calculate the size of each nodal message */
    msgOffsets = (int*) malloc(sizeof(int)*p_mashm->numNodalSubMsgs);
    //for (i = 0; i < commArrayOffset; i++) {
    //  msgOffsets[i] = -1;
    //}

    p_mashm->nodalMsgSizes = (int*) malloc(sizeof(int)*p_mashm->numNodalMsgs);

    nodeCounter = 0;
    msgOffsets[0] = 0;
    tmpMsgSize = commArray[commArrayOffset].msgSize;
    for (i = 1; i < p_mashm->numNodalSubMsgs; i++) {
      iOffset = i+commArrayOffset;
      if (commArray[iOffset].destNodeIndex == commArray[iOffset-1].destNodeIndex) {
        msgOffsets[i] = tmpMsgSize;
        tmpMsgSize += commArray[iOffset].msgSize;
      }
      else {
        printf("Okay %d, %d\n", commArray[iOffset].destNodeIndex, commArray[iOffset-1].destNodeIndex);
        printf("  Okay2 %d, %d\n", iOffset, iOffset-1);
        p_mashm->nodalMsgSizes[nodeCounter] = tmpMsgSize;

        tmpMsgSize = 667;
        msgOffsets[i] = tmpMsgSize;

        nodeCounter += 1;
      }
    }
    p_mashm->nodalMsgSizes[nodeCounter] = tmpMsgSize;
    for (i = 0; i < p_mashm->numNodalMsgs; i++) {
      printf("Rankn %d nodal msg %d size %d\n", p_mashm->rank, i, p_mashm->nodalMsgSizes[i]);

    }

    for (i = 0; i < sumNumMsgs; i++) {
      allSrcSharedMemRanks[i] = commArray[i].srcSharedMemRank;
      allDestGlobalRanks[i] = commArray[i].destGlobalRank;
      allMsgNodeIndices[i] = commArray[i].destNodeIndex;
    }

    if (p_mashm->intraComm.isMasterProc) {
      /* Print the comm Schedule */
      for (i = 0; i < sumNumMsgs; i++) {
        printf("Rankm %d: msg %i: srcRank %d, nodeIndex %d, destRank %d, size %d\n",
               p_mashm->rank, i, commArray[i].srcSharedMemRank, commArray[i].destNodeIndex,
               commArray[i].destGlobalRank, commArray[i].msgSize);
        if (i >= commArrayOffset) {
          printf("  offset %d\n", msgOffsets[i-commArrayOffset]);
          printf("    um %d, %d\n", allMsgNodeIndices[i], allMsgNodeIndices[i-1]);
        }
      }
    }

  }

  /* Now scatter the sorted list of srcSharedRanks and destGlobalRanks */
  ierr = MPI_Bcast(&sumNumMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&commArrayOffset, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&p_mashm->numNodalSubMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&p_mashm->numNodalMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&p_mashm->numCommNodes, 1, MPI_INT, 0, p_mashm->intraComm.comm);

  if (!p_mashm->intraComm.isMasterProc) {
    allSrcSharedMemRanks = (int*) malloc(sizeof(int)*sumNumMsgs);
    allDestGlobalRanks = (int*) malloc(sizeof(int)*sumNumMsgs);
    allMsgNodeIndices = (int*) malloc(sizeof(int)*sumNumMsgs);
    msgOffsets = (int*) malloc(sizeof(int)*p_mashm->numNodalSubMsgs);
    p_mashm->nodalMsgSizes = (int*) malloc(sizeof(int)*p_mashm->numNodalMsgs);
    p_mashm->uniqueNodeIndices = (int*) malloc(sizeof(int)*(p_mashm->numCommNodes));
  }

  /* Broadcast the msgOffsets and the nodalMsgsSizes */
  ierr = MPI_Bcast(allSrcSharedMemRanks, sumNumMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(allDestGlobalRanks, sumNumMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(allMsgNodeIndices, sumNumMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(msgOffsets, p_mashm->numNodalSubMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(p_mashm->nodalMsgSizes, p_mashm->numNodalMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(p_mashm->uniqueNodeIndices, p_mashm->numCommNodes, MPI_INT, 0, p_mashm->intraComm.comm);

  p_mashm->msgOffsets = (int*) malloc(sizeof(int)*p_mashm->numOrigMessages);

  /* Each rank goes through the lists, finds their message, and the appropriate offset */
  for (i = commArrayOffset; i < sumNumMsgs; i++) {
    if (allSrcSharedMemRanks[i] == p_mashm->intraComm.rank) {
      for (j = 0; j < p_mashm->numOrigMessages; j++) {
        if (p_mashm->commCollection.commArray[j].pairRank == allDestGlobalRanks[i]) {
          p_mashm->msgOffsets[j] = msgOffsets[i];
        }
      }
    }
  }
  return;
}

void p_mashmSetupAggType(struct MashmPrivate* p_mashm) {
  int numNodalMsgs;
  int numNodalMsgsOnRank;
  int i;
  int nodalMsgOffset;
  int ierr;

  MPI_Request* tmpRecvRequests;
  MPI_Request* tmpSendRequests;
  MPI_Status* tmpRecvStatuses;
  MPI_Status* tmpSendStatuses;
  int sharedRankMsgOwner;
  int globalRankMsgOwner;

  /* TODO: ensure that the commType is MASHM_COMM_MIN_AGG */
  numNodalMsgs = p_mashm->numNodalMsgs; 

  p_mashm->nodalSharedBufferSize = 0;

  p_mashm->nodalMsgOwner = (int*) malloc(sizeof(int)*numNodalMsgs);
 
  switch (p_mashm->minAggScheme) {
    case(MASHM_MIN_AGG_ROUND_ROBIN):
        p_mashm->numOwnedNodalMsgs = 0;

        for (i = 0; i < numNodalMsgs; i++) {
          p_mashm->nodalMsgOwner[i] = i % p_mashm->intraComm.size;
          if (p_mashm->nodalMsgOwner[i] == p_mashm->intraComm.rank) {
            p_mashm->nodalSharedBufferSize += p_mashm->nodalMsgSizes[i];
            p_mashm->numOwnedNodalMsgs += 1;
          }
        }
        break;
    case(MASHM_MIN_AGG_ROOT_PROC):
      /* Root process in shared memory node owns all messages */
      for (i = 0; i < numNodalMsgs; i++) {
        p_mashm->nodalMsgOwner[i] = 0;
      }
      if (p_mashm->intraComm.isMasterProc) {
        p_mashm->numOwnedNodalMsgs = numNodalMsgs;
        for (i = 0; i < numNodalMsgs; i++) {
          p_mashm->nodalSharedBufferSize += p_mashm->nodalMsgSizes[i];
        }
      }
      else {
        p_mashm->numOwnedNodalMsgs = 0;
        /* Some MPI implementations can't allocate buffers of zero size */
      }
      break;
  }

  /* Now each rank needs to figure out the nodal offsets */
  p_mashm->nodalOffsets = (int*) malloc(sizeof(int)*numNodalMsgs);
  nodalMsgOffset = 0;
  for (i = 0; i < numNodalMsgs; i++) {
    if (p_mashm->nodalMsgOwner[i] == p_mashm->intraComm.rank) {
      p_mashm->nodalOffsets[i] = nodalMsgOffset;
    }
    ierr = MPI_Bcast(&p_mashm->nodalOffsets[i], 1, MPI_INT, p_mashm->nodalMsgOwner[i], p_mashm->intraComm.comm);
    if (p_mashm->nodalMsgOwner[i] == p_mashm->intraComm.rank) {
      nodalMsgOffset += p_mashm->nodalMsgSizes[i];
    }
  }

  
  /* Now we need to exchange the rank of the nodal message owner with the rank of the receiver message owner */
  p_mashm->nodalRecvRank = (int*) malloc(sizeof(int)*numNodalMsgs);
  /* For this we need the nodal Comm again */

  MPI_Comm rankComm; 

  ierr = MPI_Comm_split(p_mashm->comm, p_mashm->intraComm.rank, p_mashm->rank, &rankComm);
 
  /* Only the nodal root is participates */
  if (p_mashm->intraComm.rank == 0) {

    tmpRecvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNodalMsgs);
    tmpSendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNodalMsgs);

    tmpRecvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNodalMsgs);
    tmpSendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNodalMsgs);

    /* in this ... the rank is the node index ... */
    /* Exchange the ranks with the */
    for (i = 0; i < numNodalMsgs; i++) {
      sharedRankMsgOwner = p_mashm->nodalMsgOwner[i];
      if (sharedRankMsgOwner == p_mashm->intraComm.rank) {
        globalRankMsgOwner = p_mashm->intraComm.parentRanksOnNode[sharedRankMsgOwner];
        ierr = MPI_Irecv(&(p_mashm->nodalRecvRank[i]), 1, MPI_INT,
                         p_mashm->uniqueNodeIndices[i+1], 1, rankComm, &(tmpRecvRequests[i])); 
      }
    }
    for (i = 0; i < numNodalMsgs; i++) {
      sharedRankMsgOwner = p_mashm->nodalMsgOwner[i];
      if (sharedRankMsgOwner == p_mashm->intraComm.rank) {
        globalRankMsgOwner = p_mashm->intraComm.parentRanksOnNode[sharedRankMsgOwner];
        ierr = MPI_Isend(&globalRankMsgOwner, 1, MPI_INT,
                         p_mashm->uniqueNodeIndices[i+1], 1, rankComm, &(tmpSendRequests[i])); 
      }
    }
    ierr = MPI_Waitall(numNodalMsgs, tmpRecvRequests, tmpSendStatuses);
    ierr = MPI_Waitall(numNodalMsgs, tmpSendRequests, tmpSendStatuses);
    free(tmpRecvRequests);
    free(tmpSendRequests);
    free(tmpRecvStatuses);
    free(tmpSendStatuses);
  }

  /* Destroy this rank comm */
  ierr = MPI_Comm_free(&rankComm);

  ierr = MPI_Bcast(p_mashm->nodalRecvRank, numNodalMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
}


void p_mashmCalcMsgIndicesMinAgg(struct MashmPrivate* p_mashm) {
  int sendNodalSharedBufferOffset;
  int nodalMsgOffset;
  int msgOffset;
  int nodalMsgDest;
  int nodalMsgIndex;
  int i, iMsg;
  int sharedBufferOffset;

  sendNodalSharedBufferOffset = 0;
  sharedBufferOffset = 0;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    nodalMsgDest = p_mashm->msgNodeIndices[iMsg];
    /* If intranode continue */
    if (nodalMsgDest == p_mashm->sharedMemIndex) {
      /* Set the shared memory pointers */
      p_mashm->sendBufferPointers[iMsg] = &(p_mashm->p_sharedSendBuffer[sharedBufferOffset]);
      p_mashm->recvBufferPointers[iMsg] = &(p_mashm->p_sharedRecvBuffer[sharedBufferOffset]);
      sharedBufferOffset = sharedBufferOffset + (p_mashm->commCollection.commArray[iMsg]).sendSize;
      p_mashm->onNodeMessage[iMsg] = true;
    }
    else {
      /* Now find the index of the nodalMsgDest */
      /* Ignore the 0th entry since it is sharedMemIndex */
      for (i = 1; i < p_mashm->numCommNodes; i++) {
        if (p_mashm->uniqueNodeIndices[i] == nodalMsgDest) {
          break; 
        }
      }
      /* Minus one because of the difference in indexing between
       *   uniqueNodeIndices and nodalOffsets */
      nodalMsgIndex = i-1;
      nodalMsgOffset = p_mashm->nodalOffsets[nodalMsgIndex];
      msgOffset = p_mashm->msgOffsets[iMsg];
      p_mashm->sendBufferPointers[iMsg] = (p_mashm->sendNodalSharedBufferIndex[nodalMsgIndex]+nodalMsgOffset + msgOffset);
      p_mashm->recvBufferPointers[iMsg] = (p_mashm->recvNodalSharedBufferIndex[nodalMsgIndex]+nodalMsgOffset + msgOffset);
      p_mashm->onNodeMessage[iMsg] = false;
    }
  }

  /* TODO: now need to exchange with comm... to get the offset of the receive buffers */
}


void p_mashmMinAggCommBegin(struct MashmPrivate* p_mashm) {
  int ierr;
  int iMsg;
  int numMsgs = p_mashm->commCollection.commArraySize;
  int i;
  double* buf;
  int msgCounter, sharedRankMsgOwner, globalRankMsgOwner;

  /* First do the Irecvs */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->numNodalMsgs; iMsg++) {
    sharedRankMsgOwner = p_mashm->nodalMsgOwner[iMsg];
    if (sharedRankMsgOwner == p_mashm->intraComm.rank) {
      msgCounter += 1;
      globalRankMsgOwner = p_mashm->intraComm.parentRanksOnNode[sharedRankMsgOwner];
      ierr = MPI_Irecv(&(p_mashm->p_recvNodalSharedBuffer[p_mashm->nodalOffsets[iMsg]]),
                       p_mashm->nodalMsgSizes[iMsg], MPI_DOUBLE,
                       p_mashm->nodalRecvRank[iMsg],
                       1, p_mashm->comm, &(p_mashm->recvRequests[msgCounter]));
      if (ierr != 0) {
        printf("Error in MPI_Irecv message%d\n", iMsg);
      }
    }
  }
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->numNodalMsgs; iMsg++) {
    sharedRankMsgOwner = p_mashm->nodalMsgOwner[iMsg];
    if (sharedRankMsgOwner == p_mashm->intraComm.rank) {
      msgCounter += 1;
      globalRankMsgOwner = p_mashm->intraComm.parentRanksOnNode[sharedRankMsgOwner];
      ierr = MPI_Isend(&(p_mashm->p_sendNodalSharedBuffer[p_mashm->nodalOffsets[iMsg]]),
                       p_mashm->nodalMsgSizes[iMsg], MPI_DOUBLE,
                       p_mashm->nodalRecvRank[iMsg],
                       1, p_mashm->comm, &(p_mashm->sendRequests[msgCounter]));
      if (ierr != 0) {
        printf("Error in MPI_Isend message%d\n", iMsg);
      }
    }
  }
}


MashmBool p_MashmIsIntraNodeRank(struct MashmPrivate* p_mashm, int pairRank) {
  int iRank;
  for (iRank = 0; iRank < p_mashm->intraComm.size; iRank++) {
    if (pairRank == p_mashm->intraComm.parentRanksOnNode[iRank]) {
      return true;
    }
  }
  return false;
}

