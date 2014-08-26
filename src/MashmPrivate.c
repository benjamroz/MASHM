
#include "MashmPrivate.h"


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

/**
 * Should be private
 */
double* p_mashmGetBufferPointer(struct MashmPrivate* p_mashm, int msgIndex, MashmSendReceive sendReceive) {
  if (sendReceive == MASHM_SEND) {
    return p_mashm->sendBufferPointers[msgIndex];
  }
  else {
    return p_mashm->recvBufferPointers[msgIndex];
  }
}


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
  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
    ierr = MPI_Win_fence(0,p_mashm->sendNodalSharedMemWindow);
    ierr = MPI_Win_fence(0,p_mashm->recvNodalSharedMemWindow);
  }
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
void p_MashmCalcNumMpiMsgs(struct MashmPrivate* p_mashm) {

  int tmpInt;
  int numNodalMessages;
  int floorNumMessagesPerRank;
  int modNumMessagesPerRank;

  /* This is usually zero */

  switch (p_mashm->commType) {
    case MASHM_COMM_STANDARD:
      p_mashm->numIntraNodeMsgs = 0;
      p_mashm->numInterNodeMsgs = p_mashm->commCollection.commArraySize;
      break;
    case MASHM_COMM_INTRA_MSG:
      p_mashm->numIntraNodeMsgs = p_MashmNumIntraNodeMsgs(p_mashm);
      p_mashm->numInterNodeMsgs = p_mashm->commCollection.commArraySize - p_mashm->numIntraNodeMsgs;
      break;
    case MASHM_COMM_INTRA_SHARED:
      p_mashm->numIntraNodeMsgs = p_MashmNumIntraNodeMsgs(p_mashm);
      p_mashm->numInterNodeMsgs = p_mashm->commCollection.commArraySize - p_mashm->numIntraNodeMsgs;
      break;
    case MASHM_COMM_MIN_AGG:
      /* Shared memory messages */
      p_mashm->numIntraNodeMsgs = p_MashmNumIntraNodeMsgs(p_mashm);
      /* Each node sends p_mashm->numSharedMemNodes MPI Messages 
       * Assign MPI messages in a round robin fashion */
      numNodalMessages = p_mashm->numSharedMemNodes;
      floorNumMessagesPerRank = numNodalMessages/p_mashm->intraComm.size;
      modNumMessagesPerRank = numNodalMessages % p_mashm->intraComm.size;

      /* All ranks have at least this many messages */
      p_mashm->numInterNodeMsgs = floorNumMessagesPerRank;

      /* Lower ranks may have one extra MPI message */
      if (p_mashm->intraComm.rank < modNumMessagesPerRank) {
        p_mashm->numInterNodeMsgs = floorNumMessagesPerRank + 1;
      }

      break;
  }
}

int p_MashmNumIntraNodeMsgs(struct MashmPrivate* p_mashm) {
  int iMsg;
  int iRank;
  int numIntraNodeMsgs;
  int msgDestRank;

  numIntraNodeMsgs = 0;
  /* Cycle through messages
   * Compare rank with ranks of nodal communicator 
   */
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    msgDestRank = (p_mashm->commCollection.commArray[iMsg]).pairRank;
    if (p_MashmIsIntraNodeRank(p_mashm, msgDestRank)) {
      numIntraNodeMsgs = numIntraNodeMsgs + 1;
    }
  }
  return numIntraNodeMsgs;
}

void p_MashmCalcMsgBufferSize(struct MashmPrivate* p_mashm) {
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
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    msgDestRank = (p_mashm->commCollection.commArray[iMsg]).pairRank;
    if ((p_mashm->commType == MASHM_COMM_INTRA_SHARED ||
         p_mashm->commType == MASHM_COMM_MIN_AGG) &&
        p_MashmIsIntraNodeRank(p_mashm, msgDestRank)) {
      sharedBufferSize = sharedBufferSize + (p_mashm->commCollection.commArray[iMsg]).sendSize;
      numIntraNodePtrs = numIntraNodePtrs + 1;

    }
    else {
      bufferSize = bufferSize + (p_mashm->commCollection.commArray[iMsg]).sendSize;
      numInterNodePtrs = numInterNodePtrs + 1;
    }
  }

  /* If MIN_AGG then all memory is shared */
#if 0
  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
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
  p_mashm->bufferSize = bufferSize;
  p_mashm->sharedBufferSize = sharedBufferSize;
  p_mashm->numInterNodePtrs = numInterNodePtrs;
  p_mashm->numIntraNodePtrs = numIntraNodePtrs;
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
  int counter;


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

  p_mashm->sendAggMsgOffsets = (int*) malloc(sizeof(int)*p_mashm->numOrigMessages);
  p_mashm->recvAggMsgOffsets = (int*) malloc(sizeof(int)*p_mashm->numOrigMessages);

  /* Each rank goes through the lists, finds their message, and the appropriate offset */
  for (i = commArrayOffset; i < sumNumMsgs; i++) {
    if (allSrcSharedMemRanks[i] == p_mashm->intraComm.rank) {
      for (j = 0; j < p_mashm->numOrigMessages; j++) {
        if (p_mashm->commCollection.commArray[j].pairRank == allDestGlobalRanks[i]) {
          p_mashm->sendAggMsgOffsets[j] = msgOffsets[i-commArrayOffset];
        }
      }
    }
  }

  /* Now exchange with sendAggMsgOffset with neighbors */

  recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numOrigMsgs);
  sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numOrigMsgs);
  recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numOrigMsgs);
  sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numOrigMsgs);

  /* Usual point to point communication */
  counter = 0;
  for (i = 0; i < numOrigMsgs; i++) {
    if (! p_MashmIsIntraNodeRank(p_mashm, (p_mashm->commCollection.commArray[i]).pairRank)) {
      ierr = MPI_Irecv(&(p_mashm->recvAggMsgOffsets[i]), 1, MPI_INT, p_mashm->commCollection.commArray[i].pairRank, 0, p_mashm->comm, &recvRequests[counter]);
      counter += 1;
    }
  }
  counter = 0;
  for (i = 0; i < numOrigMsgs; i++) {
    if (! p_MashmIsIntraNodeRank(p_mashm, (p_mashm->commCollection.commArray[i]).pairRank)) {
      ierr = MPI_Isend(&(p_mashm->sendAggMsgOffsets[i]), 1, MPI_INT, p_mashm->commCollection.commArray[i].pairRank, 0, p_mashm->comm, &sendRequests[counter]);
      counter += 1;
    }
  }

  /* Symmetric */
  ierr = MPI_Waitall(counter,recvRequests,recvStatuses); 
  ierr = MPI_Waitall(counter,sendRequests,sendStatuses); 

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
  int sendMsgOffset, recvMsgOffset;
  int nodalMsgDest;
  int nodalMsgIndex;
  int i, iMsg;
  int sharedBufferOffset;

  MPI_Request *recvRequests, *sendRequests;
  MPI_Status *recvStatuses, *sendStatuses;

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
      sendMsgOffset = p_mashm->sendAggMsgOffsets[iMsg];
      recvMsgOffset = p_mashm->recvAggMsgOffsets[iMsg];
      p_mashm->sendBufferPointers[iMsg] = (p_mashm->sendNodalSharedBufferIndex[nodalMsgIndex]+nodalMsgOffset + sendMsgOffset);
      p_mashm->recvBufferPointers[iMsg] = (p_mashm->recvNodalSharedBufferIndex[nodalMsgIndex]+nodalMsgOffset + recvMsgOffset);
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

  /* MPI_Win_fence */
  ierr = MPI_Win_fence(0,p_mashm->sendNodalSharedMemWindow);
  ierr = MPI_Win_fence(0,p_mashm->recvNodalSharedMemWindow);

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


