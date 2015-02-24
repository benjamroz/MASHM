#include "MashmPrivate.h"

#if defined(MASHM_DEBUG)
#define MASHM_DEBUG_PRINT(fmt, args...) fprintf(stderr, "MASHM_DEBUG: %s:%d:%s(): " fmt, \
   __FILE__, __LINE__, __func__, ##args)
#else
#define DEBUG_PRINT(fmt, args...) /* Don't do anything in release builds */
#endif

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
/* datatype required to use qsort */
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

  /* Sort the destination nodes */
  diff1 =  ( (*(commTuple*)a).destNodeIndex - (*(commTuple*)b).destNodeIndex );
  if (diff1 != 0) {
   return diff1;
  }
  else {
   /* If the destination nodes are the same, first sort by sender rank 
    *   Note this has implications for shared memory */
   diff2 = ( (*(commTuple*)a).srcSharedMemRank - (*(commTuple*)b).srcSharedMemRank );
   if (diff2 != 0) {
     return diff2;
   }
   else {
     /* If the sender is the same sort by destination rank */
     diff3 = ( (*(commTuple*)a).destGlobalRank - (*(commTuple*)b).destGlobalRank );
     return diff3;
   }
  }
}


void p_MashmInit(struct MashmPrivate* p_mashm, MPI_Comm in_comm) {
  int ierr;
  int sharedMemNodeRank;
  int numSharedMemNodes;
  /* Set the communicator and get the size and rank */
  p_mashm->comm = in_comm;

  ierr = MPI_Comm_size(p_mashm->comm, &(p_mashm->size));
  ierr = MPI_Comm_rank(p_mashm->comm, &(p_mashm->rank));

  if (p_mashm->rank == 0) {
    p_mashm->isMasterProc = true;
  }
  else {
    p_mashm->isMasterProc = false;
  }

  /* Initialize the intra-node subcommunicator */
  MashmIntraNodeCommInit(&(p_mashm->intraComm),p_mashm->comm);

  /* Now calculate the number of shared memory indices */
  ierr = MPI_Comm_split(p_mashm->comm, p_mashm->intraComm.rank, p_mashm->rank, &(p_mashm->rankComm));
  /* Only the nodal root is participates */
  if (p_mashm->intraComm.rank == 0) {
    ierr = MPI_Comm_size(p_mashm->rankComm, &numSharedMemNodes);
    ierr = MPI_Comm_rank(p_mashm->rankComm, &sharedMemNodeRank);
    /* The number of shared memory nodes */
    p_mashm->numSharedMemNodes = numSharedMemNodes;
    /* The index of each shared memory node */
    p_mashm->sharedMemIndex = sharedMemNodeRank;
    /*
    if (p_mashm->sharedMemIndex == 0) {
      printf("Number of shared memory nodes %d\n", p_mashm->numSharedMemNodes);
    }
    */
  }

  /* Broadcast (to the shared sub comm) the number of shared memory nodes */
  ierr = MPI_Bcast(&(p_mashm->numSharedMemNodes), 1, MPI_INT, 0, p_mashm->intraComm.comm);
  /* Broadcast (to the shared sub comm) the index of each shared memory nodes */
  ierr = MPI_Bcast(&(p_mashm->sharedMemIndex), 1, MPI_INT, 0, p_mashm->intraComm.comm);

  /* Initialize the MashmCommCollection */
  MashmCommCollectionInit(&(p_mashm->commCollection));

  p_mashm->commType = MASHM_COMM_STANDARD;

  p_mashm->isInit = true;
  p_mashm->buffersInit = false;

  p_mashm->minAggScheme = MASHM_MIN_AGG_ROUND_ROBIN;
  //p_mashm->minAggScheme = MASHM_MIN_AGG_ROOT_PROC;

  /* Default is true */
  p_mashm->cacheBlocking = true;
}

double* p_MashmGetBufferPointer(struct MashmPrivate* p_mashm, int msgIndex, MashmSendReceive sendReceive) {
  if (sendReceive == MASHM_SEND) {
    return p_mashm->sendBufferPointers[msgIndex];
  }
  else {
    return p_mashm->recvBufferPointers[msgIndex];
  }
}

void p_MashmRetireBufferPointer(struct MashmPrivate* p_mashm, double** bufPtr) {
  /* TODO: put a debug check here to ensure that the bufPtr is belongs to 
   * sendBufferPointers or recvBufferPointers 
   */
  *bufPtr = NULL;
}

void p_MashmSetupInterNodeComm(struct MashmPrivate* p_mashm) {
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
void p_MashmSetupIntraMsgComm(struct MashmPrivate* p_mashm) {
  int numIntraNodeMsgs = p_mashm->numIntraNodeMsgs;

  /* Allocate the intranode MPI data */
  p_mashm->intraRecvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  p_mashm->intraSendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numIntraNodeMsgs);
  p_mashm->intraRecvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  p_mashm->intraSendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numIntraNodeMsgs);
  
}


void p_MashmSetupIntraSharedComm(struct MashmPrivate* p_mashm) {
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
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in MPI_Win_shared_query\n");
      }
#endif
    }
  }
}


void p_MashmStandardCommBegin(struct MashmPrivate* p_mashm) {
  int ierr;
  int iMsg;
  int msgCounter;

  /* First do the Irecvs */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (! p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;

      ierr = MPI_Irecv(p_MashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE), 
                       p_mashm->commCollection.commArray[iMsg].recvSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->recvRequests[msgCounter]));
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in MPI_Irecv on message %d\n", iMsg);
      }
#endif
    }
  }
  /* Next do the Isends */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (! p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
      ierr = MPI_Isend(p_MashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND), 
                       p_mashm->commCollection.commArray[iMsg].sendSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->sendRequests[msgCounter]));
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in MPI_Isend on message %d\n", iMsg);
      }
#endif
    }
  }
}

/* @brief Finish nodal communication
 *
 * @param in_mash
 *
 * Wait for all internode communication to be completed. Here, we call the MPI_Waitall corresponding to the MPI_Irecv/MPI_Isend calls in MashmInterNodeCommBegin.
 */
inline 
void p_MashmStandardCommEnd(struct MashmPrivate* p_mashm) {
  int ierr;

  int numMsgs;

  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
    numMsgs = p_mashm->numOwnedNodalMsgs;
  }
  else {
    numMsgs = p_mashm->numInterNodeMsgs;
  }
#ifdef MASHM_DEBUG
  /* Allow MPI errors to continue */
  MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
#endif

  ierr = MPI_Waitall(numMsgs, p_mashm->recvRequests, 
                     p_mashm->recvStatuses);
#if MASMM_DEBUG
  char err_buffer[MPI_MAX_ERROR_STRING];
  int errclass, resultlen;
  if (ierr != MPI_SUCCESS) {
    //int resultlen, errclass;
    resultlen;
    MASHM_DEBUG_PRINT("Error reported in MPI_Waitall\n");
    MPI_Error_class(ierr,&errclass);
    int i;
    for (i = 0; i < numMsgs; i++) {
      MPI_Error_string((p_mashm->recvStatuses[i].MPI_ERROR),err_buffer,&resultlen);
      MASHM_DEBUG_PRINT(err_buffer);
    }
  }
#endif

  ierr = MPI_Waitall(numMsgs, p_mashm->sendRequests, 
                     p_mashm->sendStatuses);
#if MASMM_DEBUG
  char err_buffer[MPI_MAX_ERROR_STRING];
  int errclass, resultlen;
  if (ierr != MPI_SUCCESS) {
    //int resultlen, errclass;
    resultlen;
    MASHM_DEBUG_PRINT("Error reported in MPI_Waitall\n");
    MPI_Error_class(ierr,&errclass);
    int i;
    for (i = 0; i < numMsgs; i++) {
      MPI_Error_string((p_mashm->recvStatuses[i].MPI_ERROR),err_buffer,&resultlen);
      MASHM_DEBUG_PRINT(err_buffer);
    }
  }
#endif

  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
    ierr = MPI_Win_fence(MPI_MODE_NOSTORE,p_mashm->sendNodalSharedMemWindow);
    ierr = MPI_Win_fence(MPI_MODE_NOPUT,p_mashm->recvNodalSharedMemWindow);
  }

}

void p_MashmIntraMsgsCommBegin(struct MashmPrivate* p_mashm) {
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
      ierr = MPI_Irecv(p_MashmGetBufferPointer(p_mashm, iMsg, MASHM_RECEIVE), 
                       p_mashm->commCollection.commArray[iMsg].recvSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->intraRecvRequests[msgCounter]));
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in MPI_Irecv on message %d\n", iMsg);
      }
#endif
    }
  }

  /* Next do the Isends */
  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
    if (p_mashm->onNodeMessage[iMsg]) {
      msgCounter = msgCounter + 1;
#if 0
      buf = p_MashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND);
      for (i = 0; i < p_mashm->commCollection.commArray[iMsg].recvSize; i++) {
        if (buf[i] != 666 && buf[i] != 668) {
          printf("Error with send buffer %d rank %d, bufVal %f\n", iMsg, p_mashm->rank, buf[i]);
        }
      }
#endif

      ierr = MPI_Isend(p_MashmGetBufferPointer(p_mashm, iMsg, MASHM_SEND), 
                       p_mashm->commCollection.commArray[iMsg].sendSize, 
                       MPI_DOUBLE,
                       p_mashm->commCollection.commArray[iMsg].pairRank, 
                       1, p_mashm->comm, &(p_mashm->intraSendRequests[msgCounter]));
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in sending message %d\n", iMsg);
      }
#endif
    }
  }
}


void p_MashmIntraMsgsCommEnd(struct MashmPrivate* p_mashm) {
  int ierr;
 
  ierr = MPI_Waitall(p_mashm->numIntraNodeMsgs, p_mashm->intraRecvRequests, 
                     p_mashm->intraRecvStatuses);
  ierr = MPI_Waitall(p_mashm->numIntraNodeMsgs, p_mashm->intraSendRequests, 
                     p_mashm->intraSendStatuses);
}

inline
void p_MashmIntraSharedCommBegin(struct MashmPrivate* p_mashm) {
  int iMsg, msgCounter;
  int ierr;
  int i;
  int sharedBufferOffset;

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

inline
void p_MashmIntraSharedCommEnd(struct MashmPrivate* p_mashm) {
  int ierr = MPI_Win_fence(MPI_MODE_NOSTORE,p_mashm->sendSharedMemWindow);
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
      /* Calculate the number of nodal messages */


      numNodalMessages = p_mashm->numSharedMemNodes;
      floorNumMessagesPerRank = numNodalMessages/p_mashm->intraComm.size;
      modNumMessagesPerRank = numNodalMessages % p_mashm->intraComm.size;

      /* numInterNodeMsgs Invalid at this point */
      p_mashm->numInterNodeMsgs = -1;
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

  p_mashm->bufferSize = bufferSize;
  p_mashm->sharedBufferSize = sharedBufferSize;
  p_mashm->numInterNodePtrs = numInterNodePtrs;
  p_mashm->numIntraNodePtrs = numIntraNodePtrs;
}


/* 
 * Allocates contiguous shared memory data for the send and receive buffers
 */
void p_MashmAllocateSharedMemory(struct MashmPrivate* p_mashm, int bufferSize) {
  int ierr;
  int i;
  int numOrigMsgs;
  int counter, runningOffset;
  int numIntraNodeMsgs;

  MPI_Request *recvRequests, *sendRequests;
  MPI_Status *recvStatuses, *sendStatuses;

  MPI_Info info_noncontig;
  ierr = MPI_Info_create(&info_noncontig);
  ierr = MPI_Info_set(info_noncontig, "alloc_shared_noncontig", "true");

  ierr = MPI_Win_allocate_shared(sizeof(double)*bufferSize, sizeof(double),
                                 info_noncontig, p_mashm->intraComm.comm,
                                 &(p_mashm->p_sharedSendBuffer),&(p_mashm->sendSharedMemWindow));
  if (ierr != 0) {
    printf("Error in MPI_Win_allocate_shared1\n");
  }

  ierr = MPI_Win_allocate_shared(sizeof(double)*bufferSize, sizeof(double),
                                 info_noncontig, p_mashm->intraComm.comm,
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

void p_MashmAllocateSharedMemoryMinAgg(struct MashmPrivate* p_mashm) {

  int ierr;
  int allocateSize;
  int iMsg, msgCounter;
  MPI_Aint otherWinSize;
  int otherDispUnit;

  MPI_Info info_noncontig;
  ierr = MPI_Info_create(&info_noncontig);
  ierr = MPI_Info_set(info_noncontig, "alloc_shared_noncontig", "true");

  /* Force padding on cache lines */
  allocateSize = p_mashm->nodalSharedBufferSize+7;
  /* Some MPI implementations don't like size zero arrays */
  if (allocateSize == 0) {
    /* Make this a full cache-line */
    allocateSize = 8;
  }

  ierr = MPI_Win_allocate_shared(sizeof(double)*allocateSize, sizeof(double),
                                 info_noncontig, p_mashm->intraComm.comm,
                                 &(p_mashm->p_sendNodalSharedBuffer),&(p_mashm->sendNodalSharedMemWindow));
  if (ierr != 0) {
    printf("Error in MPI_Win_allocate_shared1\n");
  }

  ierr = MPI_Win_allocate_shared(sizeof(double)*allocateSize, sizeof(double),
                                 info_noncontig, p_mashm->intraComm.comm,
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

void p_MashmCalculateNodalMsgSchedule(struct MashmPrivate* p_mashm) {
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
#if 0
    for (iRank = 0; iRank < p_mashm->numSharedMemNodes; iRank++) {

      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);

      if (p_mashm->sharedMemIndex == iRank) {
        for (i = 0; i < sumNumMsgs; i++) {
          printf("Rank %d: msg %i: srcRank %d, nodeIndex %d, destRank %d, size %d\n",
                 p_mashm->rank, i, commArray[i].srcSharedMemRank, commArray[i].destNodeIndex,
                 commArray[i].destGlobalRank, commArray[i].msgSize);
        }
      }
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      fflush(stdout);
      sleep(1);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
      ierr = MPI_Barrier(p_mashm->rankComm);
    }
#endif

    commArrayOffset = sumNumMsgs;
    /* Now advance the commArray to the first non-self node */
    for (i = 0; i < sumNumMsgs; i++) {
      if (commArray[i].destNodeIndex != p_mashm->sharedMemIndex) {
        commArrayOffset = i;
        break;
      }
    }
    p_mashm->numNodalSubMsgs = sumNumMsgs - commArrayOffset;

    /* Calculate the msg sizes */
    p_mashm->allMsgSizes = 0;
    p_mashm->intraMsgSizes = 0;
    p_mashm->interMsgSizes = 0;
    p_mashm->minIntraMsgSize = 100000000;
    p_mashm->maxIntraMsgSize = 0;
    p_mashm->minInterMsgSize = 100000000;
    p_mashm->maxInterMsgSize = 0;
    p_mashm->minNodalMsgSize = 100000000;
    p_mashm->maxNodalMsgSize = 0;
    p_mashm->sumNodalMsgSize = 0;
    for (i = 0; i < commArrayOffset; i++) {
      p_mashm->allMsgSizes += commArray[i].msgSize;
      p_mashm->intraMsgSizes += commArray[i].msgSize;
      p_mashm->minIntraMsgSize = MIN(commArray[i].msgSize,p_mashm->minIntraMsgSize);
      p_mashm->maxIntraMsgSize = MAX(commArray[i].msgSize,p_mashm->maxIntraMsgSize);
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
      p_mashm->interMsgSizes += commArray[i].msgSize;
      p_mashm->allMsgSizes += commArray[i].msgSize;
      p_mashm->minInterMsgSize=MIN(commArray[i].msgSize,p_mashm->minInterMsgSize);
      p_mashm->maxInterMsgSize=MAX(commArray[i].msgSize,p_mashm->maxInterMsgSize);
    }
    p_mashm->numNodalMsgs = nodeCounter;

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


    p_mashm->uniqueNodeIndices = (int*) malloc(sizeof(int)*(p_mashm->numCommNodes));

    nodeCounter = 0;
    p_mashm->uniqueNodeIndices[nodeCounter] = commArray[0].destNodeIndex;
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

    p_mashm->nodalMsgSendSizes = (int*) malloc(sizeof(int)*p_mashm->numNodalMsgs);
    p_mashm->nodalMsgRecvSizes = (int*) malloc(sizeof(int)*p_mashm->numNodalMsgs);

    nodeCounter = 0;
    if (p_mashm->numNodalSubMsgs > 0) {
      msgOffsets[0] = 0;
    }

    tmpMsgSize = 0;
    msgOffsets[0] = 0;
    for (i = 0; i < p_mashm->numNodalSubMsgs-1; i++) {
      iOffset = i+commArrayOffset;
      tmpMsgSize += commArray[iOffset].msgSize;
      if (commArray[iOffset+1].destNodeIndex != commArray[iOffset].destNodeIndex) {
        p_mashm->nodalMsgSendSizes[nodeCounter] = tmpMsgSize;
        nodeCounter += 1;
        tmpMsgSize = 0;
      }
      else if (p_mashm->cacheBlocking && 
               commArray[iOffset+1].srcSharedMemRank != commArray[iOffset].srcSharedMemRank) {
        /* Add 7 doubles to cache block writing to shared memory */
        tmpMsgSize += 7;
      }
      msgOffsets[i+1] = tmpMsgSize;
    }
    tmpMsgSize += commArray[p_mashm->numNodalSubMsgs - 1 + commArrayOffset].msgSize;
    p_mashm->nodalMsgSendSizes[nodeCounter] = tmpMsgSize;

    if (p_mashm->cacheBlocking) {
      /* When cache blocking we need to exchange the size of the 
       *   nodal messages since they are no longer the same */
      recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*p_mashm->numNodalMsgs);
      sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*p_mashm->numNodalMsgs);
      recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*p_mashm->numNodalMsgs);
      sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*p_mashm->numNodalMsgs);

      /* Usual point to point communication */
      for (i = 0; i < p_mashm->numNodalMsgs; i++) {
        //printf("Rank %d receiving message from rank %d\n", p_mashm->sharedMemIndex, p_mashm->uniqueNodeIndices[i+1]);
        ierr = MPI_Irecv(&(p_mashm->nodalMsgRecvSizes[i]), 1, MPI_INT, p_mashm->uniqueNodeIndices[i+1], 0, p_mashm->rankComm, &recvRequests[i]);
      }
      for (i = 0; i < p_mashm->numNodalMsgs; i++) {
        //printf("Rank %d sending rank %d message\n", p_mashm->sharedMemIndex, p_mashm->uniqueNodeIndices[i+1]);
        ierr = MPI_Isend(&(p_mashm->nodalMsgSendSizes[i]), 1, MPI_INT, p_mashm->uniqueNodeIndices[i+1], 0, p_mashm->rankComm, &sendRequests[i]);
      }

      ierr = MPI_Waitall(p_mashm->numNodalMsgs,recvRequests,recvStatuses); 
      ierr = MPI_Waitall(p_mashm->numNodalMsgs,sendRequests,sendStatuses); 
      free(recvRequests);
      free(sendRequests);
      free(recvStatuses);
      free(sendStatuses);
    }
    else {
      /* Receive length of nodal messages is the same as the send lengths */
      for (i = 0; i < p_mashm->numNodalMsgs; i++) {
        p_mashm->nodalMsgRecvSizes[i] = p_mashm->nodalMsgSendSizes[i];
      }
    }

    for (i = 0; i < p_mashm->numNodalMsgs; i++) {
      p_mashm->minNodalMsgSize = MIN(p_mashm->nodalMsgSendSizes[i],p_mashm->minNodalMsgSize);
      p_mashm->maxNodalMsgSize = MAX(p_mashm->nodalMsgSendSizes[i],p_mashm->maxNodalMsgSize);
      p_mashm->sumNodalMsgSize = p_mashm->sumNodalMsgSize + p_mashm->nodalMsgSendSizes[i];
    }

    for (i = 0; i < sumNumMsgs; i++) {
      allSrcSharedMemRanks[i] = commArray[i].srcSharedMemRank;
      allDestGlobalRanks[i] = commArray[i].destGlobalRank;
      allMsgNodeIndices[i] = commArray[i].destNodeIndex;
    }

  }

  /* Now scatter the sorted list of srcSharedRanks and destGlobalRanks */
  ierr = MPI_Bcast(&sumNumMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&commArrayOffset, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&p_mashm->numNodalSubMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&p_mashm->numNodalMsgs, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(&p_mashm->numCommNodes, 1, MPI_INT, 0, p_mashm->intraComm.comm);
  p_mashm->sumNumMsgs = sumNumMsgs;

  if (!p_mashm->intraComm.isMasterProc) {
    allSrcSharedMemRanks = (int*) malloc(sizeof(int)*sumNumMsgs);
    allDestGlobalRanks = (int*) malloc(sizeof(int)*sumNumMsgs);
    allMsgNodeIndices = (int*) malloc(sizeof(int)*sumNumMsgs);
    msgOffsets = (int*) malloc(sizeof(int)*p_mashm->numNodalSubMsgs);
    p_mashm->nodalMsgSendSizes = (int*) malloc(sizeof(int)*p_mashm->numNodalMsgs);
    p_mashm->nodalMsgRecvSizes = (int*) malloc(sizeof(int)*p_mashm->numNodalMsgs);
    p_mashm->uniqueNodeIndices = (int*) malloc(sizeof(int)*(p_mashm->numCommNodes));
  }

  /* Broadcast the msgOffsets and the nodalMsgsSizes */
  ierr = MPI_Bcast(allSrcSharedMemRanks, sumNumMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(allDestGlobalRanks, sumNumMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(allMsgNodeIndices, sumNumMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(msgOffsets, p_mashm->numNodalSubMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(p_mashm->nodalMsgSendSizes, p_mashm->numNodalMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
  ierr = MPI_Bcast(p_mashm->nodalMsgRecvSizes, p_mashm->numNodalMsgs, MPI_INT, 0, p_mashm->intraComm.comm);
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

void p_MashmSetupAggType(struct MashmPrivate* p_mashm) {
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
            p_mashm->nodalSharedBufferSize += p_mashm->nodalMsgSendSizes[i];
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
          p_mashm->nodalSharedBufferSize += p_mashm->nodalMsgSendSizes[i];
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
      nodalMsgOffset += p_mashm->nodalMsgSendSizes[i];
    }
  }

  
  /* Now we need to exchange the rank of the nodal message owner with the rank of the receiver message owner */
  p_mashm->nodalRecvRank = (int*) malloc(sizeof(int)*numNodalMsgs);
  /* For this we need the nodal Comm again */

  /* Only the nodal root is participates */
  if (p_mashm->intraComm.rank == 0) {

    tmpRecvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNodalMsgs);
    tmpSendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNodalMsgs);

    tmpRecvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNodalMsgs);
    tmpSendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNodalMsgs);

    /* Buffer to write the messaged data */
    int* nodalSendData = (int*) malloc(sizeof(int)*numNodalMsgs);

    /* The rank is the node index ... */
    /* Exchange the ranks with the */
    for (i = 0; i < numNodalMsgs; i++) {
      ierr = MPI_Irecv(&(p_mashm->nodalRecvRank[i]), 1, MPI_INT,
                       p_mashm->uniqueNodeIndices[i+1], 1, p_mashm->rankComm, &(tmpRecvRequests[i])); 
      if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Irecv %d\n", ierr);
      }
    }
    for (i = 0; i < numNodalMsgs; i++) {
      sharedRankMsgOwner = p_mashm->nodalMsgOwner[i];
      nodalSendData[i] = p_mashm->intraComm.parentRanksOnNode[sharedRankMsgOwner];
      ierr = MPI_Isend(&(nodalSendData[i]), 1, MPI_INT,
                       p_mashm->uniqueNodeIndices[i+1], 1, p_mashm->rankComm, &(tmpSendRequests[i])); 
      if (ierr != MPI_SUCCESS) {
        printf("Error in MPI_Isend %d\n", ierr);
      }
    }
    ierr = MPI_Waitall(numNodalMsgs, tmpRecvRequests, tmpRecvStatuses);
    ierr = MPI_Waitall(numNodalMsgs, tmpSendRequests, tmpSendStatuses);
    free(tmpRecvRequests);
    free(tmpSendRequests);
    free(tmpRecvStatuses);
    free(tmpSendStatuses);

    free(nodalSendData);

  }

  /* Broadcast this information to all of the processes within the node */
  ierr = MPI_Bcast(p_mashm->nodalRecvRank, numNodalMsgs, MPI_INT, 0, p_mashm->intraComm.comm);

}


void p_MashmCalcMsgIndicesMinAgg(struct MashmPrivate* p_mashm) {
  int sendNodalSharedBufferOffset;
  int nodalMsgOffset;
  int sendMsgOffset, recvMsgOffset;
  int nodalMsgDest;
  int nodalMsgIndex;
  int i, iMsg;
  int sharedBufferOffset;

  int ierr, iNode, iRank;
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
#if 0  
  for (iNode = 0; iNode < p_mashm->numSharedMemNodes; iNode++) {
    if (iNode == p_mashm->sharedMemIndex) {
      if (p_mashm->intraComm.rank == 0) {
        printf("Info for node %d\n", p_mashm->sharedMemIndex);
      }
      for (iRank = 0; iRank < p_mashm->intraComm.size; iRank++) {
        if (iRank == p_mashm->intraComm.rank) {
          for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
            if ( ! p_mashm->onNodeMessage[iMsg]) {
            printf("  off node source, dest, sendOffset, recvOffset, size %d, %d, %d, %d, %d\n", 
              p_mashm->rank,  
              (p_mashm->commCollection.commArray[iMsg]).pairRank, 
              p_mashm->sendAggMsgOffsets[iMsg],
              p_mashm->recvAggMsgOffsets[iMsg],
              (p_mashm->commCollection.commArray[iMsg]).sendSize);
            }
            else {
            printf("  on node source, dest %d, %d\n", 
              p_mashm->rank,  
              (p_mashm->commCollection.commArray[iMsg]).pairRank);
            }
          }
          ierr = MPI_Barrier(p_mashm->intraComm.comm);
          ierr = MPI_Barrier(p_mashm->intraComm.comm);
          fflush(stdout);
          ierr = MPI_Barrier(p_mashm->intraComm.comm);
          ierr = MPI_Barrier(p_mashm->intraComm.comm);
        }
      }
    }
    ierr = MPI_Barrier(p_mashm->comm);
    ierr = MPI_Barrier(p_mashm->comm);
    fflush(stdout);
    ierr = MPI_Barrier(p_mashm->comm);
    ierr = MPI_Barrier(p_mashm->comm);
  }
#endif
}

inline
void p_MashmMinAggCommBegin(struct MashmPrivate* p_mashm) {
  int ierr;
  int iMsg;
  //int numMsgs = p_mashm->commCollection.commArraySize;
  //int i;
  //double* buf;
  int msgCounter, sharedRankMsgOwner;


  /* Post the receives for the next cycle */

  /* MPI_Win_fence */
  ierr = MPI_Win_fence(MPI_MODE_NOSTORE,p_mashm->recvNodalSharedMemWindow);

  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->numNodalMsgs; iMsg++) {
    sharedRankMsgOwner = p_mashm->nodalMsgOwner[iMsg];
    if (sharedRankMsgOwner == p_mashm->intraComm.rank) {
      msgCounter += 1;
      ierr = MPI_Irecv(&(p_mashm->p_recvNodalSharedBuffer[p_mashm->nodalOffsets[iMsg]]),
                       p_mashm->nodalMsgRecvSizes[iMsg], MPI_DOUBLE,
                       p_mashm->nodalRecvRank[iMsg],
                       1, p_mashm->comm, &(p_mashm->recvRequests[msgCounter]));
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in MPI_Irecv message%d\n", iMsg);
      }
#endif
    }
  }

  /* MPI_Win_fence */
  ierr = MPI_Win_fence(MPI_MODE_NOPUT,p_mashm->sendNodalSharedMemWindow);

  msgCounter = -1;
  for (iMsg = 0; iMsg < p_mashm->numNodalMsgs; iMsg++) {
    sharedRankMsgOwner = p_mashm->nodalMsgOwner[iMsg];
    if (sharedRankMsgOwner == p_mashm->intraComm.rank) {
      msgCounter += 1;
      ierr = MPI_Isend(&(p_mashm->p_sendNodalSharedBuffer[p_mashm->nodalOffsets[iMsg]]),
                       p_mashm->nodalMsgSendSizes[iMsg], MPI_DOUBLE,
                       p_mashm->nodalRecvRank[iMsg],
                       1, p_mashm->comm, &(p_mashm->sendRequests[msgCounter]));
#ifdef MASHM_DEBUG
      if (ierr != 0) {
        MASHM_DEBUG_PRINT("Error in MPI_Isend message%d\n", iMsg);
      }
#endif
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


void p_MashmFinish(struct MashmPrivate* p_mashm) {
  int regularBufferOffset, sharedBufferOffset;
  int iMsg, msgDestRank;
  int i;
  int numMessages;
  /* Allocate data for the number of messages */
  numMessages = p_mashm->commCollection.commArraySize;

  p_mashm->numOrigMessages = numMessages;

  p_mashm->sendBufferPointers = (double**) malloc(sizeof(double*)*p_mashm->numOrigMessages);
  p_mashm->recvBufferPointers = (double**) malloc(sizeof(double*)*p_mashm->numOrigMessages);
  p_mashm->onNodeMessage = (MashmBool*) malloc(sizeof(MashmBool)*p_mashm->numOrigMessages);
  p_mashm->pairRanks = (int*) malloc(sizeof(int)*numMessages);
  p_mashm->pairSharedRanks = (int*) malloc(sizeof(int)*numMessages);

  /* */
  for (i = 0; i < p_mashm->numOrigMessages; i++) {
    p_mashm->pairRanks[i] = (p_mashm->commCollection.commArray[i]).pairRank;
    p_mashm->pairSharedRanks[i] = MashmIntraNodeCommGetSharedRank(p_mashm->intraComm, p_mashm->pairRanks[i]);
  }

  /* Determine the number of intranode and internode messages and calculate the sizes */
  /* TODO: these two methods should be condensed into one -
   *       they are doing similar things 
   */
  p_MashmCalcNumMpiMsgs(p_mashm);
  p_MashmCalcMsgBufferSize(p_mashm);
  /*
  printf("Buf size %d, shared buf size %d\n", p_mashm->bufferSize,p_mashm->sharedBufferSize);
  printf("Rank %d, intra rank %d sends %d MPI messages\n", p_mashm->rank, p_mashm->intraComm.rank, p_mashm->numInterNodeMsgs);
  */
  
  /* Allocate MPI (non-shared) buffer */
  p_mashm->p_regularSendBuffer = (double*) malloc(sizeof(double)*p_mashm->bufferSize);
  p_mashm->p_regularRecvBuffer = (double*) malloc(sizeof(double)*p_mashm->bufferSize);

  /* Allocate MPI shared memory */
  if (p_mashm->commType == MASHM_COMM_INTRA_SHARED ||
      p_mashm->commType == MASHM_COMM_MIN_AGG) {
    p_MashmAllocateSharedMemory(p_mashm, p_mashm->sharedBufferSize);
  }
  p_mashm->buffersInit = true;

  /* Set pointers to the buffer and shared buffer */
  if (p_mashm->commType == MASHM_COMM_MIN_AGG) {
    p_MashmCalculateNodalMsgSchedule(p_mashm);
    p_MashmSetupAggType(p_mashm);
    p_MashmAllocateSharedMemoryMinAgg(p_mashm);
    p_MashmCalcMsgIndicesMinAgg(p_mashm);

  }
  else {
    regularBufferOffset = 0;
    sharedBufferOffset = 0;
    for (iMsg = 0; iMsg < p_mashm->commCollection.commArraySize; iMsg++) {
      msgDestRank = (p_mashm->commCollection.commArray[iMsg]).pairRank;
      if ( (p_mashm->commType == MASHM_COMM_INTRA_SHARED ||
            p_mashm->commType == MASHM_COMM_MIN_AGG) 
          && p_MashmIsIntraNodeRank(p_mashm, msgDestRank)) {
        /* Set the shared memory pointers */
        p_mashm->sendBufferPointers[iMsg] = &(p_mashm->p_sharedSendBuffer[sharedBufferOffset]);
        p_mashm->recvBufferPointers[iMsg] = &(p_mashm->p_sharedRecvBuffer[sharedBufferOffset]);
        sharedBufferOffset = sharedBufferOffset + (p_mashm->commCollection.commArray[iMsg]).sendSize;
        p_mashm->onNodeMessage[iMsg] = true;
      }
      else {
        /* Set the pointers to regular memory */
        p_mashm->sendBufferPointers[iMsg] = &(p_mashm->p_regularSendBuffer[regularBufferOffset]);
        p_mashm->recvBufferPointers[iMsg] = &(p_mashm->p_regularRecvBuffer[regularBufferOffset]);
        regularBufferOffset = regularBufferOffset + (p_mashm->commCollection.commArray[iMsg]).sendSize;
        if (p_mashm->commType == MASHM_COMM_INTRA_MSG &&
            p_MashmIsIntraNodeRank(p_mashm, msgDestRank)) {
          /* On node MPI_Isend/MPI_Irecv */
          p_mashm->onNodeMessage[iMsg] = true;
        }
        else {
          p_mashm->onNodeMessage[iMsg] = false;
        }
      }
    }
  }

  /* All communication types need to setup the data needed for MPI_Irecv/MPI_Isend */
  p_MashmSetupInterNodeComm(p_mashm);

  /* Perform some initialization for each comm method and set pointers to proper routines */
  switch (p_mashm->commType) {
    case MASHM_COMM_STANDARD:
      /*
      printf("Using standard communication scheme\n");
      */
      p_mashm->p_interNodeCommBegin = p_MashmStandardCommBegin;
      p_mashm->p_interNodeCommEnd = p_MashmStandardCommEnd;

      p_mashm->p_intraNodeCommBegin = p_MashmNullFunction;
      p_mashm->p_intraNodeCommEnd = p_MashmNullFunction;

      break;

    case MASHM_COMM_INTRA_MSG:
      /*
      printf("Using INTRA_MSG communication scheme\n");
      */
      p_MashmSetupIntraMsgComm(p_mashm);

      p_mashm->p_interNodeCommBegin = p_MashmStandardCommBegin;
      p_mashm->p_interNodeCommEnd = p_MashmStandardCommEnd;

      p_mashm->p_intraNodeCommBegin = p_MashmIntraMsgsCommBegin;
      p_mashm->p_intraNodeCommEnd = p_MashmIntraMsgsCommEnd;

      break;

    case MASHM_COMM_INTRA_SHARED:
      /*
      printf("Using INTRA_SHARED communication scheme\n");
      */
      /* Calculate shared memory indices etc. */
      p_MashmSetupIntraSharedComm(p_mashm);

      p_mashm->p_interNodeCommBegin = p_MashmStandardCommBegin;
      p_mashm->p_interNodeCommEnd = p_MashmStandardCommEnd;

      p_mashm->p_intraNodeCommBegin = p_MashmIntraSharedCommBegin;
      p_mashm->p_intraNodeCommEnd = p_MashmIntraSharedCommEnd;

      break;

    case MASHM_COMM_MIN_AGG:
      /*
      printf("Using MIN_AGG communication scheme\n");
      */
      p_MashmSetupIntraSharedComm(p_mashm);

      p_mashm->p_interNodeCommBegin = p_MashmMinAggCommBegin;
      p_mashm->p_interNodeCommEnd = p_MashmStandardCommEnd;

      p_mashm->p_intraNodeCommBegin = p_MashmIntraSharedCommBegin;
      p_mashm->p_intraNodeCommEnd = p_MashmIntraSharedCommEnd;

      break;

  }
#if 0
  if (p_mashm->intraComm.isMasterProc) {
    printf("Node %d has %d nodal messages.\n", p_mashm->sharedMemIndex, p_mashm->numNodalMsgs);
  }
  // Test if there are more messages than intra-comm ranks
  if (p_mashm->numNodalMsgs > p_mashm->intraComm.size) {
    if (p_mashm->intraComm.rank == 0) {
      printf("Node %d has more messages than ranks\n", p_mashm->sharedMemIndex);
    }
  }
  //printf("Rank %d, owns %d messages\n", p_mashm->rank, p_mashm->numOwnedNodalMsgs);
  //fflush(stdout);
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Barrier(MPI_COMM_WORLD);
  //p_MashmPrintInterNodeMessages(p_mashm);

  p_MashmPrintMessageInformation(p_mashm);
#endif
}

void p_MashmPrintMessageInformation(struct MashmPrivate* p_mashm) {

  /* Determine the number of total messages, intranode message, and internode messages
   *   Determine the min, max, and avg. msg size
   *   Determine the min, max, and avg. nodal msg. size */

  /* Reduce all of the data */ 
  int numTotalMsgs;
  int numIntraNodeMsgs, numInterNodeMsgs;
  int ierr;

  if (p_mashm->intraComm.isMasterProc) {
    /* Note that all data is for a node */
    int intraCommNumTotalMsgs = p_mashm->sumNumMsgs;
    int intraCommNumInterNodeMsgs = p_mashm->numNodalSubMsgs;
    int intraCommNumIntraNodeMsgs = intraCommNumTotalMsgs - intraCommNumInterNodeMsgs;
    int intraCommNumNodalMsgs = p_mashm->numNodalMsgs;
    int intraSizeMin, intraSizeMax, interSizeMin, interSizeMax;

    int numTotalMsgs, numInterNodeMsgs, numIntraNodeMsgs, numNodalMsgs;
    int totalAllMsgSizes, totalIntraMsgSizes, totalInterMsgSizes;
    int minNodalMsgSize, maxNodalMsgSize, sumNodalMsgSize;

    // Reduce to master rank of rankComm
    ierr = MPI_Reduce(&intraCommNumTotalMsgs, &numTotalMsgs, 1, MPI_INT, 
                      MPI_SUM, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&intraCommNumInterNodeMsgs, &numInterNodeMsgs, 1, MPI_INT, 
                      MPI_SUM, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&intraCommNumIntraNodeMsgs, &numIntraNodeMsgs, 1, MPI_INT, 
                      MPI_SUM, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&intraCommNumNodalMsgs, &numNodalMsgs, 1, MPI_INT, 
                      MPI_SUM, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->allMsgSizes), &(totalAllMsgSizes), 1, MPI_INT,
                      MPI_SUM, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->intraMsgSizes), &(totalIntraMsgSizes), 1, MPI_INT,
                      MPI_SUM, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->interMsgSizes), &(totalInterMsgSizes), 1, MPI_INT,
                      MPI_SUM, 0, p_mashm->rankComm);


    ierr = MPI_Reduce(&(p_mashm->minIntraMsgSize), &(intraSizeMin), 1, MPI_INT,
                      MPI_MIN, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->maxIntraMsgSize), &(intraSizeMax), 1, MPI_INT,
                      MPI_MAX, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->minInterMsgSize), &(interSizeMin), 1, MPI_INT,
                      MPI_MIN, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->maxInterMsgSize), &(interSizeMax), 1, MPI_INT,
                      MPI_MAX, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->minNodalMsgSize), &(minNodalMsgSize), 1, MPI_INT,
                      MPI_MIN, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->maxNodalMsgSize), &(maxNodalMsgSize), 1, MPI_INT,
                      MPI_MAX, 0, p_mashm->rankComm);
    ierr = MPI_Reduce(&(p_mashm->sumNodalMsgSize), &(sumNodalMsgSize), 1, MPI_INT,
                      MPI_SUM, 0, p_mashm->rankComm);


    // Allocate data for sizes
    if (p_mashm->rank == 0) {
      /* Note that the MASHM: is there to allow parsing this data out of stdout */
      printf("MASHM: Total number of messages %d\n", numTotalMsgs);
      printf("MASHM: Size of all messages %d B\n", totalAllMsgSizes*8);
      printf("MASHM: Number internode messages %d\n", numInterNodeMsgs);
      printf("MASHM: Size of internode messages %d B\n", totalInterMsgSizes*8);
      printf("MASHM:   Min, Max, Avg %d, %f, %d B\n", interSizeMin*8, ((double)totalInterMsgSizes*8/numInterNodeMsgs), interSizeMax*8);
      printf("MASHM: Number intranode messages %d\n", numIntraNodeMsgs);
      printf("MASHM: Size of intranode messages %d B\n", totalIntraMsgSizes*8);
      printf("MASHM:   Min, Max, Avg %d, %f, %d B\n", intraSizeMin*8, ((double) totalIntraMsgSizes)/numIntraNodeMsgs*8, intraSizeMax*8);
      printf("MASHM: Number of nodal messages %d\n", numNodalMsgs);
      printf("MASHM:   Min, Max, Avg %d, %f, %d\n", minNodalMsgSize, ((double) totalInterMsgSizes)/numNodalMsgs, maxNodalMsgSize);
    }
  }
 

}


void p_Destroy(struct MashmPrivate* p_mashm) {
  int ierr;

  /* Destroy the MashmCommCollection 
   * TODO: this should be destroyed in the finish method */
  MashmCommCollectionDestroy(&(p_mashm->commCollection));

  /* Destroy the intra-node subcommunicator */
  MashmIntraNodeCommDestroy(&(p_mashm->intraComm));

  /* Destroy the rank communicator */
  ierr = MPI_Comm_free(&(p_mashm->rankComm));

  if (p_mashm->buffersInit) {
    free(p_mashm->sendBufferPointers);
    free(p_mashm->recvBufferPointers);
    free(p_mashm->p_regularSendBuffer);
    free(p_mashm->p_regularRecvBuffer);
    free(p_mashm->onNodeMessage);
    free(p_mashm->pairRanks);
    free(p_mashm->pairSharedRanks);
  }

  /* Deallocate the shared memory window created with the
   *   call to MPI_Win_allocate_shared
   *   Note that this also frees the underlying shared memory */
  if (p_mashm->commType == MASHM_COMM_INTRA_SHARED ||
      p_mashm->commType == MASHM_COMM_MIN_AGG) {
    ierr = MPI_Win_free(&(p_mashm->sendSharedMemWindow));
    ierr = MPI_Win_free(&(p_mashm->recvSharedMemWindow));

  }

}

void p_MashmPrintInfo(const struct MashmPrivate* p_mashm) {
  int iNode;

  if (p_mashm->isMasterProc) {
    printf("Number of shared memory nodes %d\n", p_mashm->numSharedMemNodes);
  }

  for (iNode = 0; iNode < p_mashm->numSharedMemNodes; iNode++) {
    if (p_mashm->isMasterProc) {
      printf("  Node %d\n", iNode);
    }
    if (p_mashm->sharedMemIndex == iNode) {
      MashmIntraNodeCommPrintInfo(p_mashm->intraComm);
    }
  }
}

double* p_MashmGetBufferPointerForDest(const struct MashmPrivate* p_mashm, int destRank, MashmSendReceive sendReceive) {
  int iRank;

  for (iRank = 0; iRank < p_mashm->intraComm.size; iRank++) {
    if (destRank == p_mashm->intraComm.parentRanksOnNode[iRank]) {
      if (sendReceive == MASHM_SEND) {
        return p_mashm->sendBufferPointers[iRank];
      }
      else {
        return p_mashm->recvBufferPointers[iRank];
      }
    }
  }
  return NULL;
}

void p_MashmNullFunction(struct MashmPrivate* p_mashm) {}


/* Get communicator */
MPI_Comm p_MashmGetComm(const struct MashmPrivate* p_mashm) {
  return p_mashm->comm;
}

/* Get MPI communicator size */
int p_MashmGetSize(const struct MashmPrivate* p_mashm) {
  return p_mashm->size;
}

/* Get MPI communicator rank */
int p_MashmGetRank(const struct MashmPrivate* p_mashm) {
  return p_mashm->rank;
}

void p_MashmSetNumComms(struct MashmPrivate* p_mashm, int numComms) {
  MashmCommCollectionSetSize(&(p_mashm->commCollection), numComms);
}

void p_MashmSetComm(struct MashmPrivate* p_mashm, int commIndex, int pairRank, int msgSize) {
  MashmCommCollectionSetComm(&(p_mashm->commCollection), commIndex, pairRank, msgSize, msgSize);
}

/* Get the rank of communication */
int p_MashmGetCommRank(struct MashmPrivate* p_mashm, int commIndex) {
  return MashmCommCollectionGetCommRank(&(p_mashm->commCollection), commIndex);
}

/* Get the size of communication */
int p_MashmGetCommSize(struct MashmPrivate* p_mashm, int commIndex) {
  return MashmCommCollectionGetCommSize(&(p_mashm->commCollection), commIndex);
}

void p_MashmPrintCommCollection(const struct MashmPrivate* p_mashm) {
  int i;
  for (i = 0; i < p_mashm->size; i++) {
    if (i == p_mashm->rank) {
      printf("Rank %d has communication:\n", p_mashm->rank);
      MashmCommCollectionPrint(p_mashm->commCollection);
    }
  }
}

void p_MashmSetCommMethod(struct MashmPrivate* p_mashm, MashmCommType commType) {
  p_mashm->commType = commType;
}

MashmCommType p_MashmGetCommMethod(const struct MashmPrivate* p_mashm) {
  return p_mashm->commType;
}

inline
MashmBool p_MashmIsMsgOnNode(const struct MashmPrivate* p_mashm, int msgIndex) {
  return p_mashm->onNodeMessage[msgIndex];
}

inline 
void p_MashmInterNodeCommBegin(struct MashmPrivate* p_mashm) {
  p_mashm->p_interNodeCommBegin(p_mashm);
}

inline 
void p_MashmIntraNodeCommBegin(struct MashmPrivate* p_mashm) {
  p_mashm->p_intraNodeCommBegin(p_mashm);
}

inline 
void p_MashmIntraNodeCommEnd(struct MashmPrivate* p_mashm) {
  p_mashm->p_intraNodeCommEnd(p_mashm);
}

inline 
void p_MashmInterNodeCommEnd(struct MashmPrivate* p_mashm) {
  p_mashm->p_interNodeCommEnd(p_mashm);
}

void p_MashmPrintInterNodeMessages(struct MashmPrivate* p_mashm) {
  int iNode, i, iMsg, msgCounter, ierr;
  for (iNode = 0; iNode < p_mashm->numSharedMemNodes; iNode++) {
    if (iNode == p_mashm->sharedMemIndex) {
      if (p_mashm->intraComm.rank == 0) {
        printf("Nodal messages for node %d\n", p_mashm->sharedMemIndex);
      }
      //printf("  rank, bufSize %d, %d\n", p_mashm->rank, p_mashm->nodalSharedBufferSize);
      msgCounter = -1;
      for (iMsg = 0; iMsg < p_mashm->numNodalMsgs; iMsg++ ) {
        if (p_mashm->nodalMsgOwner[iMsg] == p_mashm->intraComm.rank) {
          msgCounter = msgCounter + 1;
          printf("  source, dest, size, offsets, bufsize %d, %d, %d, %d, %d\n", p_mashm->rank, p_mashm->nodalRecvRank[iMsg], p_mashm->nodalMsgSendSizes[iMsg], p_mashm->nodalOffsets[iMsg], p_mashm->nodalSharedBufferSize);
        }
        ierr = MPI_Barrier(p_mashm->intraComm.comm);
        ierr = MPI_Barrier(p_mashm->intraComm.comm);
        ierr = MPI_Barrier(p_mashm->intraComm.comm);
        ierr = MPI_Barrier(p_mashm->intraComm.comm);
      }
    }
    ierr = MPI_Barrier(p_mashm->comm);
    ierr = MPI_Barrier(p_mashm->comm);
    ierr = MPI_Barrier(p_mashm->comm);
    ierr = MPI_Barrier(p_mashm->comm);
  }
}

void p_MashmSetCacheBlocking(struct MashmPrivate* p_mashm, MashmBool blockCache) {
  p_mashm->cacheBlocking = blockCache;
}
