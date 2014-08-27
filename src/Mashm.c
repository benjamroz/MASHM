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
  p_Init(in_mashm->p);

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
    in_mashm.p->pairSharedRanks[i] = MashmIntraNodeCommGetSharedRank(in_mashm.p->intraComm, in_mashm.p->pairRanks[i]);
  }

  /* Determine the number of intranode and internode messages and calculate the sizes */
  /* TODO: these two methods should be condensed into one -
   *       they are doing similar things 
   */
  p_MashmCalcNumMpiMsgs(in_mashm.p);
  p_MashmCalcMsgBufferSize(in_mashm.p);
  /*
  printf("Buf size %d, shared buf size %d\n", in_mashm.p->bufferSize,in_mashm.p->sharedBufferSize);
  printf("Rank %d, intra rank %d sends %d MPI messages\n", in_mashm.p->rank, in_mashm.p->intraComm.rank, in_mashm.p->numInterNodeMsgs);
  */
  
  /* Allocate MPI (non-shared) buffer */
  in_mashm.p->p_regularSendBuffer = (double*) malloc(sizeof(double)*in_mashm.p->bufferSize);
  in_mashm.p->p_regularRecvBuffer = (double*) malloc(sizeof(double)*in_mashm.p->bufferSize);

  /* Allocate MPI shared memory */
  if (in_mashm.p->commType == MASHM_COMM_INTRA_SHARED ||
      in_mashm.p->commType == MASHM_COMM_MIN_AGG) {
    p_mashmAllocateSharedMemory(in_mashm.p, in_mashm.p->sharedBufferSize);
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

void MashmRetireBufferPointer(Mashm in_mashm, double** bufPtr) {
  /* Call the private routine */
  return p_mashmRetireBufferPointer(in_mashm.p, bufPtr);
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

