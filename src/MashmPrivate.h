#ifndef MASHM_PRIVATE_H
#define MASHM_PRIVATE_H

#include "mpi.h"

#include "MashmBool.h"
#include "MashmDefs.h"
#include "MashmIntraNodeComm.h"
#include "MashmCommCycle.h"

struct MashmPrivate {
  /* Root communicator */
  MPI_Comm comm;
  int size;
  int rank;
  MashmBool isMasterProc;
  /* Intranodal communicator struct */
  MashmIntraNodeComm intraComm;
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

  int* sendAggMsgOffsets;
  int* recvAggMsgOffsets;

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

  /* Nodal message data for the MIN_AGG scheme */
  double* p_sendNodalSharedBuffer;
  double* p_recvNodalSharedBuffer;
  MPI_Win sendNodalSharedMemWindow;
  MPI_Win recvNodalSharedMemWindow;
  double** sendNodalSharedBufferIndex;
  double** recvNodalSharedBufferIndex;

};


/* Set up Internode Communication */
void p_mashmSetupInterNodeComm(struct MashmPrivate* p_mashm);

/* Allocate memory for the intranode messages */
void p_mashmAllocateSharedMemory(struct MashmPrivate* p_mashm, int bufferSize);

MashmBool p_MashmIsIntraNodeRank(struct MashmPrivate* p_mashm, int pairRank);
void p_MashmCalcMsgBufferSize(struct MashmPrivate* p_mashm);

/* Set up Intranode Communication */
void p_mashmSetupIntraMsgComm(struct MashmPrivate* p_mashm);
void p_mashmSetupIntraSharedComm(struct MashmPrivate* p_mashm);

/* Setup MIN_AGG data */
void p_mashmSetupAggType(struct MashmPrivate* p_mashm);
void p_mashmAllocateSharedMemoryMinAgg(struct MashmPrivate* p_mashm);
void p_mashmCalcMsgIndicesMinAgg(struct MashmPrivate* p_mashm);
void p_MashmCalcNumMpiMsgs(struct MashmPrivate* p_mashm);

/* Figure out nodal message scheduling for MASHM_COMM_MIN_AGG */
void p_mashmCalculateNodalMsgSchedule(struct MashmPrivate* p_mashm);

/* Return the buffer pointers for the specified message */
double* p_mashmGetBufferPointer(struct MashmPrivate* p_mashm, int msgIndex, MashmSendReceive sendReceive);

/* Retire (set to null) the buffer pointers */
void p_mashmRetireBufferPointer(struct MashmPrivate* p_mashm, double** bufPtr);

/* Communication Routines for MPI_Irecv/MPI_Isend - (all but MIN_AGG) */
void p_mashmStandardCommBegin(struct MashmPrivate* p_mashm);
void p_mashmStandardCommEnd(struct MashmPrivate* p_mashm);

/* Communication routines for intranode MPI_Irecv/MPI_Isend - INTRA_MSGS only */
void p_mashmIntraMsgsCommBegin(struct MashmPrivate* p_mashm);
void p_mashmIntraMsgsCommEnd(struct MashmPrivate* p_mashm);

/* Communication routines for intranode shared memory - INTRA_SHARED and MIN_AGG */
void p_mashmIntraSharedCommBegin(struct MashmPrivate* p_mashm);
void p_mashmIntraSharedCommEnd(struct MashmPrivate* p_mashm);

/* Communication routines for single nodal messages - MIN_AGG only */
void p_mashmMinAggCommBegin(struct MashmPrivate* p_mashm);
void p_mashmMinAggCommEnd(struct MashmPrivate* p_mashm);

/* Null function to satisfy function pointer assignment */
void p_nullFunction(struct MashmPrivate* p_mashm);

/* Destroy the MashmPrivate object */
void p_Destroy(struct MashmPrivate* p_mashm);

#endif
