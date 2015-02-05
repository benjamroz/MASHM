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

  /* Across node sub-communicator (ie. intraComm.rank = 0 on all nodes) */
  MPI_Comm rankComm;

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

  int sumNumMsgs;
  int allMsgSizes;
  int intraMsgSizes;
  int interMsgSizes;

  int minIntraMsgSize;
  int maxIntraMsgSize;
  int minInterMsgSize;
  int maxInterMsgSize;
  int minNodalMsgSize;
  int maxNodalMsgSize;
  int sumNodalMsgSize;

};


/* Initialize the MashmPrivate object */
void p_MashmInit(struct MashmPrivate* p_mashm, MPI_Comm in_comm);

/* Set up Internode Communication */
void p_MashmSetupInterNodeComm(struct MashmPrivate* p_mashm);

/* Allocate memory for the intranode messages */
void p_MashmAllocateSharedMemory(struct MashmPrivate* p_mashm, int bufferSize);

MashmBool p_MashmIsIntraNodeRank(struct MashmPrivate* p_mashm, int pairRank);
void p_MashmCalcMsgBufferSize(struct MashmPrivate* p_mashm);

/* Get communicator */
MPI_Comm p_MashmGetComm(const struct MashmPrivate* p_mashm);

/* Get MPI communicator size */
int p_MashmGetSize(const struct MashmPrivate* p_mashm);

/* Get MPI communicator rank */
int p_MashmGetRank(const struct MashmPrivate* p_mashm);

/* Set up Intranode Communication */
void p_MashmSetupIntraMsgComm(struct MashmPrivate* p_mashm);
void p_MashmSetupIntraSharedComm(struct MashmPrivate* p_mashm);

/* Setup MIN_AGG data */
void p_MashmSetupAggType(struct MashmPrivate* p_mashm);
void p_MashmAllocateSharedMemoryMinAgg(struct MashmPrivate* p_mashm);
void p_MashmCalcMsgIndicesMinAgg(struct MashmPrivate* p_mashm);
void p_MashmCalcNumMpiMsgs(struct MashmPrivate* p_mashm);

/* Figure out nodal message scheduling for MASHM_COMM_MIN_AGG */
void p_MashmCalculateNodalMsgSchedule(struct MashmPrivate* p_mashm);

/* Return the buffer pointers for the specified message */
double* p_MashmGetBufferPointer(struct MashmPrivate* p_mashm, int msgIndex, MashmSendReceive sendReceive);

/* Return the buffer pointers for the specified rank and send/receive */
double* p_MashmGetBufferPointerForDest(const struct MashmPrivate* p_mashm, int destRank, MashmSendReceive sendReceive);

/* Retire (set to null) the buffer pointers */
void p_MashmRetireBufferPointer(struct MashmPrivate* p_mashm, double** bufPtr);

/* Communication Routines for MPI_Irecv/MPI_Isend - (all but MIN_AGG) */
void p_MashmStandardCommBegin(struct MashmPrivate* p_mashm);
void p_MashmStandardCommEnd(struct MashmPrivate* p_mashm);

/* Communication routines for intranode MPI_Irecv/MPI_Isend - INTRA_MSGS only */
void p_MashmIntraMsgsCommBegin(struct MashmPrivate* p_mashm);
void p_MashmIntraMsgsCommEnd(struct MashmPrivate* p_mashm);

/* Communication routines for intranode shared memory - INTRA_SHARED and MIN_AGG */
inline 
void p_MashmIntraSharedCommBegin(struct MashmPrivate* p_mashm);

inline
void p_MashmIntraSharedCommEnd(struct MashmPrivate* p_mashm);

/* Communication routines for single nodal messages - MIN_AGG only */
inline
void p_MashmMinAggCommBegin(struct MashmPrivate* p_mashm);

inline
void p_MashmMinAggCommEnd(struct MashmPrivate* p_mashm);

/* Null function to satisfy function pointer assignment */
void p_MashmNullFunction(struct MashmPrivate* p_mashm);

/* Complete the MashmPrivate struct - allocate data and set up pointers */
void p_MashmFinish(struct MashmPrivate* p_mashm);

/* Destroy the MashmPrivate object */
void p_MashmDestroy(struct MashmPrivate* p_mashm);

/* Print information about the MashmPrivate object */
void p_MashmPrintInfo(const struct MashmPrivate* p_mashm);

/* Set the number of communications */
void p_MashmSetNumComms(struct MashmPrivate* p_mashm, int numComms);

/* Set the communication data for the specified message index */
void p_MashmSetComm(struct MashmPrivate* p_mashm, int commIndex, int pairRank, int msgSize);

/* Get the rank of communication */
int p_MashmGetCommRank(struct MashmPrivate* p_mashm, int commIndex);

/* Get the size of communication */
int p_MashmGetCommSize(struct MashmPrivate* p_mashm, int commIndex);

/* Print the communication collection */
void p_MashmPrintCommCollection(const struct MashmPrivate* p_mashm);

/* Set the communication method */
//void p_MashmSetCommMethod(struct MashmPrivate* in_mashm, MashmCommType commType);
void p_MashmSetCommMethod(struct MashmPrivate* in_mashm, int commType);

/* Return the communication method */
MashmCommType p_MashmGetCommMethod(const struct MashmPrivate* p_mashm);

/* Whether the message is on node or not */
inline
MashmBool p_MashmIsMsgOnNode(const struct MashmPrivate* p_mashm, int msgIndex);

/* Begin internode communication */
inline
void p_MashmInterNodeCommBegin(struct MashmPrivate* p_mashm);

/* Begin intranode communication */
inline
void p_MashmIntraNodeCommBegin(struct MashmPrivate* p_mashm);

/* Finish intranode communication */
inline
void p_MashmIntraNodeCommEnd(struct MashmPrivate* p_mashm);

/* Finish internode communication */
inline
void p_MashmInterNodeCommEnd(struct MashmPrivate* p_mashm);

void p_MashmPrintInterNodeMessages(struct MashmPrivate* p_mashm);

void p_MashmPrintMessageInformation(struct MashmPrivate* p_mashm);
#endif
