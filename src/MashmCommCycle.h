#ifndef MASHM_COMM_CYCLE_H
#define MASHM_COMM_CYCLE_H

#include "stdio.h"
#include "stdlib.h"

#include "MashmBool.h"

/* Structure representing one point to point communication */
typedef struct {
  int pairRank;
  int sendSize;
  int recvSize;
} MashmComm;

void MashmCommInit(MashmComm* in_comm);

void MashmCommSetComm(MashmComm* in_comm, int pairRank, int sendSize, int recvSize);

int MashmCommGetPairRank(const MashmComm in_comm);

int MashmCommGetSendSize(const MashmComm in_comm);

int MashmCommGetRecvSize(const MashmComm in_comm);

void MashmCommCopy(const MashmComm in_comm, MashmComm* copy_comm);

void MashmCommPrint(const MashmComm in_comm);

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
void MashmCommCollectionInit(MashmCommCollection* commCollection);

/**
 * Set the size of the comm collection
 */
void MashmCommCollectionSetSize(MashmCommCollection* commCollection, int numComms);

void MashmCommCollectionSetComm(MashmCommCollection* commCollection, int commIndex, int pairRank, int sendSize, int recvSize);

void MashmCommCollectionPrint(const MashmCommCollection commCollection);

void MashmCommCollectionReset(MashmCommCollection* commCollection);

void MashmCommCollectionDestroy(MashmCommCollection* commCollection);


#endif 
