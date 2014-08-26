#ifndef INTRA_NODE_COMM_H
#define INTRA_NODE_COMM_H

#include "stdio.h"
#include "stdlib.h"

#include "mpi.h"

#include "MashmBool.h"

typedef struct {
  MPI_Comm comm;
  MPI_Comm parentComm;
  int size;
  int rank;
  MashmBool isMasterProc;
  int parentRank;
  int* parentRanksOnNode;
} intraNodeComm;

int intraNodeInit(intraNodeComm* intraComm, MPI_Comm in_comm);

MPI_Comm intraNodeGetComm(const intraNodeComm intraComm);

int intraNodeGetSize(const intraNodeComm intraComm);

int intraNodeGetRank(const intraNodeComm intraComm);

int intraNodeDetermineGlobalInfo(intraNodeComm* intraComm);

int intraNodeGetSharedRank(const intraNodeComm intraComm, int pairRank);

int intraNodeDetermineNodalInfo(intraNodeComm* intraComm);

void intraNodePrintInfo(const intraNodeComm intraComm);

int intraNodeDestroy(intraNodeComm* intraComm);

#endif 
