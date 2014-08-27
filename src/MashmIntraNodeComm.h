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
} MashmIntraNodeComm;

int intraNodeInit(MashmIntraNodeComm* intraComm, MPI_Comm in_comm);

MPI_Comm intraNodeGetComm(const MashmIntraNodeComm intraComm);

int intraNodeGetSize(const MashmIntraNodeComm intraComm);

int intraNodeGetRank(const MashmIntraNodeComm intraComm);

int intraNodeDetermineGlobalInfo(MashmIntraNodeComm* intraComm);

int intraNodeGetSharedRank(const MashmIntraNodeComm intraComm, int pairRank);

int intraNodeDetermineNodalInfo(MashmIntraNodeComm* intraComm);

void intraNodePrintInfo(const MashmIntraNodeComm intraComm);

int intraNodeDestroy(MashmIntraNodeComm* intraComm);

#endif 
