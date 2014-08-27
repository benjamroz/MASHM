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

int MashmIntraNodeCommInit(MashmIntraNodeComm* intraComm, MPI_Comm in_comm);

MPI_Comm MashmIntraNodeCommGetComm(const MashmIntraNodeComm intraComm);

int MashmIntraNodeCommGetSize(const MashmIntraNodeComm intraComm);

int MashmIntraNodeCommGetRank(const MashmIntraNodeComm intraComm);

int MashmIntraNodeCommGetSharedRank(const MashmIntraNodeComm intraComm, int pairRank);

void MashmIntraNodeCommPrintInfo(const MashmIntraNodeComm intraComm);

int MashmIntraNodeCommDestroy(MashmIntraNodeComm* intraComm);

#endif 
