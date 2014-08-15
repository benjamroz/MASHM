#ifndef INTRA_NODE_COMM_H
#define INTRA_NODE_COMM_H

#include "mpi.h"

#include "mashmBool.h"

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

int intraNodeDetermineNodalInfo(intraNodeComm* intraComm);

void intraNodePrintInfo(const intraNodeComm intraComm);

int intraNodeDestroy(intraNodeComm* intraComm);

#endif 
