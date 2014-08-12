#ifndef INTRA_NODE_COMM_H
#define INTRA_NODE_COMM_H

#include <mpi.h>

#include "mashmBool.h"

/* */

typedef struct {
  MPI_Comm comm;
  MPI_Comm parentComm;
  int size;
  int rank;
  mashmBool isMasterProc;
} intraNodeComm;

int init(intraNodeComm* intraComm, MPI_Comm in_comm);

MPI_Comm getComm(const intraNodeComm intraComm);

int getSize(const intraNodeComm intraComm);

int getRank(const intraNodeComm intraComm);

int determineGlobalInfo(intraNodeComm* intraComm);

int determineNodalInfo(intraNodeComm* intraComm);

void intraNodeCommPrintInfo(const intraNodeComm intraComm);

#endif 
