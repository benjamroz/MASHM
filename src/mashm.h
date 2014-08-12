#ifndef MASHM_H
#define MASHM_H

#include "mpi.h"

typedef struct MashmPrivate _p_mashm;

typedef struct {
  _p_mashm* p;
} Mashm;

int mashmInit(Mashm* in_mashm, MPI_Comm in_comm);
MPI_Comm mashmGetComm(const Mashm in_mashm);
int mashmGetSize(const Mashm in_mashm);
int mashmGetRank(const Mashm in_mashm);

void mashmAddSymComm(Mashm in_mashm, int pairRank, int msgSize);
void mashmCommFinish(Mashm in_mashm);


void mashmInterNodeCommBegin(Mashm myMashm);
void mashmIntraNodeExchange(Mashm myMashm);
void mashmInterNodeCommEnd(Mashm myMashm);

#endif 
