#ifndef MASHM_H
#define MASHM_H

#include "mpi.h"

#include "mashmBool.h"

typedef enum { MASHM_COMM_STANDARD, MASHM_COMM_INTRA_MSG, MASHM_COMM_INTRA_SHARED, MASHM_COMM_MIN_AGG } MashmCommType;

typedef enum { MASHM_SEND, MASHM_RECEIVE } MashmSendReceive;

struct MashmPrivate;
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

void mashmSetCommMethod(Mashm in_mashm, MashmCommType commType);
MashmCommType mashmGetCommMethod(Mashm in_mashm);


int mashmNumMpiMsgs(Mashm in_mashm);
int mashmNumIntraNodeMsgs(Mashm in_mashm);

void mashmCalcNumMpiMsgs(Mashm in_mashm);
MashmBool mashmIsIntraNodeRank(Mashm in_mashm, int pairRank);
void mashmCalcMsgBufferSize(Mashm in_mashm);
void mashmCalcNumConnectedNodes(Mashm in_mashm);
void mashmSetupStandardComm(Mashm in_mashm);

double* mashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive);

double* mashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive);

void mashmSetupIntraNodeComm(Mashm in_mashm);
void mashmSetupInterNodeComm(Mashm in_mashm);

void mashmStandardCommBegin(_p_mashm* p_mashm);
void mashmStandardCommEnd(_p_mashm* p_mashm);

void mashmIntraMsgsCommBegin(Mashm in_mashm);
void mashmIntraMsgsCommEnd(Mashm in_mashm);

void mashmIntraNodeCommBegin(Mashm in_mashm);
void mashmIntraNodeCommEnd(Mashm in_mashm);

double* p_mashmGetBufferPointer(_p_mashm* p_mashm, int msgIndex, MashmSendReceive sendReceive);
void p_nullFunction(_p_mashm* p_mashm);
#endif 
