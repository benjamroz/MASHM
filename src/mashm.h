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

/* Set up Intranode Communication */
void mashmSetupIntraNodeComm(Mashm in_mashm);

/* Set up Internode Communication */
void mashmSetupInterNodeComm(Mashm in_mashm);

/* Internode communication */
void mashmInterNodeCommBegin(Mashm in_mashm);
void mashmInterNodeCommEnd(Mashm in_mashm);

/* Intranode communication */
void mashmIntraNodeCommBegin(Mashm in_mashm);
void mashmIntraNodeCommEnd(Mashm in_mashm);

/* Private routines */
void p_mashmStandardCommBegin(_p_mashm* p_mashm);
void p_mashmStandardCommEnd(_p_mashm* p_mashm);

void p_mashmIntraMsgsCommBegin(_p_mashm* p_mashm);
void p_mashmIntraMsgsCommEnd(_p_mashm* p_mashm);

void p_mashmIntraSharedCommBegin(_p_mashm* p_mashm);
void p_mashmIntraSharedCommEnd(_p_mashm* p_mashm);

void p_mashmMinAggCommBegin(_p_mashm* p_mashm);
void p_mashmMinAggCommEnd(_p_mashm* p_mashm);

double* p_mashmGetBufferPointer(_p_mashm* p_mashm, int msgIndex, MashmSendReceive sendReceive);

void p_nullFunction(_p_mashm* p_mashm);

#endif 
