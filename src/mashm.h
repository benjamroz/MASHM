#ifndef MASHM_H
#define MASHM_H

#include "mpi.h"

#include "mashmBool.h"

#include "cmake_fortran_c_interface.h"

typedef enum { MASHM_COMM_STANDARD, MASHM_COMM_INTRA_MSG, MASHM_COMM_INTRA_SHARED, MASHM_COMM_MIN_AGG } MashmCommType;

typedef enum { MASHM_SEND, MASHM_RECEIVE } MashmSendReceive;

struct MashmPrivate;
typedef struct MashmPrivate _p_mashm;

typedef struct {
  _p_mashm* p;
} Mashm;


/* Initialize the Mashm type */
void mashmInit(Mashm* in_mashm, MPI_Comm in_comm);

/* Special binding for initializing using a Fortran MPI Communicator
 *   This gets mapped to mashmInit in Fortran */
#define mashmInitF2C FCI_GLOBAL(mashminit,MASHMINIT)
void mashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm); 

MPI_Comm mashmGetComm(const Mashm in_mashm);
int mashmGetSize(const Mashm in_mashm);
int mashmGetRank(const Mashm in_mashm);

#define mashmAddSymComm FCI_GLOBAL(mashmaddsymcomm,MASHMADDSYMCOMM)
void mashmAddSymComm(Mashm in_mashm, int pairRank, int msgSize);


#define mashmCommFinish FCI_GLOBAL(mashmcommfinish,MASHMCOMMFINISH)
void mashmCommFinish(Mashm in_mashm);

#define mashmPrintInfo FCI_GLOBAL(mashmprintinfo,MASHMPRINTINFO)
void mashmPrintInfo(const Mashm in_mashm);

#define mashmSetCommMethod FCI_GLOBAL(mashmsetcommmethod,MASHMSETCOMMMETHOD)
void mashmSetCommMethod(Mashm in_mashm, MashmCommType commType);


MashmCommType mashmGetCommMethod(Mashm in_mashm);


int mashmNumMpiMsgs(Mashm in_mashm);
int mashmNumIntraNodeMsgs(Mashm in_mashm);

void mashmCalcNumMpiMsgs(Mashm in_mashm);
MashmBool mashmIsIntraNodeRank(Mashm in_mashm, int pairRank);
void mashmCalcMsgBufferSize(Mashm in_mashm);
void mashmCalcNumConnectedNodes(Mashm in_mashm);
void mashmSetupStandardComm(Mashm in_mashm);


#define mashmGetBufferPointer FCI_GLOBAL(mashmgetbufferpointer,MASHMGETBUFFERPOINTER)
double* mashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive);

double* mashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive);

/* Set up Intranode Communication */
void mashmSetupIntraNodeComm(Mashm in_mashm);

/* Set up Internode Communication */
void mashmSetupInterNodeComm(Mashm in_mashm);

/* Internode communication */
#define mashmInterNodeCommBegin FCI_GLOBAL(mashminternodecommbegin,MASHMINTERNODECOMMBEGIN)
void mashmInterNodeCommBegin(Mashm in_mashm);

#define mashmInterNodeCommEnd FCI_GLOBAL(mashminternodecommend,MASHMINTERNODECOMMEND)
void mashmInterNodeCommEnd(Mashm in_mashm);

/* Intranode communication */
#define mashmIntraNodeCommBegin FCI_GLOBAL(mashmintranodecommbegin,MASHMINTRANODECOMMBEGIN)
void mashmIntraNodeCommBegin(Mashm in_mashm);
#define mashmIntraNodeCommEnd FCI_GLOBAL(mashmintranodecommend,MASHMINTRANODECOMMEND)
void mashmIntraNodeCommEnd(Mashm in_mashm);

#define mashmDestroy FCI_GLOBAL(mashmdestroy,MASHMDESTROY)
void mashmDestroy(Mashm* in_mashm);

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
