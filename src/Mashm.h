#ifndef MASHM_H
#define MASHM_H

#include "mpi.h"

#include "MashmBool.h"

#include "cmake_fortran_c_interface.h"

typedef enum { MASHM_COMM_STANDARD, MASHM_COMM_INTRA_MSG, MASHM_COMM_INTRA_SHARED, MASHM_COMM_MIN_AGG } MashmCommType;

typedef enum { MASHM_SEND, MASHM_RECEIVE } MashmSendReceive;

struct MashmPrivate;
typedef struct MashmPrivate _p_mashm;

typedef struct {
  _p_mashm* p;
} Mashm;


/* Initialize the Mashm type */
void MashmInit(Mashm* in_mashm, MPI_Comm in_comm);

/* Special binding for initializing using a Fortran MPI Communicator
 *   This gets mapped to mashmInit in Fortran */
#define MashmInitF2C FCI_GLOBAL(mashminit,MASHMINIT)
void MashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm); 

MPI_Comm MashmGetComm(const Mashm in_mashm);
int MashmGetSize(const Mashm in_mashm);
int MashmGetRank(const Mashm in_mashm);

#define MashmAddSymComm FCI_GLOBAL(mashmaddsymcomm,MASHMADDSYMCOMM)
void MashmAddSymComm(Mashm in_mashm, int pairRank, int msgSize);


#define MashmCommFinish FCI_GLOBAL(mashmcommfinish,MASHMCOMMFINISH)
void MashmCommFinish(Mashm in_mashm);

#define MashmPrintInfo FCI_GLOBAL(mashmprintinfo,MASHMPRINTINFO)
void MashmPrintInfo(const Mashm in_mashm);

#define MashmIsMsgOnNode FCI_GLOBAL(mashmismsgonnode,MASHMISMSGONNODE)
MashmBool MashmIsMsgOnNode(Mashm in_mashm, int msgIndex);

#define MashmSetCommMethod FCI_GLOBAL(mashmsetcommmethod,MASHMSETCOMMMETHOD)
void MashmSetCommMethod(Mashm in_mashm, MashmCommType commType);


MashmCommType MashmGetCommMethod(Mashm in_mashm);


int MashmNumMpiMsgs(Mashm in_mashm);
int MashmNumIntraNodeMsgs(Mashm in_mashm);

void MashmCalcNumMpiMsgs(Mashm in_mashm);
MashmBool MashmIsIntraNodeRank(Mashm in_mashm, int pairRank);
void MashmCalcMsgBufferSize(Mashm in_mashm);
void MashmCalcNumConnectedNodes(Mashm in_mashm);
void MashmSetupStandardComm(Mashm in_mashm);


#define MashmGetBufferPointer FCI_GLOBAL(mashmgetbufferpointer,MASHMGETBUFFERPOINTER)
double* MashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive);

double* MashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive);

/* Internode communication */
#define MashmInterNodeCommBegin FCI_GLOBAL(mashminternodecommbegin,MASHMINTERNODECOMMBEGIN)
void MashmInterNodeCommBegin(Mashm in_mashm);

#define MashmInterNodeCommEnd FCI_GLOBAL(mashminternodecommend,MASHMINTERNODECOMMEND)
void MashmInterNodeCommEnd(Mashm in_mashm);

/* Intranode communication */
#define MashmIntraNodeCommBegin FCI_GLOBAL(mashmintranodecommbegin,MASHMINTRANODECOMMBEGIN)
void MashmIntraNodeCommBegin(Mashm in_mashm);
#define MashmIntraNodeCommEnd FCI_GLOBAL(mashmintranodecommend,MASHMINTRANODECOMMEND)
void MashmIntraNodeCommEnd(Mashm in_mashm);

#define MashmDestroy FCI_GLOBAL(mashmdestroy,MASHMDESTROY)
void MashmDestroy(Mashm* in_mashm);

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

void p_mashmAllocateSharedMemory(struct MashmPrivate* p_mashm, int bufferSize);

/* Set up Intranode Communication */
void p_mashmSetupIntraMsgComm(struct MashmPrivate* p_mashm);
void p_mashmSetupIntraSharedComm(struct MashmPrivate* p_mashm);

/* Set up Internode Communication */
void p_mashmSetupInterNodeComm(_p_mashm* p_mashm);

#endif 
