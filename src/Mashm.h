#ifndef MASHM_H
#define MASHM_H

#include "mpi.h"

#include "MashmBool.h"
#include "MashmDefs.h"

#include "cmake_fortran_c_interface.h"

typedef struct {
  struct MashmPrivate* p;
} Mashm;

/* Initialize the Mashm type */
void MashmInit(Mashm* in_mashm, MPI_Comm in_comm);

/* Special binding for initializing using a Fortran MPI Communicator
 *   This gets mapped to mashmInit in Fortran */
#define MashmInitF2C FCI_GLOBAL(mashminit,MASHMINIT)
void MashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm); 

#define MashmGetComm FCI_GLOBAL(mashmgetcomm,MASHMGETCOMM)
MPI_Comm MashmGetComm(const Mashm in_mashm);

#define MashmGetSize FCI_GLOBAL(mashmgetsize,MASHMGETSIZE)
int MashmGetSize(const Mashm in_mashm);


#define MashmGetRank FCI_GLOBAL(MashmGetRank,MashmGetRank)
int MashmGetRank(const Mashm in_mashm);

#define MashmSetComm FCI_GLOBAL(mashmsetcomm,MASHMSETCOMM)
void MashmSetComm(Mashm in_mashm, int commIndex, int pairRank, int msgSize);

#define MashmCommFinish FCI_GLOBAL(mashmcommfinish,MASHMCOMMFINISH)
void MashmCommFinish(Mashm in_mashm);

#define MashmPrintInfo FCI_GLOBAL(mashmprintinfo,MASHMPRINTINFO)
void MashmPrintInfo(const Mashm in_mashm);

#define MashmIsMsgOnNode FCI_GLOBAL(mashmismsgonnode,MASHMISMSGONNODE)
MashmBool MashmIsMsgOnNode(Mashm in_mashm, int msgIndex);

#define MashmSetCommMethod FCI_GLOBAL(mashmsetcommmethod,MASHMSETCOMMMETHOD)
void MashmSetCommMethod(Mashm in_mashm, MashmCommType commType);

#define MashmSetNumComms FCI_GLOBAL(mashmsetnumcomms,MASHMSETNUMCOMMS)
void MashmSetNumComms(Mashm in_mashm, int numNeighbors);

#define MashmGetCommMethod FCI_GLOBAL(mashmgetcommmethod,MASHMGETCOMMMETHOD)
MashmCommType MashmGetCommMethod(Mashm in_mashm);

#define MashmNumMpiMsgs FCI_GLOBAL(mashmnummpimsgs,MASHMNUMMPIMSGS)
int MashmNumMpiMsgs(Mashm in_mashm);

#define MashmNumIntraNodeMsgs FCI_GLOBAL(mashmnumintranodemsgs,MASHMNUMINTRANODEMSGS)
int MashmNumIntraNodeMsgs(Mashm in_mashm);

#define MashmIsIntraNodeRank FCI_GLOBAL(mashmisintranoderank,MASHMISINTRANODERANK)
MashmBool MashmIsIntraNodeRank(Mashm in_mashm, int pairRank);

#define MashmCalcMsgBufferSize FCI_GLOBAL(mashmcalcmsgbuffersize,MASHMCALCMSGBUFFERSIZE)
void MashmCalcMsgBufferSize(Mashm in_mashm);

#define MashmCalcNumConnectedNodes FCI_GLOBAL(mashmcalcnumconnectednodes,MASHMCALCNUMCONNECTEDNODES)
void MashmCalcNumConnectedNodes(Mashm in_mashm);

#define MashmSetupStandardComm FCI_GLOBAL(mashmsetupstandardcomm,MASHMSETUPSTANDARDCOMM)
void MashmSetupStandardComm(Mashm in_mashm);

#define MashmGetBufferPointer FCI_GLOBAL(mashmgetbufferpointer,MASHMGETBUFFERPOINTER)
double* MashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive);

#define MashmGetBufferPointerForDest FCI_GLOBAL(mashmgetbufferpointerfordest,MASHMGETBUFFERPOINTERFORDEST)
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

#endif 
