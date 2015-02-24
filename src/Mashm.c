#include "stddef.h"

#include "mpi.h"

#include "MashmPrivate.h"
#include "Mashm.h"


/* Method to convert a Fortran MPI_Comm (integer) to  
 *   a C MPI_Comm type */
void MashmInitF2C(Mashm* in_mashm, MPI_Fint f_comm) {
  MPI_Comm c_comm;
  c_comm = MPI_Comm_f2c(f_comm);
  MashmInit(in_mashm, c_comm);
}

void MashmInit(Mashm* in_mashm, MPI_Comm in_comm) {

  /* Allocate the MashmPrivate structure */
  in_mashm->p = (struct MashmPrivate*) malloc(sizeof(struct MashmPrivate));

  /* Initialize the MashmPrivate structure */
  p_MashmInit(in_mashm->p, in_comm);

}

MPI_Comm MashmGetComm(const Mashm in_mashm) {
  return p_MashmGetComm(in_mashm.p);
}

int MashmGetSize(const Mashm in_mashm) {
  return p_MashmGetSize(in_mashm.p);
}

int MashmGetRank(const Mashm in_mashm) {
  return p_MashmGetRank(in_mashm.p);
}

void MashmPrintInfo(const Mashm in_mashm) {
  p_MashmPrintInfo(in_mashm.p);
}

void MashmSetNumComms(Mashm in_mashm, int numComms) {
  p_MashmSetNumComms(in_mashm.p, numComms);
}

void MashmSetComm(Mashm in_mashm, int commIndex, int pairRank, int msgSize) {
  p_MashmSetComm(in_mashm.p, commIndex, pairRank, msgSize);
}

int MashmGetCommRank(Mashm in_mashm, int commIndex) {
  return p_MashmGetCommRank(in_mashm.p, commIndex);
}

int MashmGetCommSize(Mashm in_mashm, int commIndex) {
  return p_MashmGetCommSize(in_mashm.p, commIndex);
}

void MashmPrintCommCollection(const Mashm in_mashm) {
  p_MashmPrintCommCollection(in_mashm.p);
}

void MashmSetCommMethod(Mashm in_mashm, MashmCommType commType) {
  p_MashmSetCommMethod(in_mashm.p, commType);
}

MashmCommType MashmGetCommMethod(const Mashm in_mashm) {
  return p_MashmGetCommMethod(in_mashm.p);
}

void MashmCommFinish(Mashm in_mashm) {
  p_MashmFinish(in_mashm.p);
}

inline
MashmBool MashmIsMsgIntraNodal(Mashm in_mashm, int msgIndex) {
  return p_MashmIsMsgOnNode(in_mashm.p, msgIndex);
}

double* MashmGetBufferPointer(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive) {
  /* Call the private routine */
  /* Debug
  printf("pointer, msgIndex, sendReceive %p,%d,%d\n", p_MashmGetBufferPointer(in_mashm.p, msgIndex, sendReceive),msgIndex, sendReceive);
  */
  return p_MashmGetBufferPointer(in_mashm.p, msgIndex, sendReceive);
}

void MashmGetBufferPointer2(Mashm in_mashm, int msgIndex, MashmSendReceive sendReceive, double** outPtr) {
  /* Call the private routine */
  /*
  printf("pointer, msgIndex, sendReceive %p,%d,%d\n", p_MashmGetBufferPointer(in_mashm.p, msgIndex, sendReceive),msgIndex, sendReceive);
  */
  *outPtr = p_MashmGetBufferPointer(in_mashm.p, msgIndex, sendReceive);
}

void MashmRetireBufferPointer(Mashm in_mashm, double** bufPtr) {
  /* Call the private routine */
  return p_MashmRetireBufferPointer(in_mashm.p, bufPtr);
}

double* MashmGetBufferPointerForDest(Mashm in_mashm, int destRank, MashmSendReceive sendReceive) {
  return p_MashmGetBufferPointerForDest(in_mashm.p, destRank, sendReceive);
}

void MashmDestroy(Mashm* in_mashm) {

  /* Destroy the MashmPrivate data */
  p_Destroy(in_mashm->p);

  /* Free the pointer to the MashmPrivate data */
  free(in_mashm->p);

}

inline
void MashmInterNodeCommBegin(Mashm in_mashm) {
  p_MashmInterNodeCommBegin(in_mashm.p);
}

inline
void MashmIntraNodeCommBegin(Mashm in_mashm) {
  p_MashmIntraNodeCommBegin(in_mashm.p);
}

inline
void MashmIntraNodeCommEnd(Mashm in_mashm) {
  p_MashmIntraNodeCommEnd(in_mashm.p);
}

inline
void MashmInterNodeCommEnd(Mashm in_mashm) {
  p_MashmInterNodeCommEnd(in_mashm.p);
}

inline
MashmBool MashmIsIntraNodeRank(const Mashm in_mashm, int pairRank) {
  return p_MashmIsIntraNodeRank(in_mashm.p, pairRank);
}

void MashmPrintMessageStats(const Mashm in_mashm) {
  p_MashmPrintMessageInformation(in_mashm.p);
}
