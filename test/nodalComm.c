#include "stdio.h"

#include <mpi.h>

#include "mashm.h"

int main(int argc, char** argv) {
  int ierr;
  Mashm myMashm;

  ierr = MPI_Init(&argc, &argv);

  /* Initialize the MASHM object */
  mashmInit(&myMashm, MPI_COMM_WORLD);

  mashmPrintInfo(myMashm);

  /* Add a communication */
  mashmAddSymComm(&myMashm, dest, size );

  /* Perform precalculation */
  mashmCommFinish(&myMashm);

  /* Allocate shared memory data */
  //shmMem = mashmAllocateSharedData(&myMashm);

  /* Fill shared memory */


  /* Send/Receive cycles */
  //mashmInterNodeCommBegin(myMashm);
  //mashmIntraNodeExchange(myMashm);
  //mashmInterNodeCommEnd(myMashm);


  //printf("Rank %d of size %d\n", mashmGetRank(myMashm), mashmGetSize(myMashm));


}
