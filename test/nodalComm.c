#include "stdio.h"

#include <mpi.h>

#include "mashm.h"

int main(int argc, char** argv) {
  int ierr;
  Mashm myMashm;

  ierr = MPI_Init(&argc, &argv);

  mashmInit(&myMashm, MPI_COMM_WORLD);

  mashmPrintInfo(myMashm);

  //printf("Rank %d of size %d\n", mashmGetRank(myMashm), mashmGetSize(myMashm));


}
