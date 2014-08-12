#include "stdio.h"

#include <mpi.h>

#include "mashm.h"

#include "decomp2d.h"

int main(int argc, char** argv) {
  int ierr;
  int rank, numProcs;
  int numElems;
  int m,n;
  int mIndex, nIndex;
  int i,j;
  int owner;
  int iRank;
  int calcRank;

  /* Communication graph */
  int* elements;
  int* neighbors;
  int* msgSizes;
  int numNeighbors;

  Mashm myMashm;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  /* Initialize the MASHM object */
  mashmInit(&myMashm, MPI_COMM_WORLD);

  m = 10;
  n = 9;
  numElems = decomp2dRectNumElements(m, n, rank, numProcs);

  printf("rank %d has %d elements\n", rank, numElems);

  if (rank == 0) {
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        owner = decomp2dGetOwnerRank(i,j,numProcs,m,n);
        printf("Element %d,%d owned by proc %d\n", i,j, owner);
      }
    }

  }

  decomp2dCreateGraph(m, n, rank, numProcs, &elements, &neighbors, &msgSizes, &numNeighbors);

  for (iRank = 0; iRank < numProcs; iRank++) {
    if (iRank == rank) {
      printf("Process %d owns the folloing elements:\n", rank);
      for (i = 0; i < numElems; i++) {
        printf("  elem %d\n", elements[i]);
        // Calc
        //mIndex = elements[i] / n;
        //nIndex = elements[i] % n;
        //calcRank = decomp2dGetOwnerRank(mIndex, nIndex, numProcs, m, n);
        //printf("  elem %d, calc %d, m %d, n %d\n", elements[i], calcRank, mIndex, nIndex);
      }
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  ierr = MPI_Barrier(MPI_COMM_WORLD);


  for (iRank = 0; iRank < numProcs; iRank++) {
    if (iRank == rank) {
      printf("Process %d has %d mpi neighbors\n", rank, numNeighbors);
      for (i = 0; i < numNeighbors; i++) {
        printf("Process %d communicates with process %d size %d\n", rank, neighbors[i], msgSizes[i]);
      }
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  mashmPrintInfo(myMashm);

  /* Add a communication */
  //mashmAddSymComm(myMashm, dest, size );

  /* Perform precalculation */
  mashmCommFinish(myMashm);

  /* Allocate shared memory data */
  //shmMem = mashmAllocateSharedData(&myMashm);

  /* Fill shared memory */


  /* Send/Receive cycles */
  //mashmInterNodeCommBegin(myMashm);
  //mashmIntraNodeExchange(myMashm);
  //mashmInterNodeCommEnd(myMashm);


  //printf("Rank %d of size %d\n", mashmGetRank(myMashm), mashmGetSize(myMashm));

  decomp2dDestroyGraph(&neighbors, &msgSizes);

}
