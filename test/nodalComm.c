#include "stdio.h"

#include "mpi.h"

#include "Mashm.h"

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

  /* MPI Communication boilerplate */
  MPI_Request* recvRequests;
  MPI_Request* sendRequests;
  MPI_Status* recvStatuses;
  MPI_Status* sendStatuses;

  /* Send, Recv buffers */
  int** recvBuffers;
  int** sendBuffers;

  MPI_Comm commWorld;
  Mashm myMashm;

  /* */
  double** mashmSendBufferPtrs;
  double** mashmRecvBufferPtrs;

  ierr = MPI_Init(NULL, NULL);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  m = 10;
  n = 10;

#if 0
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
#endif

  decomp2dCreateGraph(m, n, rank, numProcs, &numElems, &elements, &neighbors, &msgSizes, &numNeighbors);

  /* Print element to process map */
  for (iRank = 0; iRank < numProcs; iRank++) {
    if (iRank == rank) {
      printf("Process %d owns the following elements:\n", rank);
      for (i = 0; i < numElems; i++) {
        printf("  elem %d\n", elements[i]);
      }
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  ierr = MPI_Barrier(MPI_COMM_WORLD);

  /* Print communication graph */
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

  /* Allocate send and receive buffers */
  recvBuffers = (int**) malloc(sizeof(int*)*numNeighbors);
  sendBuffers = (int**) malloc(sizeof(int*)*numNeighbors);

  recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNeighbors);
  sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNeighbors);
  recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNeighbors);
  sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNeighbors);

  /* Symmetric */
  for (i = 0; i < numNeighbors; i++) {
    recvBuffers[i] = (int*) malloc(sizeof(int)*msgSizes[i]);
    sendBuffers[i] = (int*) malloc(sizeof(int)*msgSizes[i]);
  }

  /* Now fill the buffers */

  /* Usual point to point communication */
  for (i = 0; i < numNeighbors; i++) {
    ierr = MPI_Irecv(recvBuffers[i], msgSizes[i], MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &recvRequests[i]);
  }
  for (i = 0; i < numNeighbors; i++) {
    ierr = MPI_Isend(recvBuffers[i], msgSizes[i], MPI_INT, neighbors[i], 0, MPI_COMM_WORLD, &sendRequests[i]);
  }

  ierr = MPI_Waitall(numNeighbors,recvRequests,recvStatuses); 
  ierr = MPI_Waitall(numNeighbors,sendRequests,sendStatuses); 

  /****************************************************************
   * Okay that was a lot of setup
   * Now use the MASHM library
   ****************************************************************/

  /* Initialize the MASHM object */

  commWorld = MPI_COMM_WORLD;
  MashmInit(&myMashm, commWorld);

  /* Print nodal comm info */
  MashmPrintInfo(myMashm);

  /* Add communications calculated above */
  for (i = 0; i < numNeighbors; i++) {
    MashmAddSymComm(myMashm, neighbors[i], msgSizes[i]);
  }

  //MashmSetCommMethod(myMashm, MASHM_COMM_STANDARD);
  //MashmSetCommMethod(myMashm, MASHM_COMM_INTRA_MSG);
  //MashmSetCommMethod(myMashm, MASHM_COMM_INTRA_SHARED);
  MashmSetCommMethod(myMashm, MASHM_COMM_MIN_AGG);

  /* Perform precalculation */
  MashmCommFinish(myMashm);
  return 0;
  /* Retrieve pointers for buffers */
  mashmSendBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);
  mashmRecvBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);
  for (i = 0; i < numNeighbors; i++) {
    mashmSendBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_SEND);
    mashmRecvBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE);
  }

 
  /*************************************************************
   * Now perform communication 
   ************************************************************/
  /* Fill internode buffers */
  for (i = 0; i < numNeighbors; i++) {
    if (! MashmIsMsgOnNode(myMashm, i)) {
      /* Fill individual buffer */
    }
  }

  /* Send internode messages */
  MashmInterNodeCommBegin(myMashm);

  /* Messages sent and receives posted 
   * Can asynchronously do work on nodal data 
   */
  for (i = 0; i < numNeighbors; i++) {
    if (MashmIsMsgOnNode(myMashm, i)) {
      /* Fill individual buffer */
    }
  }

  /* Send intranode messages */
  MashmIntraNodeCommBegin(myMashm);
  MashmIntraNodeCommEnd(myMashm);
  /* At this stage you have completed the intra-node communication */

  /* Asynchronously do work on nodal data */
  for (i = 0; i < numNeighbors; i++) {
    if (MashmIsMsgOnNode(myMashm, i)) {
      /* Unpack individual buffer */
    }
  }

  /* Now wait on nodal messages */
  MashmInterNodeCommEnd(myMashm);

  for (i = 0; i < numNeighbors; i++) {
    if (! MashmIsMsgOnNode(myMashm, i)) {
      /* Unpack individual buffer */
    }
  }
  /* Destroy the Mashm object */
  MashmDestroy(&myMashm);

  decomp2dDestroyGraph(&neighbors, &msgSizes);

  ierr = MPI_Finalize();

  return 0;
}
