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
  int sumMsgSizes;
  int numNeighbors;

  /* MPI Communication boilerplate */
  MPI_Request* recvRequests;
  MPI_Request* sendRequests;
  MPI_Status* recvStatuses;
  MPI_Status* sendStatuses;

  /* Send, Recv buffers */
  double** recvBuffers;
  double** sendBuffers;

  MPI_Comm commWorld;
  Mashm myMashm;

  /* */
  double** mashmSendBufferPtrs;
  double** mashmRecvBufferPtrs;

  int commMethodInt;
  MashmCommType commType;

  int counter;
  double* origBuffer;
  double* mashmData;
  MashmBool testFailed = false;

  ierr = MPI_Init(&argc,&argv);
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
  /* Determine the total size of all messages to be sent */
  sumMsgSizes = 0;
  for (i = 0; i < numNeighbors; i++) {
    sumMsgSizes += msgSizes[i]; 
  }

  /* Allocate send and receive buffers */
  recvBuffers = (double**) malloc(sizeof(double*)*numNeighbors);
  sendBuffers = (double**) malloc(sizeof(double*)*numNeighbors);

  recvRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNeighbors);
  sendRequests = (MPI_Request*) malloc(sizeof(MPI_Request)*numNeighbors);
  recvStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNeighbors);
  sendStatuses = (MPI_Status*) malloc(sizeof(MPI_Status)*numNeighbors);

  /* Symmetric */
  for (i = 0; i < numNeighbors; i++) {
    recvBuffers[i] = (double*) malloc(sizeof(double)*msgSizes[i]);
    sendBuffers[i] = (double*) malloc(sizeof(double)*msgSizes[i]);
  }

  /* Now fill the buffers */
  for (i = 0; i < numNeighbors; i++) {
    for (j = 0; j < msgSizes[i]; j++) {
      sendBuffers[i][j] = rank*msgSizes[i]+j;
    }
  }

  /* Usual point to point communication */
  for (i = 0; i < numNeighbors; i++) {
    ierr = MPI_Irecv(recvBuffers[i], msgSizes[i], MPI_DOUBLE, neighbors[i], 0, MPI_COMM_WORLD, &recvRequests[i]);
  }
  for (i = 0; i < numNeighbors; i++) {
    ierr = MPI_Isend(sendBuffers[i], msgSizes[i], MPI_DOUBLE, neighbors[i], 0, MPI_COMM_WORLD, &sendRequests[i]);
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
  MashmSetNumComms(myMashm, numNeighbors);
  for (i = 0; i < numNeighbors; i++) {
    MashmSetComm(myMashm, i, neighbors[i], msgSizes[i]);
  }

  if (argc > 1) {
    commMethodInt = atoi(argv[1]);
    if (commMethodInt > 4 || commMethodInt < 0) {
      printf("Error: commMethodInt >= 0 <= 3. User supplied %s\n", argv[1]);
    }
  }
  else {
    commMethodInt = 0;
  }

  switch (commMethodInt) {
    case (0):
      printf("Using commType = MASHM_COMM_STANDARD\n");
      commType = MASHM_COMM_STANDARD;
      break;
    case (1): 
      printf("Using commType = MASHM_COMM_INTRA_MSG\n");
      commType = MASHM_COMM_INTRA_MSG;
      break;
    case (2): 
      printf("Using commType = MASHM_COMM_INTRA_SHARED\n");
      commType = MASHM_COMM_INTRA_SHARED;
      break;
    case (3): 
      printf("Using commType = MASHM_COMM_MIN_AGG\n");
      commType = MASHM_COMM_MIN_AGG;
      break;
  }

  MashmSetCommMethod(myMashm, commType);

  /* Perform precalculation */
  MashmCommFinish(myMashm);
  /* Retrieve pointers for buffers */
  mashmSendBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);
  mashmRecvBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);


  /* Initialize the data to be sent */
  for (i = 0; i < numNeighbors; i++) {
    mashmSendBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_SEND);
    mashmRecvBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE);
    for (j = 0; j < msgSizes[i]; j++) {
      mashmSendBufferPtrs[i][j] = rank*msgSizes[i]+j;
    }
  }

 
  /*************************************************************
   * Now perform communication 
   ************************************************************/

  for (i = 0; i < numNeighbors; i++) {
    /* Fill individual buffer */
    for (j = 0; j < msgSizes[i]; j++) {
    }
  }


  /* Fill internode buffers */
  for (i = 0; i < numNeighbors; i++) {
    if (! MashmIsMsgOnNode(myMashm, i)) {
    }
  }

  /* Send internode messages */
  MashmInterNodeCommBegin(myMashm);

  /* Messages sent and receives posted 
   * Can asynchronously do work on nodal data 
   */
  for (i = 0; i < numNeighbors; i++) {
    if (MashmIsMsgOnNode(myMashm, i)) {
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


  /****************************************************************
   * All communication has been performed
   * 
   * Now compare unpacked buffers
   */
  origBuffer = (double*) malloc(sizeof(double)*sumMsgSizes);
  for (i = 0; i < sumMsgSizes; i++) {
    origBuffer[i] = 666.0;
  }

  counter = 0;
  for (i = 0; i < numNeighbors; i++) {
    for (j = 0; j < msgSizes[i]; j++) {
      origBuffer[counter] = recvBuffers[i][j];
      counter += 1;
    }
  }

  mashmData = (double*) malloc(sizeof(double)*sumMsgSizes);
  for (i = 0; i < sumMsgSizes; i++) {
    mashmData[i] = 666.0;
  }


  counter = 0;
  for (i = 0; i < numNeighbors; i++) {
    for (j = 0; j < msgSizes[i]; j++) {
      mashmData[counter] = mashmRecvBufferPtrs[i][j];
      counter += 1;
    }
  }

  /* Test that the original buffer and the Mashm data are the same */
  for (iRank = 0; iRank < numProcs; iRank++) {
    if (iRank == rank) {
      printf("Process %d (orig, mashm, diff): \n", rank);
      for (i = 0; i < sumMsgSizes; i++) {
        printf(  "%f, %f, %f\n", origBuffer[i], mashmData[i], origBuffer[i] - mashmData[i]);
        /* Non-equivalence to zero is okay here */
        if (origBuffer[i] - mashmData[i] != 0.0) {
          printf("  Difference is significant! Test failed.\n");
          testFailed = true;
        }
      }

    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  /* Destroy the Mashm object */
  MashmDestroy(&myMashm);

  decomp2dDestroyGraph(&neighbors, &msgSizes);

  /* Hacky failure for now */
  if (testFailed) {
    return -1;
  }
  
  ierr = MPI_Finalize();

  return 0;
}
