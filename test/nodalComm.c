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
  int testFailedInt, numTestsFailed;
  int* msgOffsets;
  int offset;

  ierr = MPI_Init(&argc,&argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  m = 10;
  n = 10;
  setbuf(stdout, NULL);
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

  msgOffsets = (int*) malloc(sizeof(int)*numNeighbors);
  msgOffsets[0] = 0;
  for (i = 1; i < numNeighbors; i++) {
    msgOffsets[i] = msgOffsets[i-1] + msgSizes[i-1];
  }

  /* Print element to process map */
  /*
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
  */

  /* Print communication graph */
  /*
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
  */

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
  origBuffer = (double*) malloc(sizeof(double)*sumMsgSizes);
  for (i = 0; i < numNeighbors; i++) {
    offset = msgOffsets[i];
    for (j = 0; j < msgSizes[i]; j++) {
      origBuffer[offset + j] = recvBuffers[i][j];
    }
  }
  /****************************************************************
   * Okay that was a lot of setup
   * Now use the MASHM library
   ****************************************************************/

  /* Initialize the MASHM object */
  if (rank == 0) printf("Initializing Mashm object\n");
  MashmInit(&myMashm, MPI_COMM_WORLD);

  /* Print nodal comm info */
  MashmPrintInfo(myMashm);

  /* Add communications calculated above */
  if (rank == 0) printf("Setting message information\n");
  MashmSetNumComms(myMashm, numNeighbors);

  for (i = 0; i < numNeighbors; i++) {
    MashmSetComm(myMashm, i, neighbors[i], msgSizes[i]);
  }

  if (argc > 1) {
    commMethodInt = atoi(argv[1]);
    if (commMethodInt > 4 || commMethodInt < 0) {
      printf("Error: commMethodInt should be >= 0 and <= 3. User supplied %s\n", argv[1]);
    }
  }
  else {
    commMethodInt = 0;
  }

  switch (commMethodInt) {
    case (0):
      if (rank == 0) printf("Mashm using commType = MASHM_COMM_STANDARD\n");
      commType = MASHM_COMM_STANDARD;
      break;
    case (1): 
      if (rank == 0) printf("Mashm using commType = MASHM_COMM_INTRA_MSG\n");
      commType = MASHM_COMM_INTRA_MSG;
      break;
    case (2): 
      if (rank == 0) printf("Mashm using commType = MASHM_COMM_INTRA_SHARED\n");
      commType = MASHM_COMM_INTRA_SHARED;
      break;
    case (3): 
      if (rank == 0) printf("Mashm using commType = MASHM_COMM_MIN_AGG\n");
      commType = MASHM_COMM_MIN_AGG;
      break;
  }

  MashmSetCommMethod(myMashm, commType);

  /* Perform precalculation */
  if (rank == 0) printf("Calling MashmCommFinish.\n");
  MashmCommFinish(myMashm);

  /* Print the comm collection */
  /*
  MashmPrintCommCollection(myMashm);
  */

  /* Retrieve pointers for buffers */
  mashmSendBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);
  mashmRecvBufferPtrs = (double**) malloc(sizeof(double*)*numNeighbors);


  if (rank == 0) printf("Getting buffer pointers.\n");
  /* Initialize the data to be sent */
  for (i = 0; i < numNeighbors; i++) {
    mashmSendBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_SEND);
    mashmRecvBufferPtrs[i] = MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE);
  }

 
  /*************************************************************
   * Now perform communication 
   ************************************************************/

  /* Fill internode buffers */
  for (i = 0; i < numNeighbors; i++) {
    if (! MashmIsMsgIntraNodal(myMashm, i)) {
      for (j = 0; j < msgSizes[i]; j++) {
        mashmSendBufferPtrs[i][j] = rank*msgSizes[i]+j;
      }
    }
  }

  if (rank == 0) printf("Mashm Internode communication begin.\n");
  /* Send internode messages */
  MashmInterNodeCommBegin(myMashm);

  /* Messages sent and receives posted 
   * Can asynchronously do work on nodal data 
   */
  for (i = 0; i < numNeighbors; i++) {
    if (MashmIsMsgIntraNodal(myMashm, i)) {
      for (j = 0; j < msgSizes[i]; j++) {
        mashmSendBufferPtrs[i][j] = rank*msgSizes[i]+j;
      }
    }
  }

  /* Send intranode messages */
  if (rank == 0) printf("Mashm Intranode communication begin.\n");
  MashmIntraNodeCommBegin(myMashm);

  /* Asynchronously do some computation */
  mashmData = (double*) malloc(sizeof(double)*sumMsgSizes);

  if (rank == 0) printf("Mashm Intranode communication end.\n");
  MashmIntraNodeCommEnd(myMashm);
  /* At this stage you have completed the intra-node communication */

  /* Asynchronously do work on nodal data */
  for (i = 0; i < numNeighbors; i++) {
    if (MashmIsMsgIntraNodal(myMashm, i)) {
      /* Unpack individual buffer */
      offset = msgOffsets[i];
      for (j = 0; j < msgSizes[i]; j++) {
        mashmData[offset+j] = mashmRecvBufferPtrs[i][j];
      }
    }
  }

  /* Now wait on nodal messages */
  if (rank == 0) printf("Mashm Internode communication end.\n");
  MashmInterNodeCommEnd(myMashm);

  for (i = 0; i < numNeighbors; i++) {
    if (! MashmIsMsgIntraNodal(myMashm, i)) {
      /* Unpack individual buffer */
      offset = msgOffsets[i];
      for (j = 0; j < msgSizes[i]; j++) {
        mashmData[offset+j] = mashmRecvBufferPtrs[i][j];
      }
    }
  }


  /****************************************************************
   * All communication has been performed
   * 
   * Now compare unpacked buffers
   ****************************************************************/

  if (rank == 0) printf("Testing Mashm data against original.\n");
  /* Test that the original buffer and the Mashm data are the same */
  for (iRank = 0; iRank < numProcs; iRank++) {
    if (iRank == rank) {
      /* printf("Process %d (orig, mashm, diff): \n", rank); */
      counter = 0;
      for (i = 0; i < numNeighbors; i++) {
        /* printf("  Message from %d to process %d:\n", rank, neighbors[i]); */
        for (j = 0; j < msgSizes[i]; j++) {
          /* printf(  "%f, %f, %f\n", origBuffer[counter], mashmData[counter], origBuffer[counter] - mashmData[counter]); */
          if (origBuffer[counter] - mashmData[counter] != 0.0) {
            /* printf("  Difference is significant! Test failed.\n"); */
            testFailed = true;
          }
          counter += 1;
        }
      }
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  if (testFailed) {
    testFailedInt = 1;
  }
  else {
    testFailedInt = 0;
  }

  ierr = MPI_Reduce(&testFailedInt, &numTestsFailed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) printf("Retiring Mashm buffer pointers.\n");
  /* Retire the Mashm buffer pointers */
  for (i = 0; i < numNeighbors; i++) {
    MashmRetireBufferPointer(myMashm, &(mashmSendBufferPtrs[i]));
    MashmRetireBufferPointer(myMashm, &(mashmRecvBufferPtrs[i]));
  }

  /* Destroy the Mashm object */
  if (rank == 0) printf("Calling Mashm Destroy.\n");
  MashmDestroy(&myMashm);

  decomp2dDestroyGraph(&neighbors, &msgSizes);


  free(mashmSendBufferPtrs);
  free(mashmRecvBufferPtrs);
 
  free(msgOffsets);

  ierr = MPI_Finalize();

  if (rank != 0) {
    return 0;
  }
  else if (numTestsFailed != 0) {
    printf("Test Failed: buffers are different.\n");
    return -1;
  } 
  else {
    printf("Test Passed: buffers are identical.\n");
    return 0;
  }
}
