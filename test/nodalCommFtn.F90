

program nodalCommFtn
#include "mashmf.h"
use mashm_mod
use mpi
implicit none
Mashm :: myMashm
integer :: ierr
integer :: myComm
integer :: i, numNeighbors
real*8, allocatable :: mashmSendBufferPtrs(:)
real*8, allocatable :: mashmRecvBufferPtrs(:)
integer :: neighbors, msgSizes

call MPI_Init(ierr)

myComm = MPI_COMM_WORLD
call mashmInit(myMashm, MPI_COMM_WORLD)

!/* Print nodal comm info */
call mashmPrintInfo(myMashm)

!/* Add communications calculated above */
do i = 1, numNeighbors
  call mashmAddSymComm(myMashm, neighbors(i), msgSizes(i))
enddo

!/* Perform precalculation */
call mashmCommFinish(myMashm)
!/* Retrieve pointers for buffers */
allocate(mashmSendBufferPtrs(numNeighbors))
allocate(mashmRecvBufferPtrs(numNeighbors))

!do i = 1, numNeighbors
!  mashmSendBufferPtrs(i) = mashmGetBufferPointer(myMashm, i, MASHM_SEND)
!  mashmRecvBufferPtrs(i) = mashmGetBufferPointer(myMashm, i, MASHM_RECEIVE)
!enddo

!/* Fill buffers */

!/*************************************************************
! * Now perform communication 
! ************************************************************/

!/* Send internode messages */
call mashmInterNodeCommBegin(myMashm)

!/* Messages sent and receives posted 
! * Can asynchronously do work on nodal data 
! */
!/* Send intranode messages */
call mashmIntraNodeCommBegin(myMashm)
call mashmIntraNodeCommEnd(myMashm)
!/* At this stage you have completed the intra-node communication */

!/* Asynchronously do work on nodal data */

!/* Now wait on nodal messages */
call mashmInterNodeCommEnd(myMashm)

!/* Destroy the Mashm object */
call mashmDestroy(myMashm)

!decomp2dDestroyGraph(neighbors, &msgSizes)

call MPI_Finalize(ierr)

end program nodalCommFtn
