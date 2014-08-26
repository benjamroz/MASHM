

program nodalCommFtn
#include "Mashmf.h"
use Mashm_mod
use mpi
implicit none
Mashm :: myMashm
integer :: ierr
integer :: myComm
integer :: i, numNeighbors
real*8, allocatable :: mashmSendBufferPtrs(:)
real*8, allocatable :: mashmRecvBufferPtrs(:)
integer, allocatable :: neighbors(:), msgSizes(:)

call MPI_Init(ierr)

myComm = MPI_COMM_WORLD
call MashmInit(myMashm, MPI_COMM_WORLD)

!/* Print nodal comm info */
call MashmPrintInfo(myMashm)

call MashmSetNumComms(myMashm, numNeighbors)

!/* Add communications calculated above */
do i = 1, numNeighbors
  call MashmSetComm(myMashm, i, neighbors(i), msgSizes(i))
enddo

!/* Perform precalculation */
call MashmCommFinish(myMashm)
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
call MashmInterNodeCommBegin(myMashm)

!/* Messages sent and receives posted 
! * Can asynchronously do work on nodal data 
! */
!/* Send intranode messages */
call MashmIntraNodeCommBegin(myMashm)
call MashmIntraNodeCommEnd(myMashm)
!/* At this stage you have completed the intra-node communication */

!/* Asynchronously do work on nodal data */

!/* Now wait on nodal messages */
call MashmInterNodeCommEnd(myMashm)

!/* Destroy the Mashm object */
call MashmDestroy(myMashm)

!decomp2dDestroyGraph(neighbors, &msgSizes)

call MPI_Finalize(ierr)

end program nodalCommFtn
