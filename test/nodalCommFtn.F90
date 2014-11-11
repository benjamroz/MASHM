module decomp2d_mod
  interface

    subroutine decomp2dCreateGraph(m, n, rank, numProcs, numElems, elements, neighbors, msgSizes, numNeighbors)
      use, intrinsic :: iso_c_binding
      integer, value :: m, n, rank, numProcs
      integer :: numElems
      type(c_ptr) :: elements
      type(c_ptr) :: neighbors
      type(c_ptr) :: msgSizes
      integer :: numNeighbors
    end subroutine

    subroutine decomp2dDestroyGraph(neighbors, msgSizes)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: neighbors
      type(c_ptr) :: msgSizes
    end subroutine

  end interface
end module

program nodalCommFtn
#include "Mashmf.h"
use iso_c_binding
use Mashm_mod
use mpi
use decomp2d_mod
implicit none
Mashm :: myMashm
integer :: ierr
integer :: i, numNeighbors
type(MashmBufferPointer), allocatable :: mashmSendBufferPtrs(:)
type(MashmBufferPointer), allocatable :: mashmRecvBufferPtrs(:)
integer, pointer :: neighbors(:), msgSizes(:)
integer :: m, n, numElems
integer :: rank, numProcs
!integer, pointer :: elements(:)
integer :: j, sumMsgSizes
integer, allocatable :: sendStatuses(:,:)
integer, allocatable :: recvStatuses(:,:)
integer, allocatable :: sendRequests(:)
integer, allocatable :: recvRequests(:)
real*8, allocatable :: origBuffer(:)
real*8, allocatable :: mashmData(:)
integer :: tag
type(c_ptr) :: cptrElements, cptrNeighbors, cptrMsgSizes
integer :: iRank
!MashmPointerArr, allocatable :: recvBuffers(:), sendBuffers(:)
type(MashmPointer1d), pointer :: recvBuffers(:), sendBuffers(:)
integer :: offset
integer, allocatable :: msgOffsets(:)
logical :: testFailed
integer :: testFailedInt, numTestsFailed

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

m = 10
n = 10

call decomp2dCreateGraph(m, n, rank, numProcs, numElems, cptrElements, cptrNeighbors, cptrMsgSizes, numNeighbors)

!call c_f_pointer(cptrElements, elements, shape=[numNeighbors])
call c_f_pointer(cptrNeighbors, neighbors, (/numNeighbors/))
call c_f_pointer(cptrMsgSizes, msgSizes, (/numNeighbors/))


allocate(msgOffsets(numNeighbors))
msgOffsets(1) = 0;
do i = 2, numNeighbors
  msgOffsets(i) = msgOffsets(i-1) + msgSizes(i-1)
enddo

sumMsgSizes = 0;
do i = 1, numNeighbors
  sumMsgSizes = sumMsgSizes + msgSizes(i)
enddo

#if 0
do iRank = 0, numProcs - 1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (rank == iRank) then
    write(*,*) "Rank ", rank
    do i = 1, numNeighbors
      write(*,*) "  msg ", i, " neighbors ", neighbors(i), " msgSizes ", msgSizes(i)
    enddo
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
enddo
#endif

!/* Allocate send and receive buffers */
allocate(recvBuffers(numNeighbors))
allocate(sendBuffers(numNeighbors))

allocate(recvRequests(numNeighbors))
allocate(sendRequests(numNeighbors))

allocate(recvStatuses(MPI_STATUS_SIZE,numNeighbors))
allocate(sendStatuses(MPI_STATUS_SIZE,numNeighbors))

!/* Symmetric */
do i = 1, numNeighbors
  allocate(recvBuffers(i)%p(msgSizes(i)))
  allocate(sendBuffers(i)%p(msgSizes(i)))
enddo

!/* Now fill the buffers */
do i = 1, numNeighbors
  do j = 1, msgSizes(i)
    sendBuffers(i)%p(j) = rank*msgSizes(i)+j;
  enddo
enddo


!/* Usual point to point communication */
tag = 0
do i = 1, numNeighbors
  call MPI_Irecv(recvBuffers(i)%p(1), msgSizes(i), MPI_REAL8, neighbors(i), tag, MPI_COMM_WORLD, recvRequests(i),ierr)
  if (ierr .ne. MPI_SUCCESS) then
    print *, "Rank ", rank, " error on irecv"
  endif
enddo

do i = 1, numNeighbors
  call MPI_Isend(sendBuffers(i)%p(1), msgSizes(i), MPI_REAL8, neighbors(i), tag, MPI_COMM_WORLD, sendRequests(i),ierr)
  if (ierr .ne. MPI_SUCCESS) then
    print *, "Rank ", rank, " error on isend"
  endif
enddo
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)

call MPI_Waitall(numNeighbors,sendRequests,sendStatuses, ierr)
if (ierr .ne. MPI_SUCCESS) then
  print *, "Rank ", rank, " error on MPI_Waitall send"
endif
call MPI_Waitall(numNeighbors,recvRequests,recvStatuses, ierr) 
if (ierr .ne. MPI_SUCCESS) then
  print *, "Rank ", rank, " error on MPI_Waitall recv"
endif

allocate(origBuffer(sumMsgSizes))
do i = 1, numNeighbors
  offset = msgOffsets(i)
  do j = 1, msgSizes(i)
    origBuffer(offset+j) = recvBuffers(i)%p(j)
  enddo
enddo

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)

call MashmInit(myMashm, MPI_COMM_WORLD)

! Print nodal comm info
call MashmPrintInfo(myMashm)

call MashmSetNumComms(myMashm, numNeighbors)

! Add communications calculated above
do i = 1, numNeighbors
  ! Fortran to C indexing 
  call MashmSetComm(myMashm, i, neighbors(i), msgSizes(i))
enddo

! Perform precalculation
call MashmCommFinish(myMashm)

!call MashmPrintCommCollection(myMashm)

! Retrieve pointers for buffers
allocate(mashmSendBufferPtrs(numNeighbors))
allocate(mashmRecvBufferPtrs(numNeighbors))

do i = 1, numNeighbors
  call MashmGetBufferPointer(myMashm, i, MASHM_SEND, mashmSendBufferPtrs(i), msgSizes(i))
  call MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE, mashmRecvBufferPtrs(i), msgSizes(i))
enddo
! Fill buffers

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
!************************************************************
! * Now perform communication 
! ************************************************************/
!  /* Fill internode buffers */
do i = 1, numNeighbors
  if (.not. MashmIsMsgOnNode(myMashm, i))  then
    do j = 1, msgSizes(i)
      mashmSendBufferPtrs(i)%p(j) = rank*msgSizes(i)+j;
    enddo
  endif
enddo

! Send internode messages
if (rank == 0) print *, "Mashm Internode communication begin."
call MashmInterNodeCommBegin(myMashm)

do i = 1, numNeighbors
  if (MashmIsMsgOnNode(myMashm, i)) then
    do j = 1, msgSizes(i)
      mashmSendBufferPtrs(i)%p(j) = rank*msgSizes(i)+j
    enddo
  endif
enddo

! Send intranode messages
if (rank == 0) print *, "Mashm Intranode communication begin."
call MashmIntraNodeCommBegin(myMashm)

allocate(mashmData(sumMsgSizes))

if (rank == 0) print *, "Mashm Intranode communication end."
call MashmIntraNodeCommEnd(myMashm)

! At this stage you have completed the intra-node communication
do i = 1, numNeighbors
  if (MashmIsMsgOnNode(myMashm, i)) then
    offset = msgOffsets(i)
    do j = 1, msgSizes(i)
      mashmData(offset+j) = mashmRecvBufferPtrs(i)%p(j)
    enddo
  endif
enddo

! Asynchronously do work on nodal data

! Now wait on nodal messages
if (rank == 0) print *, "Mashm Internode communication end."
call MashmInterNodeCommEnd(myMashm)

do i = 1, numNeighbors
  if (.not. MashmIsMsgOnNode(myMashm, i)) then
    offset = msgOffsets(i)
    do j = 1, msgSizes(i)
      mashmData(offset+j) = mashmRecvBufferPtrs(i)%p(j)
    enddo
  endif
enddo

call MPI_Barrier(MPI_COMM_WORLD, ierr)

testFailed = .false.
do i = 1, sumMsgSizes
  if (origBuffer(i) .ne. mashmData(i)) then
    print *, "Error entry ", i, " different with value ", origBuffer(i) - mashmData(i)
    testFailed = .true.
  endif
enddo

call MPI_Barrier(MPI_COMM_WORLD, ierr)

if (testFailed) then
  testFailedInt = 1
else
  testFailedInt = 0
endif

call MPI_Reduce(testFailedInt, numTestsFailed, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

if (rank == 0) then
  if (numTestsFailed .ne. 0) then
    print *, "Test Failed: buffers are different"
  else
    print *, "Test Passed: buffers are identical"
  endif
endif 

! Destroy the Mashm object
call MashmDestroy(myMashm)

! Destroy the graph data
call decomp2dDestroyGraph(cptrNeighbors, cptrMsgSizes)

call MPI_Finalize(ierr)

end program nodalCommFtn


