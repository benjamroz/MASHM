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
integer :: j, sumMsgSizes, counter
integer, allocatable :: sendStatuses(:,:)
integer, allocatable :: recvStatuses(:,:)
integer, allocatable :: sendRequests(:)
integer, allocatable :: recvRequests(:)
real*8, allocatable :: origBuffer(:)
integer :: tag
type(c_ptr) :: cptrElements, cptrNeighbors, cptrMsgSizes
integer :: iRank
!MashmPointerArr, allocatable :: recvBuffers(:), sendBuffers(:)
type(MashmPointer1d), pointer :: recvBuffers(:), sendBuffers(:)

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

m = 10
n = 10

call decomp2dCreateGraph(m, n, rank, numProcs, numElems, cptrElements, cptrNeighbors, cptrMsgSizes, numNeighbors)

!call c_f_pointer(cptrElements, elements, shape=[numNeighbors])
call c_f_pointer(cptrNeighbors, neighbors, (/numNeighbors/))
call c_f_pointer(cptrMsgSizes, msgSizes, (/numNeighbors/))


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

sumMsgSizes = 0
do i = 1, numNeighbors
  sumMsgSizes = sumMsgSizes + msgSizes(i) 
enddo

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
    !sendBuffers(i)%p(j) = rank*msgSizes(i)+j;
    sendBuffers(i)%p(j) = -(rank+1)
  enddo
enddo

#if 0
do iRank = 0, numProcs - 1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (rank == iRank) then
    write(*,*) "Rank ", rank
    do i = 1, numNeighbors
      write(*,*) "  send msg ", i, " data ", sendBuffers(i)%p
    enddo
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call flush(6)
enddo
#endif

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

#if 0
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)
do iRank = 0, numProcs - 1
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  if (rank == iRank) then
    write(*,*) "Rank ", rank
    do i = 1, numNeighbors
      write(*,*) "  receive msg ", i, " data ", recvBuffers(i)%p
    enddo
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call flush(6)
enddo
#endif


allocate(origBuffer(sumMsgSizes))
print *, "origBuffer = "
counter = 1
do i = 1, numNeighbors
  do j = 1, msgSizes(i)
    origBuffer(counter) = recvBuffers(i)%p(j)
    print *, "  ", origBuffer(counter)
    counter = counter + 1
  enddo
enddo


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

call MashmPrintCommCollection(myMashm)

! Retrieve pointers for buffers
allocate(mashmSendBufferPtrs(numNeighbors))
allocate(mashmRecvBufferPtrs(numNeighbors))

do i = 1, numNeighbors
  call MashmGetBufferPointer(myMashm, i, MASHM_SEND, mashmSendBufferPtrs(i)) 
  call MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE, mashmRecvBufferPtrs(i)) 
enddo

! Fill buffers

!************************************************************
! * Now perform communication 
! ************************************************************/

! Send internode messages
call MashmInterNodeCommBegin(myMashm)

! Messages sent and receives posted 
! * Can asynchronously do work on nodal data 
!
! Send intranode messages
call MashmIntraNodeCommBegin(myMashm)
call MashmIntraNodeCommEnd(myMashm)
! At this stage you have completed the intra-node communication

! Asynchronously do work on nodal data

! Now wait on nodal messages
call MashmInterNodeCommEnd(myMashm)

! Destroy the Mashm object
call MashmDestroy(myMashm)

! Destroy the graph data
call decomp2dDestroyGraph(cptrNeighbors, cptrMsgSizes)

call MPI_Finalize(ierr)
end program nodalCommFtn


