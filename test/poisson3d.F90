module arrayOfPointers_mod
  type arrayOfPointers
    real*8, allocatable :: p(:)
  end type arrayOfPointers
end module 


module commCycle

contains

function evalOperator(inData, ix, iy, iz) result(outData)
implicit none
real*8, intent(in) :: inData(:,:,:)
integer, intent(in) :: ix, iy, iz

real*8 :: outData

integer :: i, j, k
real*8 :: coefs(-1:1,-1:1,-1:1)

! TODO change this to 27 point stencil
coefs(0,0,-1) = 1
coefs(0,-1,0) = 1
coefs(-1,0,0) = 1
coefs(0,0,0) = -8
coefs(1,0,0) = 1
coefs(0,1,0) = 1
coefs(0,0,1) = 1
outData = 0
do k = -1, 1
  do j = -1, 1
    do i = -1, 1
      outData = outData + coefs(i,j,k)*inData(i,j,k)
    enddo
  enddo
enddo

end function


subroutine haloExchange

end subroutine

subroutine relaxation(inData, outData, gridIndicesStart, gridIndicesEnd)
implicit none
real*8, intent(out) :: outData(:,:,:)
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)

integer :: ix, iy, iz

do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      outData(ix, iy, iz) = evalOperator(inData, ix, iy, iz)
    enddo
  enddo
enddo
end subroutine

subroutine packData(inData, gridIndicesStart, gridIndicesEnd,  &
                    msgOffsets, msgLengths, packBuffer, packDir)
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
integer, intent(in) :: msgOffsets(:), msgLengths(:)
real*8, intent(out) :: packBuffer(:)
integer, intent(in) :: packDir(-1:1,-1:1,-1:1)


integer :: packIndex, packOffset

! pack top
packIndex = packDir(0,0,1)
if (packIndex > 0) then
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      packBuffer(packOffset) = inData(ix,iy,iz)
      packOffset = packOffset + 1
    enddo
  enddo
endif


! pack bottom
packIndex = packDir(0,0,-1)
if (packIndex > 0) then
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      packBuffer(packOffset) = inData(ix,iy,iz)
      packOffset = packOffset + 1
    enddo
  enddo
endif


! pack south
packIndex = packDir(0,-1,0)
if (packIndex > 0) then
  iy = gridIndicesStart(2)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      packBuffer(packOffset) = inData(ix,iy,iz)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack north
packIndex = packDir(0,1,0)
if (packIndex > 0) then
  iy = gridIndicesEnd(2)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      packBuffer(packOffset) = inData(ix,iy,iz)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack east
packIndex = packDir(-1,0,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do iy = gridIndicesStart(2), gridIndicesEnd(2)
      packBuffer(packOffset) = inData(ix,iy,iz)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack west
packIndex = packDir(1,0,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do iy = gridIndicesStart(2), gridIndicesEnd(2)
      packBuffer(packOffset) = inData(ix,iy,iz)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! Pack bottom edges
! Bottom east
packIndex = packDir(-1,0,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Bottom west
packIndex = packDir(1,0,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Bottom south
packIndex = packDir(0,-1,-1)
if (packIndex > 0) then
  iy = gridIndicesStart(2)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Bottom north
packIndex = packDir(0,1,-1)
if (packIndex > 0) then
  iy = gridIndicesEnd(2)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif


! Pack top edges
! Bottom east
packIndex = packDir(-1,0,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Bottom west
packIndex = packDir(1,0,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Bottom south
packIndex = packDir(0,-1,1)
if (packIndex > 0) then
  iy = gridIndicesStart(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Bottom north
packIndex = packDir(0,1,1)
if (packIndex > 0) then
  iy = gridIndicesEnd(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    packBuffer(packOffset) = inData(ix,iy,iz)
    packOffset = packOffset + 1
  enddo
endif

! Pack corners
! Bottom south east
packIndex = packDir(-1,-1,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iy = gridIndicesStart(2)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

packIndex = packDir(1,-1,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iy = gridIndicesStart(2)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

packIndex = packDir(-1,1,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iy = gridIndicesEnd(2)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

packIndex = packDir(1,1,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iy = gridIndicesEnd(2)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

! Top edges
packIndex = packDir(-1,-1,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iy = gridIndicesStart(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

packIndex = packDir(1,-1,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iy = gridIndicesStart(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

packIndex = packDir(-1,1,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iy = gridIndicesEnd(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

packIndex = packDir(1,1,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iy = gridIndicesEnd(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  packBuffer(packOffset) = inData(ix,iy,iz)
endif

end subroutine

subroutine unpackData(inData, gridIndicesStart, gridIndicesEnd,  &
                    msgOffsets, msgLengths, unpackBuffer, packDir)
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                                gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                               gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
integer, intent(in) :: msgOffsets(:), msgLengths(:)
real*8, intent(in) :: unpackBuffer(:)
integer, intent(in) :: packDir(-1:1,-1:1,-1:1)


integer :: packIndex, packOffset

! pack top
packIndex = packDir(0,0,1)
if (packIndex > 0) then
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      !packBuffer(packOffset) = inData(ix,iy,iz)
      inData(ix,iy,iz) = unpackBuffer(packOffset)
      packOffset = packOffset + 1
    enddo
  enddo
endif


! pack bottom
packIndex = packDir(0,0,-1)
if (packIndex > 0) then
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      !packBuffer(packOffset) = inData(ix,iy,iz)
      inData(ix,iy,iz) = unpackBuffer(packOffset)
      packOffset = packOffset + 1
    enddo
  enddo
endif


! pack south
packIndex = packDir(0,-1,0)
if (packIndex > 0) then
  iy = gridIndicesStart(2) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      !packBuffer(packOffset) = inData(ix,iy,iz)
      inData(ix,iy,iz) = unpackBuffer(packOffset)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack north
packIndex = packDir(0,1,0)
if (packIndex > 0) then
  iy = gridIndicesEnd(2) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      !packBuffer(packOffset) = inData(ix,iy,iz)
      inData(ix,iy,iz) = unpackBuffer(packOffset)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack east
packIndex = packDir(-1,0,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do iy = gridIndicesStart(2), gridIndicesEnd(2)
      !packBuffer(packOffset) = inData(ix,iy,iz)
      inData(ix,iy,iz) = unpackBuffer(packOffset)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack west
packIndex = packDir(1,0,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do iy = gridIndicesStart(2), gridIndicesEnd(2)
      !packBuffer(packOffset) = inData(ix,iy,iz)
      inData(ix,iy,iz) = unpackBuffer(packOffset)
      packOffset = packOffset + 1
    enddo
  enddo
endif

! Pack bottom edges
! Bottom east
packIndex = packDir(-1,0,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Bottom west
packIndex = packDir(1,0,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Bottom south
packIndex = packDir(0,-1,-1)
if (packIndex > 0) then
  iy = gridIndicesStart(2) - 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Bottom north
packIndex = packDir(0,1,-1)
if (packIndex > 0) then
  iy = gridIndicesEnd(2) + 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif


! Pack top edges
! Bottom east
packIndex = packDir(-1,0,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Bottom west
packIndex = packDir(1,0,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Bottom south
packIndex = packDir(0,-1,1)
if (packIndex > 0) then
  iy = gridIndicesStart(2) - 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Bottom north
packIndex = packDir(0,1,1)
if (packIndex > 0) then
  iy = gridIndicesEnd(2) + 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    !packBuffer(packOffset) = inData(ix,iy,iz)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
endif

! Pack corners
! Bottom south east
packIndex = packDir(-1,-1,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iy = gridIndicesStart(2) - 1 
  iz = gridIndicesStart(3) - 1 
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

packIndex = packDir(1,-1,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iy = gridIndicesStart(2) - 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

packIndex = packDir(-1,1,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iy = gridIndicesEnd(2) + 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

packIndex = packDir(1,1,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iy = gridIndicesEnd(2) + 1
  iz = gridIndicesStart(3) - 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

! Top corners
packIndex = packDir(-1,-1,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iy = gridIndicesStart(2) - 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

packIndex = packDir(1,-1,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iy = gridIndicesStart(2) - 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

packIndex = packDir(-1,1,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iy = gridIndicesEnd(2) + 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif

packIndex = packDir(1,1,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iy = gridIndicesEnd(2) + 1
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  !packBuffer(packOffset) = inData(ix,iy,iz)
  inData(ix,iy,iz) = unpackBuffer(packOffset)
endif


end subroutine

end module commCycle

program nodalCommFtn
#include "Mashmf.h"
use iso_c_binding
use Mashm_mod
use mpi
use arrayOfPointers_mod
use grid_data
use commCycle
implicit none

integer :: ierr
integer :: rank, numProcs
integer :: numElems
integer, pointer :: gridIndicesStart(:), gridIndicesEnd(:)
integer :: ix, iy, iz, iRank
integer :: elemRank, rankCounter
integer, allocatable :: elemRankCounter(:)
integer :: numMessages
integer, allocatable :: neighborRanks(:), msgSizes(:)
logical :: found
real*8, allocatable :: domain(:,:,:)
real*8, allocatable :: outData(:,:,:)
integer :: totalMessageSize
real*8, allocatable :: packBuffer(:), unpackBuffer(:)
integer, allocatable :: msgOffsets(:)
integer :: offsetValue, i, j, k
integer :: msgDirIndex(-1:1,-1:1,-1:1)

integer :: errorCode, errorLen
character(256) :: errorString
integer, allocatable :: sendRequest(:), recvRequest(:)
integer, allocatable :: sendStatus(:,:), recvStatus(:,:)

Mashm :: myMashm

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

! Read namelist
call read_grid_data_namelist(MPI_COMM_WORLD)


! Get the grid decomposition
call grid_3d_decomp_num_elements(rank, numProcs)

call grid_3d_get_indices(numElems, gridIndicesStart, gridIndicesEnd)


call determineCommSchedule(rank, numProcs, gridIndicesStart, gridIndicesEnd, numMessages, msgDirIndex, &
                           msgSizes, msgOffsets, neighborRanks)

allocate(sendRequest(numMessages))
allocate(recvRequest(numMessages))
allocate(sendStatus(MPI_STATUS_SIZE,numMessages))
allocate(recvStatus(MPI_STATUS_SIZE,numMessages))

totalMessageSize=msgOffsets(numMessages+1)
allocate(packBuffer(totalMessageSize))
allocate(unpackBuffer(totalMessageSize))

! Initialize domain data
allocate(domain(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                gridIndicesStart(3)-1:gridIndicesEnd(3)+1))

! Initialize to zero
domain = 0
offsetValue = rank &
            * (gridIndicesEnd(3)-gridIndicesStart(3)+2) &
            * (gridIndicesEnd(2)-gridIndicesStart(2)+2) &
            * (gridIndicesEnd(1)-gridIndicesStart(1)+2)
do k = gridIndicesStart(3)-1, gridIndicesEnd(3)+1
  do j = gridIndicesStart(2)-1, gridIndicesEnd(2)+1
    do i = gridIndicesStart(1)-1, gridIndicesEnd(1)+1
      domain(i,j,k) = offsetValue &
                      + k*(gridIndicesEnd(2)-gridIndicesStart(2)+2)*(gridIndicesEnd(1)-gridIndicesStart(1)+2) &
                      + j*(gridIndicesEnd(1)-gridIndicesStart(1)+2) + i
    enddo
  enddo
enddo

! The pack is incorrect
call packData(domain, gridIndicesStart, gridIndicesEnd,  &
              msgOffsets, msgSizes, packBuffer, msgDirIndex)

! Need messaging loop
!  Should have all of the data that we need
!  Allocate a receive buffer - pass that into unpack

do i = 1, numMessages
  call MPI_Isend(packBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,sendRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
    errorcode = ierr
    call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
    print *,'Error after call to MPI_Isend: ',errorstring
  endif
end do    ! icycle

!==================================================
!  Post the Receives
!==================================================
do i = 1, numMessages
  call MPI_Irecv(unpackBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,recvRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
     errorcode = ierr
     call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
     print *,'Error after call to MPI_Irecv: ',errorstring
  endif
end do    ! icycle


call MPI_Waitall(numMessages,sendRequest,sendStatus,ierr)
call MPI_Waitall(numMessages,recvRequest,recvStatus,ierr)

call unpackData(domain, gridIndicesStart, gridIndicesEnd,  &
              msgOffsets, msgSizes, unpackBuffer, msgDirIndex)


call MPI_Finalize(ierr)

end program nodalCommFtn

