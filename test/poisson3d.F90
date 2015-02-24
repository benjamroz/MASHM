#define DEBUG 0
!#define PRINT_ITER
#define PI 3.1415926535897932384626

module arrayOfPointers_mod
  type arrayOfPointers
    real*8, allocatable :: p(:)
  end type arrayOfPointers
end module 


module commCycle

integer, allocatable :: sendRequest(:), recvRequest(:)
integer, allocatable :: sendStatus(:,:), recvStatus(:,:)

real*8 :: coefs(-1:1,-1:1,-1:1)
real*8 :: rhsFactor

real*8 :: coefXX, coefYY, coefZZ, coefXY, coefXZ, coefYZ

integer :: globNx, globNy, globNz
real*8 :: globDx, globDy, globDz

contains

subroutine setOperator(nu1,nu2,nu3,dx,dy,dz,angleAlpha,angleBeta)
implicit none
real*8, intent(in) :: dx, dy, dz
real*8, intent(in) :: nu1, nu2, nu3
real*8, intent(in) :: angleAlpha, angleBeta

real*8 :: tmpGrad(3,2,3,3)
real*8 :: derivXX(-1:1,-1:1,-1:1)
real*8 :: derivYY(-1:1,-1:1,-1:1)
real*8 :: derivZZ(-1:1,-1:1,-1:1)
real*8 :: derivXY(-1:1,-1:1,-1:1)
real*8 :: derivXZ(-1:1,-1:1,-1:1)
real*8 :: derivYZ(-1:1,-1:1,-1:1)
real*8 :: factor
derivXX = 0
derivXX(-1,0,0) = 1./(dx**2)
derivXX(0,0,0) = -2./(dx**2)
derivXX(1,0,0) = 1./(dx**2)

derivYY = 0
derivYY(0,-1,0) = 1./(dy**2)
derivYY(0,0,0) = -2./(dy**2)
derivYY(0,1,0) = 1./(dy**2)

derivZZ = 0
derivZZ(0,0,-1) = 1./(dz**2)
derivZZ(0,0,0) = -2./(dz**2)
derivZZ(0,0,1) = 1./(dz**2)

derivXY = 0
derivXY(-1,-1,0) = 1./(4*dx*dy)
derivXY(-1,1,0) = -1./(4*dx*dy)
derivXY(1,-1,0) = -1./(4*dx*dy)
derivXY(1,1,0) = 1./(4*dx*dy)

derivXZ = 0
derivXZ(-1,0,-1) = 1./(4*dx*dz)
derivXZ(-1,0,1) = -1./(4*dx*dz)
derivXZ(1,0,-1) = -1./(4*dx*dz)
derivXZ(1,0,1) = 1./(4*dx*dz)

derivYZ = 0
derivYZ(0,-1,-1) = 1./(4*dy*dz)
derivYZ(0,-1,1) = -1./(4*dy*dz)
derivYZ(0,1,-1) = -1./(4*dy*dz)
derivYZ(0,1,1) = 1./(4*dy*dz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The following defines the three dimensional tensor for general anisotropic
!  diffusion.
!
! L = \nabla \cdot P_{3D}^* nuT P_{3D} \nabla
! 
! where
!
!  P_{3D} = 
!   | cos(beta) 0 -sin(beta) | | cos(alpha) -sin(alpha) 0 |
!   |      0    1      0     | | sin(alpha)  cos(alpha) 0 |
!   | sin(beta) 0  cos(beta) | |     0           0      1 |
!
! and
! 
!  nuT = 
!   | nu1  0   0  |
!   |  0  nu2  0  |
!   |  0   0  nu3 |
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The coefficients for each of the quadratic derivative terms are the 
!  following
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
coefXX = nu1*dcos(angleBeta)**2*dcos(angleAlpha)**2 + &
         nu2*dsin(angleAlpha)**2 + &
         nu3*dsin(angleBeta)**2*dcos(angleAlpha)**2 
coefYY = nu1*dcos(angleBeta)**2*dsin(angleAlpha)**2 + &
         nu2*dcos(angleAlpha)**2 + &
         nu3*dsin(angleBeta)**2*dsin(angleAlpha)**2
coefZZ = nu1*dsin(angleBeta)**2 + &
         nu3*dcos(angleBeta)**2
coefXY = -2*nu1*dcos(angleBeta)**2*dcos(angleAlpha)*dsin(angleAlpha) + &
         -2*nu2*dsin(angleAlpha)*dcos(angleAlpha) + &
         -2*nu3*dsin(angleBeta)**2*dcos(angleAlpha)*dsin(angleAlpha)
coefXZ = -2*nu1*dcos(angleBeta)*dsin(angleBeta)*dcos(angleAlpha) + &
          2*nu3*dsin(angleBeta)*dcos(angleBeta)*dcos(angleAlpha)
coefYZ =  2*nu1*dcos(angleBeta)*dsin(angleBeta)*dsin(angleAlpha) + &
         -2*nu3*dsin(angleBeta)*dcos(angleBeta)*dsin(angleAlpha)

coefs = coefXX * derivXX &
      + coefYY * derivYY &
      + coefZZ * derivZZ &
      + coefXY * derivXY &
      + coefXZ * derivXZ &
      + coefYZ * derivYZ

coefs = -coefs*(dx**2 * dy**2 * dz**2)

factor = -(dx**2 * dy**2 * dz**2)/coefs(0,0,0)
! RHS factor for relaxation
rhsFactor = -factor

coefs = coefs/coefs(0,0,0)

!print *, "coefs(0,0,0) = ", coefs(0,0,0)
coefs = -coefs

! To ease the relaxation algorithm
coefs(0,0,0) = 0.0
!print *, "coefs(1,0,0) = ", coefs(1,0,0)
!print *, "coefs(-1,0,0) = ", coefs(-1,0,0)
!print *, "rhsFactor = ", rhsFactor

end subroutine

subroutine relaxation(inData, outData, rhs, gridIndicesStart, gridIndicesEnd)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(inout) :: outData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: rhs(gridIndicesStart(1):gridIndicesEnd(1), &
                          gridIndicesStart(2):gridIndicesEnd(2), &
                          gridIndicesStart(3):gridIndicesEnd(3))

integer :: ix, iy, iz
integer :: i, j, k
real*8 :: relaxOmega

relaxOmega = 1.0
!relaxOmega = 2.0/3.0

do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      outData(ix,iy,iz) = (1.0-relaxOmega)*inData(ix,iy,iz) + relaxOmega*rhsFactor*rhs(ix,iy,iz)
      do k = -1, 1
        do j = -1, 1
          do i = -1, 1
            outData(ix,iy,iz) = outData(ix,iy,iz) + relaxOmega*coefs(i,j,k)*inData(ix+i,iy+j,iz+k)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
end subroutine

subroutine communication(sendBuffer, recvBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)
use mpi
use gptl
implicit none
real*8, intent(in) :: sendBuffer(:)
real*8, intent(out) :: recvBuffer(:)
integer, intent(in) :: numMessages
integer, intent(in) :: neighborRanks(:)
integer, intent(in) :: msgSizes(:)
integer, intent(in) :: msgOffsets(:)

integer :: i

integer :: ierr, errorCode, errorLen
character(256) :: errorString
integer :: gptlError

do i = 1, numMessages
  call MPI_Irecv(recvBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,recvRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
     errorcode = ierr
     call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
     print *,'Error after call to MPI_Irecv: ',errorstring
  endif
end do

do i = 1, numMessages
  call MPI_Isend(sendBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,sendRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
    errorcode = ierr
    call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
    print *,'Error after call to MPI_Isend: ',errorstring
  endif
end do

gptlError = gptlstart('communication1_waitall')
call MPI_Waitall(numMessages,recvRequest,recvStatus,ierr)
call MPI_Waitall(numMessages,sendRequest,sendStatus,ierr)
gptlError = gptlstop('communication1_waitall')

end subroutine 

subroutine communication_notimers(sendBuffer, recvBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)
use mpi
use gptl
implicit none
real*8, intent(in) :: sendBuffer(:)
real*8, intent(out) :: recvBuffer(:)
integer, intent(in) :: numMessages
integer, intent(in) :: neighborRanks(:)
integer, intent(in) :: msgSizes(:)
integer, intent(in) :: msgOffsets(:)

integer :: i

integer :: ierr, errorCode, errorLen
character(256) :: errorString
integer :: gptlError

do i = 1, numMessages
  call MPI_Irecv(recvBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,recvRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
     errorcode = ierr
     call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
     print *,'Error after call to MPI_Irecv: ',errorstring
  endif
end do

do i = 1, numMessages
  call MPI_Isend(sendBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,sendRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
    errorcode = ierr
    call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
    print *,'Error after call to MPI_Isend: ',errorstring
  endif
end do

call MPI_Waitall(numMessages,recvRequest,recvStatus,ierr)
call MPI_Waitall(numMessages,sendRequest,sendStatus,ierr)

end subroutine 


subroutine setupComm(numMessages)
use mpi
implicit none
integer, intent(in) :: numMessages
allocate(sendRequest(numMessages))
allocate(recvRequest(numMessages))
allocate(sendStatus(MPI_STATUS_SIZE,numMessages))
allocate(recvStatus(MPI_STATUS_SIZE,numMessages))
end subroutine

subroutine setSolution(solution, gridIndicesStart, gridIndicesEnd, dx, dy, dz, &
                       nu1, nu2, nu3)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(out) :: solution(gridIndicesStart(1):gridIndicesEnd(1), &
                                gridIndicesStart(2):gridIndicesEnd(2), &
                                gridIndicesStart(3):gridIndicesEnd(3))
real*8, intent(in) :: dx, dy, dz
real*8, intent(in) :: nu1, nu2, nu3

integer :: i, j, k
real*8 :: x, y, z

do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      x = i*dx
      y = j*dy
      z = k*dz
      !solution(i,j,k) = sin(x)*sin(y)*sin(z)
      solution(i,j,k) = 0.0
    enddo
  enddo
enddo
end subroutine 

subroutine setRhs(rhs, gridIndicesStart, gridIndicesEnd, dx, dy, dz, &
                  nu1, nu2, nu3)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(out) :: rhs(gridIndicesStart(1):gridIndicesEnd(1), &
                           gridIndicesStart(2):gridIndicesEnd(2), &
                           gridIndicesStart(3):gridIndicesEnd(3))
real*8, intent(in) :: dx, dy, dz
real*8, intent(in) :: nu1, nu2, nu3

integer :: i, j, k
real*8 :: x, y, z

do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      x = i*dx
      y = j*dy
      z = k*dz
      !rhs(i,j,k) = 3.0*sin(x)*sin(y)*sin(z)
      !rhs(i,j,k) = (nu1 + nu2 + nu3)*sin(x)*sin(y)*sin(z)
      rhs(i,j,k) = 0.0
    enddo
  enddo
enddo

end subroutine

subroutine calcL2Norm(numSoln, solution, gridIndicesStart, gridIndicesEnd, globalL2Norm, globalMaxNorm, totalNumCells)
use mpi
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: numSoln(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                              gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                              gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: solution(gridIndicesStart(1):gridIndicesEnd(1), &
                               gridIndicesStart(2):gridIndicesEnd(2), &
                               gridIndicesStart(3):gridIndicesEnd(3))
real*8, intent(out) :: globalL2Norm, globalMaxNorm
integer, intent(in) :: totalNumCells

real*8 :: localL2Norm, localMaxNorm
real*8 :: newMaxNorm

integer :: i, j, k
integer :: ierr

localL2Norm = 0.0
localMaxNorm = 0.0
do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      localL2Norm = localL2Norm + (numSoln(i,j,k) - solution(i,j,k))**2
      localMaxNorm = maxval( (/ localMaxNorm, dabs( numSoln(i,j,k) - solution(i,j,k) ) /) )
    enddo
  enddo
enddo

! MPI_Allreduce 
call MPI_Allreduce(localL2Norm, globalL2Norm, 1, MPI_REAL8, MPI_SUM, &
                   MPI_COMM_WORLD, ierr)
call MPI_Allreduce(localMaxNorm, globalMaxNorm, 1, MPI_REAL8, MPI_MAX, &
                   MPI_COMM_WORLD, ierr)

globalL2Norm = dsqrt(globalL2Norm/totalNumCells)
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
use Mashm_type
use gptl
use poisson3dPack
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
real*8, allocatable :: tmpDomain(:,:,:)
real*8, allocatable :: solution(:,:,:)
real*8, allocatable :: rhs(:,:,:)
integer :: totalMessageSize
real*8, allocatable :: packBuffer(:), unpackBuffer(:)
integer, allocatable :: msgOffsets(:)
integer :: offsetValue, i, j, k
integer :: msgDirIndex(-1:1,-1:1,-1:1)

integer :: nx, ny, nz
real*8 :: dx, dy, dz
real*8 :: angleAlpha, angleBeta
real*8 :: nu1, nu2, nu3
integer :: iIter, numIters
real*8 :: residualL2, residualMax
integer :: totalNumCells
integer :: counter 
type(Mashm) :: myMashm
!MashmCommType :: commMethod
integer(kind=c_int) :: commMethod
type(MashmBufferPointer), allocatable :: mashmSendBufferPtrs(:)
type(MashmBufferPointer), allocatable :: mashmRecvBufferPtrs(:)

integer, allocatable :: msgDirIndex2(:)
integer :: msgIndex
integer :: packInt3(3)

integer :: gptlError, gptlRet

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

! Read namelist
call read_grid_data_namelist(MPI_COMM_WORLD)

! Set gptl to nanotime
gptlError = gptlsetutr(gptlnanotime)

! Initialize gptl
gptlError = gptlinitialize()

! Get the grid decomposition
call grid_3d_decomp_num_elements(rank, numProcs)

call grid_3d_get_indices(numElems, gridIndicesStart, gridIndicesEnd)


call determineCommSchedule(rank, numProcs, gridIndicesStart, gridIndicesEnd, numMessages, msgDirIndex, &
                           msgSizes, msgOffsets, neighborRanks)
call setupComm(numMessages)

!print *, "rank ", rank, ", numMessages = ", numMessages

allocate(msgDirIndex2(numMessages))
do k = -1,1
  do j = -1,1
    do i = -1,1
      msgIndex = msgDirIndex(i,j,k)
      if (msgIndex > 0) then
        packInt3(1) = i
        packInt3(2) = j
        packInt3(3) = k
        msgDirIndex2(msgIndex) = packDirInt( packInt3 )
      endif
    enddo
  enddo
enddo


totalMessageSize=msgOffsets(numMessages+1)

allocate(packBuffer(totalMessageSize))
allocate(unpackBuffer(totalMessageSize))

! Initialize domain data
allocate(domain(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                gridIndicesStart(3)-1:gridIndicesEnd(3)+1))

allocate(tmpDomain(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                gridIndicesStart(3)-1:gridIndicesEnd(3)+1))

allocate(solution(gridIndicesStart(1):gridIndicesEnd(1), &
                  gridIndicesStart(2):gridIndicesEnd(2), &
                  gridIndicesStart(3):gridIndicesEnd(3)))

allocate(rhs(gridIndicesStart(1):gridIndicesEnd(1), &
             gridIndicesStart(2):gridIndicesEnd(2), &
             gridIndicesStart(3):gridIndicesEnd(3)))

nx = get_num_cells(1)
ny = get_num_cells(2)
nz = get_num_cells(3)
totalNumCells = nx*ny*nz

globNx = nx
globNy = ny
globNz = nz

dx = get_grid_length(1)/(nx + 1)
dy = get_grid_length(2)/(ny + 1)
dz = get_grid_length(3)/(nz + 1)

globDx = dx
globDy = dy
globDz = dz

angleAlpha = 3.14/8.0
angleBeta = 3.14/8.0

nu1 = 1.0
nu2 = 2.0
nu3 = 5.0

call setRhs(rhs,gridIndicesStart,gridIndicesEnd,dx,dy,dz,nu1,nu2,nu3)
call setSolution(solution,gridIndicesStart,gridIndicesEnd,dx,dy,dz,nu1,nu2,nu3)

call setOperator(nu1,nu2,nu3,dx,dy,dz,angleAlpha,angleBeta)

!if (rank == 0) print *, "coefs = ", coefs

! Set the numerical solution arrays to zero
!   Note that this will provide Dirichlet zero boundary conditions
domain = 0.0d0
tmpDomain = 0.0d0
do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      domain(i,j,k) = dsin(10*i*dx)*dsin(10*j*dy)*dsin(10*k*dz)
    enddo
  enddo
enddo

call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)

if (rank == 0) print *, "Initial difference", residualL2, residualMax

numIters = 10000


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  "Warmup"
!  Perform 2 cycles to ensure communication initialization is not timed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iIter = 1, 2, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, packBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call communication_notimers(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, unpackBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  do i = 1, numMessages
    call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, packBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call communication_notimers(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)

  do i = 1, numMessages
    call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, unpackBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MPI_Barrier(MPI_COMM_WORLD, ierr)
gptlError = gptlstart('method1')

do iIter = 1, numIters, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, packBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  gptlError = gptlstart('communication1')
  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)
  gptlError = gptlstop('communication1')

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, unpackBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

#ifdef PRINT_ITER
  call calcL2Norm(tmpDomain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)
  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
#endif

  do i = 1, numMessages
    call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, packBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  gptlError = gptlstart('communication1')
  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)
  gptlError = gptlstop('communication1')


  do i = 1, numMessages
    call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, unpackBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

#ifdef PRINT_ITER
  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)
  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax
#endif

enddo
gptlError = gptlstop('method1')

#ifndef PRINT_ITER
call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)
if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                        residualMax
#endif

call MashmInit(myMashm, MPI_COMM_WORLD)

! Print nodal comm info
!call MashmPrintInfo(myMashm)

call MashmSetNumComms(myMashm, numMessages)

! Add communications calculated above
do i = 1, numMessages
  ! Fortran to C indexing 
  call MashmSetComm(myMashm, i, neighborRanks(i), msgSizes(i))
enddo

! Choose the communitcation method
commMethod = MASHM_COMM_MIN_AGG

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
if (rank == 0) print *, "Calling finish"
flush(6)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)

call MashmSetCommMethod(myMashm, commMethod)

call MashmSetCacheBlocking(myMashm, .true.)

! Perform precalculation
call MashmCommFinish(myMashm)

! Print nodal messages statistics
call MashmPrintMessageStats(myMashm)
!call MashmPrintCommCollection(myMashm)

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
if (rank == 0) print *, "Returned from finish"
flush(6)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
! Retrieve pointers for buffers
allocate(mashmSendBufferPtrs(numMessages))
allocate(mashmRecvBufferPtrs(numMessages))

do i = 1, numMessages
  call MashmGetBufferPointer(myMashm, i, MASHM_SEND, mashmSendBufferPtrs(i))
  call MashmGetBufferPointer(myMashm, i, MASHM_RECEIVE, mashmRecvBufferPtrs(i))
enddo
! Reset the solution
tmpDomain = 0.0
do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      domain(i,j,k) = dsin(10*i*dx)*dsin(10*j*dy)*dsin(10*k*dz)
    enddo
  enddo
enddo

call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)

if (rank == 0) print *, "Now running MASHM"
if (rank == 0) print *, "Initial difference", residualL2, residualMax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  "Warmup"
!  Perform 2 cycles to ensure communication initialization is not timed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iIter = 1, 2, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call MashmInterNodeCommBegin(myMashm)
  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)
  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  do i = 1, numMessages
    call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call MashmInterNodeCommBegin(myMashm)
  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)
  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call MPI_Barrier(MPI_COMM_WORLD, ierr)
gptlError = gptlstart('method2')
do iIter = 1, numIters, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  gptlError = gptlstart('communication2')
  call MashmInterNodeCommBegin(myMashm)
  gptlError = gptlstart('communication2_intra')
  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)
  gptlError = gptlstop('communication2_intra')
  gptlError = gptlstart('communication2_waitall')
  call MashmInterNodeCommEnd(myMashm)
  gptlError = gptlstop('communication2_waitall')
  gptlError = gptlstop('communication2')

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

#ifdef PRINT_ITER
  call calcL2Norm(tmpDomain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
#endif

  do i = 1, numMessages
    call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  gptlError = gptlstart('communication2')
  call MashmInterNodeCommBegin(myMashm)
  gptlError = gptlstart('communication2_intra')
  call MashmIntraNodeCommBegin(myMashm)
  call MashmIntraNodeCommEnd(myMashm)
  gptlError = gptlstop('communication2_intra')

  gptlError = gptlstart('communication2_waitall')
  call MashmInterNodeCommEnd(myMashm)
  gptlError = gptlstop('communication2_waitall')
  gptlError = gptlstop('communication2')
  do i = 1, numMessages
    call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

#ifdef PRINT_ITER
  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax
#endif

enddo
gptlError = gptlstop('method2')

#ifndef PRINT_ITER
call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)
if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                        residualMax
#endif

! Reset the solution
tmpDomain = 0.0
do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      domain(i,j,k) = dsin(10*i*dx)*dsin(10*j*dy)*dsin(10*k*dz)
    enddo
  enddo
enddo

call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)

if (rank == 0) print *, "Now running MASHM with asynchronous packing/unpacking"
if (rank == 0) print *, "Initial difference", residualL2, residualMax


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  "Warmup"
!  Perform 2 cycles to ensure communication initialization is not timed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iIter = 1, 2, 2

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommBegin(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmIntraNodeCommBegin(myMashm)
  call MashmIntraNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommBegin(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Stride two 
call MPI_Barrier(MPI_COMM_WORLD, ierr)
gptlError = gptlstart('method3')
do iIter = 1, numIters, 2

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  gptlError = gptlstart('communication3_inter')
  gptlError = gptlstart('communication3')
  call MashmInterNodeCommBegin(myMashm)
  gptlError = gptlstop('communication3')
  gptlError = gptlstop('communication3_inter')

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  gptlError = gptlstart('communication3_intra')
  gptlError = gptlstart('communication3')
  call MashmIntraNodeCommBegin(myMashm)
  call MashmIntraNodeCommEnd(myMashm)
  gptlError = gptlstop('communication3')
  gptlError = gptlstop('communication3_intra')

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  gptlError = gptlstart('communication3_waitall')
  gptlError = gptlstart('communication3_inter')
  gptlError = gptlstart('communication3')
  call MashmInterNodeCommEnd(myMashm)
  gptlError = gptlstop('communication3')
  gptlError = gptlstop('communication3_inter')
  gptlError = gptlstop('communication3_waitall')

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

#ifdef PRINT_ITER
  call calcL2Norm(tmpDomain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)
  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
#endif

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  gptlError = gptlstart('communication3_inter')
  gptlError = gptlstart('communication3')
  call MashmInterNodeCommBegin(myMashm)
  gptlError = gptlstop('communication3')
  gptlError = gptlstop('communication3_inter')

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  gptlError = gptlstart('communication3_intra')
  gptlError = gptlstart('communication3')
  call MashmIntraNodeCommBegin(myMashm)
  call MashmIntraNodeCommEnd(myMashm)
  gptlError = gptlstop('communication3')
  gptlError = gptlstop('communication3_intra')

  do i = 1, numMessages
    if (MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  gptlError = gptlstart('communication3_waitall')
  gptlError = gptlstart('communication3_inter')
  gptlError = gptlstart('communication3')
  call MashmInterNodeCommEnd(myMashm)
  gptlError = gptlstop('communication3')
  gptlError = gptlstop('communication3_inter')
  gptlError = gptlstop('communication3_waitall')

  do i = 1, numMessages
    if (.not. MashmIsMsgIntraNodal(myMashm, i)) then
      call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

#ifdef PRINT_ITER
  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)
  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax
#endif

enddo
gptlError = gptlstop('method3')

#ifndef PRINT_ITER
call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)
if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                        residualMax
#endif

gptlError = gptlpr_summary_file(MPI_COMM_WORLD,'poisson3d.timing')
!gptlError = gptlpr_file('poisson3d.timing')
! Restore (nullify) the Mashm access pointers
do i = 1, numMessages
  call MashmRetireBufferPointer(myMashm, mashmSendBufferPtrs(i))
  call MashmRetireBufferPointer(myMashm, mashmRecvBufferPtrs(i))
enddo

! Destroy the Mashm object
call MashmDestroy(myMashm)


deallocate(packBuffer)
deallocate(unpackBuffer)

deallocate(domain)
deallocate(tmpDomain)
deallocate(solution)
deallocate(rhs)
call MPI_Finalize(ierr)

end program nodalCommFtn

