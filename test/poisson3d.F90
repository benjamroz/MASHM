module arrayOfPointers_mod
  type arrayOfPointers
    real*8, allocatable :: p(:)
  end type arrayOfPointers
end module 


module commCycle

integer, allocatable :: sendRequest(:), recvRequest(:)
integer, allocatable :: sendStatus(:,:), recvStatus(:,:)

real*8 :: coefs(-1:1,-1:1,-1:1)

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
derivXY(-1,-1,0) = 1./(dx*dy)
derivXY(-1,1,0) = -1./(dx*dy)
derivXY(1,-1,0) = -1./(dx*dy)
derivXY(1,1,0) = 1./(dx*dy)

derivXZ = 0
derivXZ(-1,0,-1) = 1./(dx*dz)
derivXZ(-1,0,1) = -1./(dx*dz)
derivXZ(1,0,-1) = -1./(dx*dz)
derivXZ(1,0,1) = 1./(dx*dz)

derivYZ = 0
derivYZ(0,-1,-1) = 1./(dy*dz)
derivYZ(0,-1,1) = -1./(dy*dz)
derivYZ(0,1,-1) = -1./(dy*dz)
derivYZ(0,1,1) = 1./(dy*dz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       
!      | nu1  0   0  | | cos(beta) 0 -sin(beta) | | cos(alpha) -sin(alpha) 0 |
!  K = |  0  nu2  0  | |      0    1      0     | | sin(alpha)  cos(alpha) 0 |
!      |  0   0  nu3 | | sin(beta) 0  cos(beta) | |     0           0      1 |
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
coefs = nu1 * dcos(angleBeta) * dcos(angleAlpha) * derivXX &
      + nu2 * dcos(angleAlpha) * derivYY &
      + nu3 * dcos(angleBeta) * derivZZ &
      + (nu2 * dsin(angleAlpha) - nu1 * dcos(angleBeta)*dcos(angleAlpha)) * derivXY &
      + (nu3 * dsin(angleBeta)*dcos(angleAlpha) - nu1 * dsin(angleBeta)) * derivXZ &
      + (-nu3 * dsin(angleBeta)*dsin(angleAlpha)) * derivYZ
! \nabla \cdot : [-1, 1

end subroutine

function evalOperator(inData, ix, iy, iz) result(outData)
implicit none
real*8, intent(in) :: inData(-1:1,-1:1,-1:1)
integer, intent(in) :: ix, iy, iz

real*8 :: outData

integer :: i, j, k

! TODO change this to 27 point stencil
outData = 0
do k = -1, 1
  do j = -1, 1
    do i = -1, 1
      outData = outData + coefs(i,j,k)*inData(i,j,k)
    enddo
  enddo
enddo

end function

subroutine relaxation(inData, outData, gridIndicesStart, gridIndicesEnd)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: outData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)

integer :: ix, iy, iz
integer :: i, j, k
real*8 :: relaxDt

relaxDt = 0.1

do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      outData(ix,iy,iz) = inData(ix,iy,iz)
      do k = -1, 1
        do j = -1, 1
          do i = -1, 1
            outData(ix,iy,iz) = outData(ix,iy,iz) + relaxDt*coefs(i,j,k)*inData(ix+i,iy+j,iz+k)
          enddo
        enddo
      enddo
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

subroutine communication(sendBuffer, recvBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)
use mpi
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

do i = 1, numMessages
  call MPI_Isend(sendBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,sendRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
    errorcode = ierr
    call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
    print *,'Error after call to MPI_Isend: ',errorstring
  endif
end do

do i = 1, numMessages
  call MPI_Irecv(recvBuffer(msgOffsets(i)+1),msgSizes(i),MPI_REAL8,neighborRanks(i),10,MPI_COMM_WORLD,recvRequest(i),ierr)
  if(ierr .ne. MPI_SUCCESS) then
     errorcode = ierr
     call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
     print *,'Error after call to MPI_Irecv: ',errorstring
  endif
end do

call MPI_Waitall(numMessages,sendRequest,sendStatus,ierr)
call MPI_Waitall(numMessages,recvRequest,recvStatus,ierr)

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
real*8, allocatable :: tmpDomain(:,:,:)
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
call setupComm(numMessages)

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

nx = get_num_cells(1)
ny = get_num_cells(2)
nz = get_num_cells(3)

dx = 6.28/nx
dy = 6.28/ny
dz = 6.28/nz
angleAlpha = 3.14/2.0
angleBeta = 3.14/2.0
nu1 = 1.0
nu2 = 10.0
nu3 = 1.0

call setOperator(nu1,nu2,nu3,dx,dy,dz,angleAlpha,angleBeta)

! The pack is incorrect

numIters = 100

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
do iIter = 1, numIters
  if (rank == 0) print *, "running iter ", iIter
  ! Do two relaxations per iteration (to keep the data in domain)

  call packData(domain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, packBuffer, msgDirIndex)

  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)


  call unpackData(domain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, unpackBuffer, msgDirIndex)

  call relaxation(domain, tmpDomain, gridIndicesStart, gridIndicesEnd)

  call packData(tmpDomain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, packBuffer, msgDirIndex)

  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)


  call unpackData(tmpDomain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, unpackBuffer, msgDirIndex)

  call relaxation(tmpDomain, domain, gridIndicesStart, gridIndicesEnd)

enddo


call MPI_Finalize(ierr)

end program nodalCommFtn

