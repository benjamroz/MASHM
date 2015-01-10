#define DEBUG 0

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

! Faces
integer, parameter :: bottomVal = 4
integer, parameter :: topVal = 22
integer, parameter :: northVal = 16
integer, parameter :: southVal = 10
integer, parameter :: eastVal = 14
integer, parameter :: westVal = 12

! Edges
integer, parameter :: bsVal = 1
integer, parameter :: bnVal = 7
integer, parameter :: beVal = 5
integer, parameter :: bwVal = 3

integer, parameter :: tsVal = 19
integer, parameter :: tnVal = 25
integer, parameter :: teVal = 23
integer, parameter :: twVal = 21

integer, parameter :: swVal = 9
integer, parameter :: seVal = 11
integer, parameter :: nwVal = 15
integer, parameter :: neVal = 17

contains

function packDirInt(packDir3d) result (packInt)
  implicit none
  integer, intent(in) :: packDir3d(3)
  
  integer :: packInt

  packInt = (packDir3d(1)+1) + 3*(packDir3d(2)+1) + 9*(packDir3d(3)+1)

end function

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

subroutine packData(inData, gridIndicesStart, gridIndicesEnd,  &
                    msgOffsets, msgLengths, packBuffer, packDir)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
integer, intent(in) :: msgOffsets(:), msgLengths(:)
real*8, intent(out) :: packBuffer(:)
integer, intent(in) :: packDir(-1:1,-1:1,-1:1)
integer :: ix, iy, iz
integer :: packIndex, packOffset

! pack top
packIndex = packDir(0,0,1)
if (packIndex > 0) then
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy 
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
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
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
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
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
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
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
      packOffset = packOffset + 1
    enddo
  enddo
endif

! pack west
packIndex = packDir(-1,0,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
      packOffset = packOffset + 1
    enddo
  enddo
endif

! Pack east
packIndex = packDir(1,0,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
      packOffset = packOffset + 1
    enddo
  enddo
endif

! Pack Middle Edges
! South West
packIndex = packDir(-1,-1,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iy = gridIndicesStart(2)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! South East
packIndex = packDir(1,-1,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iy = gridIndicesStart(2)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! North West
packIndex = packDir(-1,1,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iy = gridIndicesEnd(2)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! North East
packIndex = packDir(1,1,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iy = gridIndicesEnd(2)
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif


! Pack bottom edges
! Bottom West
packIndex = packDir(-1,0,-1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! Bottom East
packIndex = packDir(1,0,-1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iz = gridIndicesStart(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
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
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
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
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif


! Pack top edges
! Top West
packIndex = packDir(-1,0,1)
if (packIndex > 0) then
  ix = gridIndicesStart(1)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! Top East
packIndex = packDir(1,0,1)
if (packIndex > 0) then
  ix = gridIndicesEnd(1)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! Top South
packIndex = packDir(0,-1,1)
if (packIndex > 0) then
  iy = gridIndicesStart(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif

! Top North
packIndex = packDir(0,1,1)
if (packIndex > 0) then
  iy = gridIndicesEnd(2)
  iz = gridIndicesEnd(3)
  packOffset = msgOffsets(packIndex) + 1
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
      packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
      packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
endif
end subroutine

subroutine unpackData(inData, gridIndicesStart, gridIndicesEnd,  &
                    msgOffsets, msgLengths, unpackBuffer, packDir)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                                gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                               gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
integer, intent(in) :: msgOffsets(:), msgLengths(:)
real*8, intent(in) :: unpackBuffer(:)
integer, intent(in) :: packDir(-1:1,-1:1,-1:1)
integer :: ix, iy, iz

integer :: packIndex, packOffset

! pack top
packIndex = packDir(0,0,1)
if (packIndex > 0) then
  iz = gridIndicesEnd(3) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    do ix = gridIndicesStart(1), gridIndicesEnd(1)
      !packBuffer(packOffset) = inData(ix,iy,iz)
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
      packOffset = packOffset + 1
    enddo
  enddo
endif

! Pack Middle Edges
! South West
packIndex = packDir(-1,-1,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iy = gridIndicesStart(2) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    !packBuffer(packOffset) = inData(ix,iy,iz)
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
    packOffset = packOffset + 1
  enddo
endif

! South East
packIndex = packDir(1,-1,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iy = gridIndicesStart(2) - 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    !packBuffer(packOffset) = inData(ix,iy,iz)
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
    packOffset = packOffset + 1
  enddo
endif

! North West
packIndex = packDir(-1,1,0)
if (packIndex > 0) then
  ix = gridIndicesStart(1) - 1
  iy = gridIndicesEnd(2) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    !packBuffer(packOffset) = inData(ix,iy,iz)
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
    packOffset = packOffset + 1
  enddo
endif

! North East
packIndex = packDir(1,1,0)
if (packIndex > 0) then
  ix = gridIndicesEnd(1) + 1
  iy = gridIndicesEnd(2) + 1
  packOffset = msgOffsets(packIndex) + 1
  do iz = gridIndicesStart(3), gridIndicesEnd(3)
    !packBuffer(packOffset) = inData(ix,iy,iz)
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
    packOffset = packOffset + 1
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
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
#if DEBUG
      if (unpackBuffer(packOffset) .ne. (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy) then
        print *, "Error in unpack, stopping"
        stop -1
      endif
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#else
      inData(ix,iy,iz) = unpackBuffer(packOffset)
#endif
    packOffset = packOffset + 1
  enddo
endif
end subroutine

subroutine packData2(inData, gridIndicesStart, gridIndicesEnd,  &
                     packBuffer, packIndex)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer, intent(in) :: packIndex
integer :: ix, iy, iz

select case (packIndex)
  case (topVal)
    call packDataTop(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (bottomVal)
    call packDataBottom(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (eastVal)
    call packDataEast(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (westVal)
    call packDataWest(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (northVal)
    call packDataNorth(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (southVal)
    call packDataSouth(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (tnVal)
    call packDataTN(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (tsVal)
    call packDataTS(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (teVal)
    call packDataTE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (twVal)
    call packDataTW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (bnVal)
    call packDataBN(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (bsVal)
    call packDataBS(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (beVal)
    call packDataBE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (bwVal)
    call packDataBW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (swVal)
    call packDataSW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (seVal)
    call packDataSE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (nwVal)
    call packDataNW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  case (neVal)
    call packDataNE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
  end select
end subroutine

subroutine packDataTop(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! pack top
iz = gridIndicesEnd(3)
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy 
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
enddo
end subroutine


subroutine packDataBottom(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! pack bottom
iz = gridIndicesStart(3)
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine packDataSouth(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! pack south
iy = gridIndicesStart(2)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine packDataNorth(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! pack north
iy = gridIndicesEnd(2)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine packDataWest(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Pack West
ix = gridIndicesStart(1)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine packDataEast(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Pack East
ix = gridIndicesEnd(1)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

! Pack Middle Edges
subroutine packDataSW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! South West
ix = gridIndicesStart(1)
iy = gridIndicesStart(2)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataSE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! South East
ix = gridIndicesEnd(1)
iy = gridIndicesStart(2)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataNW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! North West
ix = gridIndicesStart(1)
iy = gridIndicesEnd(2)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataNE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! North East
ix = gridIndicesEnd(1)
iy = gridIndicesEnd(2)
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine


! Pack bottom edges
subroutine packDataBW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Bottom West
ix = gridIndicesStart(1)
iz = gridIndicesStart(3)
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataBE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Bottom East
ix = gridIndicesEnd(1)
iz = gridIndicesStart(3)
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataBS(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Bottom south
iy = gridIndicesStart(2)
iz = gridIndicesStart(3)
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataBN(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Bottom north
iy = gridIndicesEnd(2)
iz = gridIndicesStart(3)
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

! Pack top edges
subroutine packDataTW(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Top West
ix = gridIndicesStart(1)
iz = gridIndicesEnd(3)
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataTE(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Top East
ix = gridIndicesEnd(1)
iz = gridIndicesEnd(3)
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataTS(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Top south
iy = gridIndicesStart(2)
iz = gridIndicesEnd(3)
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine packDataTN(inData, gridIndicesStart, gridIndicesEnd, packBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(in) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(out) :: packBuffer(:)
integer :: ix, iy, iz
integer :: packOffset
! Top north
iy = gridIndicesEnd(2)
iz = gridIndicesEnd(3)
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
#if DEBUG
    packBuffer(packOffset) = (ix - 1) + (iy - 1) * globNx + iz * globNx * globNy
#else
    packBuffer(packOffset) = inData(ix,iy,iz)
#endif
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackData2(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer, packIndex)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer, intent(in) :: packIndex
integer :: ix, iy, iz

select case (packIndex)
  case (topVal)
    call unpackDataTop(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (bottomVal)
    call unpackDataBottom(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (eastVal)
    call unpackDataEast(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (westVal)
    call unpackDataWest(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (northVal)
    call unpackDataNorth(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (southVal)
    call unpackDataSouth(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (tnVal)
    call unpackDataTN(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (tsVal)
    call unpackDataTS(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (teVal)
    call unpackDataTE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (twVal)
    call unpackDataTW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (bnVal)
    call unpackDataBN(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (bsVal)
    call unpackDataBS(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (beVal)
    call unpackDataBE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (bwVal)
    call unpackDataBW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (swVal)
    call unpackDataSW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (seVal)
    call unpackDataSE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (nwVal)
    call unpackDataNW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  case (neVal)
    call unpackDataNE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
  end select
end subroutine

subroutine unpackDataTop(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! pack top
iz = gridIndicesEnd(3) + 1
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
enddo
end subroutine


subroutine unpackDataBottom(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! pack bottom
iz = gridIndicesStart(3) - 1
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine unpackDataSouth(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! pack south
iy = gridIndicesStart(2) - 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine unpackDataNorth(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! pack north
iy = gridIndicesEnd(2) + 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine unpackDataWest(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Pack West
ix = gridIndicesStart(1) - 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

subroutine unpackDataEast(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Pack East
ix = gridIndicesEnd(1) + 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
  do iy = gridIndicesStart(2), gridIndicesEnd(2)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
    packOffset = packOffset + 1
  enddo
enddo
end subroutine

! Pack Middle Edges
subroutine unpackDataSW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! South West
ix = gridIndicesStart(1) - 1
iy = gridIndicesStart(2) - 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataSE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! South East
ix = gridIndicesEnd(1) + 1
iy = gridIndicesStart(2) - 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataNW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! North West
ix = gridIndicesStart(1) - 1
iy = gridIndicesEnd(2) + 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataNE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! North East
ix = gridIndicesEnd(1) + 1
iy = gridIndicesEnd(2) + 1
packOffset = 1
do iz = gridIndicesStart(3), gridIndicesEnd(3)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine


! Pack bottom edges
subroutine unpackDataBW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Bottom West
ix = gridIndicesStart(1) - 1
iz = gridIndicesStart(3) - 1
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataBE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Bottom East
ix = gridIndicesEnd(1) + 1
iz = gridIndicesStart(3) - 1
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataBS(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Bottom south
iy = gridIndicesStart(2) - 1
iz = gridIndicesStart(3) - 1
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataBN(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Bottom north
iy = gridIndicesEnd(2) + 1
iz = gridIndicesStart(3) - 1
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

! Pack top edges
subroutine unpackDataTW(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Top West
ix = gridIndicesStart(1) - 1
iz = gridIndicesEnd(3) + 1
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataTE(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Top East
ix = gridIndicesEnd(1) + 1
iz = gridIndicesEnd(3) + 1
packOffset = 1
do iy = gridIndicesStart(2), gridIndicesEnd(2)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataTS(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Top south
iy = gridIndicesStart(2) - 1
iz = gridIndicesEnd(3) + 1
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
end subroutine

subroutine unpackDataTN(inData, gridIndicesStart, gridIndicesEnd, unpackBuffer)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(inout) :: inData(gridIndicesStart(1)-1:gridIndicesEnd(1)+1, &
                             gridIndicesStart(2)-1:gridIndicesEnd(2)+1, &
                             gridIndicesStart(3)-1:gridIndicesEnd(3)+1)
real*8, intent(in) :: unpackBuffer(:)
integer :: ix, iy, iz
integer :: packOffset

! Top north
iy = gridIndicesEnd(2) + 1
iz = gridIndicesEnd(3) + 1
packOffset = 1
do ix = gridIndicesStart(1), gridIndicesEnd(1)
    inData(ix,iy,iz) = unpackBuffer(packOffset)
  packOffset = packOffset + 1
enddo
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
Mashm :: myMashm
MashmCommType :: commMethod
type(MashmBufferPointer), allocatable :: mashmSendBufferPtrs(:)
type(MashmBufferPointer), allocatable :: mashmRecvBufferPtrs(:)

integer, allocatable :: msgDirIndex2(:)
integer :: msgIndex
integer :: packInt3(3)

real :: rate 
integer :: cr, cm, clockStart, clockEnd
real :: time1, time2, time3

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

! Initialize system_clock for timers
call system_clock(count_rate=cr)
call system_clock(count_max=cm)
rate = REAL(cr)

call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)

if (rank == 0) print *, "Initial difference", residualL2, residualMax

numIters = 100

! Stride two 
call system_clock(clockStart)
do iIter = 1, numIters, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, packBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, unpackBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(tmpDomain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
  do i = 1, numMessages
    call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, packBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)


  do i = 1, numMessages
    call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, unpackBuffer(msgOffsets(i)+1:), msgDirIndex2(i))
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax

enddo
call system_clock(clockEnd)
time1 = (clockEnd - clockStart)/rate

call MashmInit(myMashm, MPI_COMM_WORLD)

! Print nodal comm info
call MashmPrintInfo(myMashm)

call MashmSetNumComms(myMashm, numMessages)

! Add communications calculated above
do i = 1, numMessages
  ! Fortran to C indexing 
  call MashmSetComm(myMashm, i, neighborRanks(i), msgSizes(i))
enddo

! Choose the communitcation method
!commMethod = MASHM_COMM_INTRA_SHARED
!commMethod = MASHM_COMM_INTRA_MSG
!commMethod = MASHM_COMM_STANDARD
commMethod = MASHM_COMM_MIN_AGG

call MashmSetCommMethod(myMashm, commMethod)

! Perform precalculation
call MashmCommFinish(myMashm)

!call MashmPrintCommCollection(myMashm)

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

! Stride two 
call system_clock(clockStart)
do iIter = 1, numIters, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  !call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)
  call MashmInterNodeCommBegin(myMashm)
  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)
  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(tmpDomain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
  do i = 1, numMessages
    call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  !call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)

  call MashmInterNodeCommBegin(myMashm)
  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)
  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax

enddo
call system_clock(clockEnd)
time2 = (clockEnd - clockStart)/rate


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

! Stride two 
call system_clock(clockStart)
do iIter = 1, numIters, 2

  do i = 1, numMessages
    call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  do i = 1, numMessages
    if (.not. MashmIsMsgOnNode(myMashm, i)) then
      call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  !call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)

  call MashmInterNodeCommBegin(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgOnNode(myMashm, i)) then
      call packData2(domain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgOnNode(myMashm, i)) then
      call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (.not. MashmIsMsgOnNode(myMashm, i)) then
      call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  do i = 1, numMessages
    call unpackData2(domain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
  enddo

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(tmpDomain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
  do i = 1, numMessages
    if (.not. MashmIsMsgOnNode(myMashm, i)) then
      call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommBegin(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgOnNode(myMashm, i)) then
      call packData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmSendBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmIntraNodeCommBegin(myMashm)

  call MashmIntraNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (MashmIsMsgOnNode(myMashm, i)) then
      call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call MashmInterNodeCommEnd(myMashm)

  do i = 1, numMessages
    if (.not. MashmIsMsgOnNode(myMashm, i)) then
      call unpackData2(tmpDomain, gridIndicesStart, gridIndicesEnd, mashmRecvBufferPtrs(i)%p, msgDirIndex2(i))
    endif
  enddo

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax

enddo
call system_clock(clockEnd)
time3 = (clockEnd - clockStart)/rate



! Restore (nullify) the Mashm access pointers
do i = 1, numMessages
  call MashmRetireBufferPointer(myMashm, mashmSendBufferPtrs(i))
  call MashmRetireBufferPointer(myMashm, mashmRecvBufferPtrs(i))
enddo

! Destroy the Mashm object
call MashmDestroy(myMashm)

if (rank == 0) then
  print *, "Number of dofs per rank = ", totalNumCells/numProcs
  print *, "system_clock rate ",rate
  print *, "Time of original method", time1
  print *, "Time of mashm method", time2
  print *, "Time of mashm method with overlap", time3
endif

deallocate(packBuffer)
deallocate(unpackBuffer)

deallocate(domain)
deallocate(tmpDomain)
deallocate(solution)
deallocate(rhs)

call MPI_Finalize(ierr)

end program nodalCommFtn

