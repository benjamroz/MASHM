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

factor = -(dx**2 * dy**2 * dz**2)/(2.0*dy**2*dz**2 + 2.0*dx**2*dz**2 + 2.0*dx**2*dy**2)

derivXX = derivXX*factor
derivYY = derivYY*factor
derivZZ = derivZZ*factor
derivXY = derivXY*factor
derivXZ = derivXZ*factor
derivYZ = derivYZ*factor

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

!coefXX = nu1
!coefYY = nu2
!coefZZ = nu3
!coefXY = 0.0
!coefXZ = 0.0
!coefYZ = 0.0

coefs = coefXX * derivXX &
      + coefYY * derivYY &
      + coefZZ * derivZZ &
      + coefXY * derivXY &
      + coefXZ * derivXZ &
      + coefYZ * derivYZ

print *, "coefs(0,0,0) = ", coefs(0,0,0)
!coefs(0,0,0) = 0.0
coefs = -coefs
coefs(0,0,0) = 0.0
print *, "coefs(1,0,0) = ", coefs(1,0,0)
print *, "coefs(-1,0,0) = ", coefs(-1,0,0)
!coefs(0,0,0) = -1.0
rhsFactor = -factor
print *, "rhsFactor = ", rhsFactor
! \nabla \cdot : [-1, 1
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
#if 0
      !outData(ix,iy,iz) = (1.0-relaxOmega)*inData(ix,iy,iz) + (relaxOmega/6.0)*(globDx**2)*rhs(ix,iy,iz)
      !outData(ix,iy,iz) = outData(ix,iy,iz) + relaxOmega/6*( &
      !                    inData(ix-1,iy,iz) + inData(ix+1,iy,iz) + &
      !                    inData(ix,iy-1,iz) + inData(ix,iy+1,iz) + &
      !                    inData(ix,iy,iz-1) + inData(ix+1,iy,iz+1))
      outData(ix,iy,iz) = (1.0/6.0)*(globDx*globDx)*rhs(ix,iy,iz) + 1.0/6.0*( &
                          inData(ix-1,iy,iz) + inData(ix+1,iy,iz) + &
                          inData(ix,iy-1,iz) + inData(ix,iy+1,iz) + &
                          inData(ix,iy,iz-1) + inData(ix,iy,iz+1))
#else
      outData(ix,iy,iz) = (1.0-relaxOmega)*inData(ix,iy,iz) + relaxOmega*rhsFactor*rhs(ix,iy,iz)
      do k = -1, 1
        do j = -1, 1
          do i = -1, 1
            outData(ix,iy,iz) = outData(ix,iy,iz) + relaxOmega*coefs(i,j,k)*inData(ix+i,iy+j,iz+k)
          enddo
        enddo
      enddo
#endif
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

! pack east
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

! pack west
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
#if 0
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
! Bottom east
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

! Bottom west
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
! Bottom east
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

! Bottom west
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

! Bottom south
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

! Bottom north
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
#endif
#if 0
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
#endif
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
#if 0
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
#endif
#if 0
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
#endif
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

subroutine setSolution(solution, gridIndicesStart, gridIndicesEnd, dx, dy, dz)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(out) :: solution(gridIndicesStart(1):gridIndicesEnd(1), &
                                gridIndicesStart(2):gridIndicesEnd(2), &
                                gridIndicesStart(3):gridIndicesEnd(3))
real*8, intent(in) :: dx, dy, dz

integer :: i, j, k
real*8 :: x, y, z

do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      !solution(i,j,k) = dsin(i*dx)
      !solution(i,j,k) = 0.0
      x = i*dx
      y = j*dy
      z = k*dz
      solution(i,j,k) = sin(x)*sin(y)*sin(z)
      !solution(i,j,k) = 0.0
      !solution(i,j,k) = (x**2 - x**4)*(y**2 - y**4)*(z**2 - z**4)
      !solution(i,j,k) = x*(1-x) + y*(1-y) + z*(1-z)
      !solution(i,j,k) = x*(1-x)
    enddo
  enddo
enddo
end subroutine 

subroutine setRhs(rhs, gridIndicesStart, gridIndicesEnd, dx, dy, dz)
implicit none
integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)
real*8, intent(out) :: rhs(gridIndicesStart(1):gridIndicesEnd(1), &
                           gridIndicesStart(2):gridIndicesEnd(2), &
                           gridIndicesStart(3):gridIndicesEnd(3))
real*8, intent(in) :: dx, dy, dz

integer :: i, j, k
real*8 :: x, y, z

do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      !rhs(i,j,k) = -(1.0+1.0+1.0) * dsin(i*dx) * dsin(j*dy) * dsin(k*dz)
      !rhs(i,j,k) = dsin(i*dx)
      !rhs(i,j,k) = 0.0
      x = i*dx
      y = j*dy
      z = k*dz
      !rhs(i,j,k) = 3.0*dsin(x)*dsin(y)*dsin(z)
      rhs(i,j,k) = 3.0*sin(x)*sin(y)*sin(z)
      !rhs(i,j,k) = 0.0
      !rhs(i,j,k) = -(2 - 12*x**2)*(y**2 - y**4)*(z**2 - z**4) &
      !           - (2 - 12*y**2)*(x**2 - x**4)*(z**2 - z**4) &
      !           - (2 - 12*z**2)*(x**2 - x**4)*(y**2 - y**4)
      !rhs(i,j,k) = 6
      !rhs(i,j,k) = 2
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
integer :: iSave, jSave, kSave

!integer :: counter
!counter = 1
localL2Norm = 0.0
localMaxNorm = 0.0
do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      localL2Norm = localL2Norm + (numSoln(i,j,k) - solution(i,j,k))**2
      !localMaxNorm = maxval( (/ localMaxNorm, dabs( numSoln(i,j,k) - solution(i,j,k) ) /) )
      newMaxNorm = maxval( (/ localMaxNorm, dabs( numSoln(i,j,k) - solution(i,j,k) ) /) )
      if (newMaxNorm > localMaxNorm) then
        iSave = i
        jSave = j
        kSave = k
      endif
      localMaxNorm = newMaxNorm
      !localL2Norm = localL2Norm + (solution(i,j,k))**2
      !print *, "localL2Norm = ", localL2Norm, ", ", counter
      !counter = counter + 1
    enddo
  enddo
enddo

! MPI_Allreduce 
call MPI_Allreduce(localL2Norm, globalL2Norm, 1, MPI_REAL8, MPI_SUM, &
                   MPI_COMM_WORLD, ierr)
call MPI_Allreduce(localMaxNorm, globalMaxNorm, 1, MPI_REAL8, MPI_MAX, &
                   MPI_COMM_WORLD, ierr)

globalL2Norm = dsqrt(globalL2Norm/totalNumCells)
#if 0
print *, "iSave,jSave,kSave = ", iSave, jSave, kSave

print *, "maxvals of edges:"
if (gridIndicesStart(1) == 1) then
  print *, maxval(numSoln(gridIndicesStart(1) - 1,:,:))
endif
if (gridIndicesEnd(1) == globNx) then
  print *, maxval(numSoln(gridIndicesEnd(1) + 1,:,:))
endif
if (gridIndicesStart(2) == 1) then
  print *, maxval(numSoln(gridIndicesStart(2) - 1,:,:))
endif
if (gridIndicesEnd(2) == globNy) then
  print *, maxval(numSoln(gridIndicesEnd(2) + 1,:,:))
endif
if (gridIndicesStart(3) == 1) then
  print *, maxval(numSoln(gridIndicesStart(3) - 1,:,:))
endif
if (gridIndicesEnd(3) == globNz) then
  print *, maxval(numSoln(gridIndicesEnd(3) + 1,:,:))
endif
#endif
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

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)

totalMessageSize=msgOffsets(numMessages+1)
if (rank == 2) then
  print *, "Rank ", rank, " message sizes offsets"
  do i = 1, numMessages
    print *, "i ", i, " rank ", neighborRanks(i), " msgSize ", msgSizes(i), " begin,end", &
             msgOffsets(i:i+1)
  enddo
  print *, "msgDirIndex:"
  counter = 1
  do iz = -1, 1
    do iy = -1, 1
      do ix = -1, 1
        print *, "  ", counter, "  ", msgDirIndex(ix,iy,iz)
        counter = counter + 1
      enddo
    enddo
  enddo
endif
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
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
if (rank == 0) print *, "nx,ny,nz = ", nx, ny, nz
if (rank == 0) print *, "dx,dy,dz = ", dx, dy, dz
if (rank == 0) print *, "lx, ly, lz = ", get_grid_length(1), get_grid_length(2), get_grid_length(3)
if (rank == 0) print *, "gx begin, end:", gridIndicesStart(1), gridIndicesEnd(1)
if (rank == 0) print *, "gy begin, end:", gridIndicesStart(2), gridIndicesEnd(2)
if (rank == 0) print *, "gz begin, end:", gridIndicesStart(3), gridIndicesEnd(3)
if (rank == 1) print *, "gx begin, end:", gridIndicesStart(1), gridIndicesEnd(1)
if (rank == 1) print *, "gy begin, end:", gridIndicesStart(2), gridIndicesEnd(2)
if (rank == 1) print *, "gz begin, end:", gridIndicesStart(3), gridIndicesEnd(3)

globDx = dx
globDy = dy
globDz = dz

!angleAlpha = 3.14/2.0
!angleBeta = 3.14/2.0
angleAlpha = 0.0
angleBeta = 0.0
nu1 = 1.0
nu2 = 1.0
nu3 = 1.0

! Initialize to zero

call setRhs(rhs,gridIndicesStart,gridIndicesEnd,dx,dy,dz)
call setSolution(solution,gridIndicesStart,gridIndicesEnd,dx,dy,dz)


domain = 0.0d0
tmpDomain = 0.0d0
do k = gridIndicesStart(3), gridIndicesEnd(3)
  do j = gridIndicesStart(2), gridIndicesEnd(2)
    do i = gridIndicesStart(1), gridIndicesEnd(1)
      !domain(i,j,k) = dsin(i*dx)*dsin(j*dy)*dsin(k*dz)+2*dsin(6*i*dx)*dsin(6*j*dy)*dsin(6*k*dz)
      domain(i,j,k) = dsin(10*i*dx)*dsin(10*j*dy)*dsin(10*k*dz)
      !domain(i,j,k) = dsin(i*dx)*dsin(j*dy)*dsin(k*dz)
      !domain(i,j,k) = 1.0
      !domain(i,j,k) = 0.33333*(dsin(i*dx) + dsin(6*i*dx) + dsin(32*i*dx))
      !domain(i,j,k) = dsin(3*i*dx)*dsin(3*
      !domain(i,j,k) = dsin(10*i*dz)
      !domain(i,j,k) = dsin(10*i*dx*PI)
      !domain(i,j,k) = 0.0
    enddo
  enddo
enddo

!call setSolution(domain(gridIndicesStart(1):gridIndicesEnd(1),  &
!                        gridIndicesStart(2):gridIndicesEnd(2),  &
!                        gridIndicesStart(3):gridIndicesEnd(3)), & 
!                        gridIndicesStart,gridIndicesEnd,dx,dy,dz)


call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                totalNumCells)
if (rank == 0) print *, "Initial difference", residualL2, residualMax

call setOperator(nu1,nu2,nu3,dx,dy,dz,angleAlpha,angleBeta)

if (rank == 0) print *, "coefs = ", coefs
! The pack is incorrect

numIters = 1000
!numIters = 1

call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call flush(6)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)

packBuffer = 666
unpackBuffer = 666

do iIter = 1, numIters, 2

  call packData(domain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, packBuffer, msgDirIndex)

  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)

  call unpackData(domain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, unpackBuffer, msgDirIndex)

  call relaxation(domain, tmpDomain, rhs, gridIndicesStart, gridIndicesEnd)

  domain = tmpDomain
  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter, " residual ", residualL2, &
                          residualMax
  cycle
  call packData(tmpDomain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, packBuffer, msgDirIndex)

  call communication(packBuffer, unpackBuffer, numMessages, neighborRanks, msgSizes, msgOffsets)


  call unpackData(tmpDomain, gridIndicesStart, gridIndicesEnd,  &
                msgOffsets, msgSizes, unpackBuffer, msgDirIndex)

  call relaxation(tmpDomain, domain, rhs, gridIndicesStart, gridIndicesEnd)

  call calcL2Norm(domain, solution, gridIndicesStart, gridIndicesEnd, residualL2, residualMax, &
                  totalNumCells)

  if (rank == 0) print *, "running iter ", iIter + 1, " residual ", residualL2, &
                          residualMax

enddo

call MPI_Finalize(ierr)

end program nodalCommFtn

