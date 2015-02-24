module poisson3dPack

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

end module
