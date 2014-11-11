
#define NP 4

#define NLEV 26

module dimensions
  integer, parameter :: np = NP
  integer, parameter :: nlev = NLEV
end module dimensions

module element_data
use dimensions
implicit none
type :: element_t
  real*8 :: vIn(np,np,nlev)
  real*8 :: vOut(np,np,nlev)
end type

end module element_data

module pack
implicit none

public :: packSides
public :: unpackSides

public :: packElements
public :: unpackElements
contains

subroutine packSides(elements,thrElemStart,thrElemEnd,buffer)
  use element_data
  use dimensions
  implicit none
  type(element_t), intent(in) :: elements(:)
  integer, intent(in) :: thrElemStart, thrElemEnd
  real*8, intent(out) :: buffer(:)

  integer :: iElem, k, j, i
  integer :: elemPackSize, elemOffset, dirOffset, offset

  elemPackSize = np**2*nlev
  do iElem = thrElemStart,thrElemEnd
    elemOffset = (iElem - 1) * elemPackSize
    do k = 1,nlev

      ! South
      dirOffset = 0
      offset = elemOffset + dirOffset + (k-1)*np
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(i,1,k)
      enddo

      ! West
      dirOffset = np*nlev
      offset = elemOffset + dirOffset + (k-1)*np
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(i,1,k)
      enddo

      ! East
      dirOffset = 2*np*nlev
      offset = elemOffset + dirOffset + (k-1)*np
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(np,i,k)
      enddo

      ! North
      dirOffset = 3*np*nlev
      offset = elemOffset + dirOffset + (k-1)*np
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(i,np,k)
      enddo

    enddo
  enddo
end subroutine packSides

subroutine unpackSides(elements,thrElemStart,thrElemEnd,buffer)
  use element_data
  use dimensions
  implicit none
  type(element_t), intent(inout) :: elements(:)
  integer, intent(in) :: thrElemStart, thrElemEnd
  real*8, intent(in) :: buffer(:)

  integer :: iElem, k, j, i
  integer :: elemPackSize, elemOffset, dirOffset, kOffset, offset

  elemPackSize = np**2*nlev

  do iElem = thrElemStart,thrElemEnd
    elemOffset = (iElem - 1) * elemPackSize
    do k = 1,nlev
      kOffset = 4*np*(k-1)

      ! South
      dirOffset = 0
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(i,1,k) = buffer(offset+i)
      enddo

      ! West
      dirOffset = np*nlev
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(i,1,k) = buffer(offset+i)
      enddo

      ! East
      dirOffset = 2*np*nlev
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(np,i,k) = buffer(offset+i)
      enddo

      ! North
      dirOffset = 3*np*nlev
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(i,np,k) = buffer(offset+i)
      enddo

    enddo
  enddo
end subroutine unpackSides


subroutine packElements(elements,thrElemStart,thrElemEnd,buffer)
  use element_data
  use dimensions
  implicit none
  type(element_t), intent(in) :: elements(:)
  integer, intent(in) :: thrElemStart, thrElemEnd
  real*8, intent(out) :: buffer(:)

  integer :: iElem, k, j, i
  integer :: elemPackSize, elemOffset, dirOffset, kOffset, offset

  elemPackSize = 4*np*nlev
  do iElem = thrElemStart,thrElemEnd
    elemOffset = (iElem - 1) * elemPackSize
    do k = 1,nlev
      kOffset = 4*np*(k-1)

      ! South
      dirOffset = 0
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(i,1,k)
      enddo

      ! West
      dirOffset = np
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(1,i,k)
      enddo

      ! East
      dirOffset = 2*np
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(np,i,k)
      enddo

      ! North
      dirOffset = 3*np
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        buffer(offset+i) = elements(iElem)%vIn(i,np,k)
      enddo

    enddo
  enddo
end subroutine packElements

subroutine unpackElements(elements,thrElemStart,thrElemEnd,buffer)
  use element_data
  use dimensions
  implicit none
  type(element_t), intent(inout) :: elements(:)
  integer, intent(in) :: thrElemStart, thrElemEnd
  real*8, intent(in) :: buffer(:)

  integer :: iElem, k, j, i
  integer :: elemPackSize, elemOffset, dirOffset, offset, kOffset

  elemPackSize = 4*np*nlev

  do iElem = thrElemStart,thrElemEnd
    elemOffset = (iElem - 1) * elemPackSize
    do k = 1,nlev
      kOffset = 4*np*(k-1)

      ! South
      dirOffset = 0
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(i,1,k) = buffer(offset+i)
      enddo

      ! West
      dirOffset = np
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(1,i,k) = buffer(offset+i)
      enddo

      ! East
      dirOffset = 2*np
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(np,i,k) = buffer(offset+i)
      enddo

      ! North
      dirOffset = 3*np
      offset = elemOffset + kOffset + dirOffset
      do i = 1, np
        elements(iElem)%vOut(i,np,k) = buffer(offset+i)
      enddo

    enddo
  enddo
end subroutine unpackElements


end module pack





program align
use element_data
use dimensions
use pack
use omp_lib
implicit none

type(element_t), pointer :: elements(:)
real*8, pointer :: buffer(:)

integer :: numElems
integer :: numThreads
integer :: iElem, k, j, i
integer :: thrElemStart, thrElemEnd
integer :: bufSize
integer :: ithr, elemPackSize, elemOffset, dirOffset, offset
integer :: errorCounter

call read_namelist(numElems, numThreads)

! Get the number of threads

! Determine the number of elements
allocate(elements(numElems))

! Get the number of threads

! Create a linear buffer
!   Calculate the size of the buffer
bufSize = 16*numElems*nlev
allocate(buffer(bufSize))
buffer = -1
#ifdef _OPENMP
!$OMP PARALLEL NUM_THREADS(numThreads), DEFAULT(SHARED), PRIVATE(ithr,thrElemStart,thrElemEnd)
ithr=omp_get_thread_num()
thrElemStart = (numElems/numThreads)*ithr + 1
if (ithr == numThreads - 1) then
  thrElemEnd = numElems
else
  thrElemEnd = (numElems/numThreads)*(ithr+1)
endif
#else
ithr=1
thrElemStart = 1
thrElemEnd = numElems
#endif
print *, "ithr ", ithr, " owns elements ", thrElemStart, " - ", thrElemEnd

do iElem = thrElemStart, thrElemEnd
  elements(iElem)%vIn = -2.0
  elements(iElem)%vOut = -3.0
enddo 

do iElem = thrElemStart, thrElemEnd
  do k = 1, nlev
    do j = 1, np
      do i = 1, np
        elements(iElem)%vIn(i,j,k) = (iElem-1)*np**2*nlev + (k-1)*np**2 + (j-1)*np + i
      enddo
    enddo
  enddo
enddo

! have threads index into that buffer

! Pack
call packElements(elements,thrElemStart,thrElemEnd,buffer)
!$OMP BARRIER
! Unpack
call unpackElements(elements,thrElemStart,thrElemEnd,buffer)
!$OMP BARRIER
!$OMP END PARALLEL

! Check the data is the same
errorCounter = 0
do iElem = 1, numElems
  do k = 1,nlev
    do j = 1, np
      do i = 1, np
        ! Only check the edges
        if (mod(i,np - 1) .ne. 1 .and. mod(j,np - 1) .ne. 1) then
          cycle
        endif
        if (elements(iElem)%vIn(i,j,k) .ne. elements(iElem)%vOut(i,j,k)) then
          write(*,*) "Difference of ", elements(iElem)%vIn(i,j,k), elements(iElem)%vOut(i,j,k), &
                     " at (iElem,i,j,k) ", iElem, i, j, k
          errorCounter = errorCounter + 1
        endif
      enddo
    enddo
  enddo
enddo

if (errorCounter > 0) then
  print *, "There were errors"
else
  print *, "There were no errors"
endif
end program align

subroutine read_namelist(numElems, numThreads)
!use, intrinsic :: iso_fortran_env, only : input_unit => stdin
implicit none
integer, intent(out) :: numElems, numThreads

namelist /align_data/ numElems, numThreads

! Set default values
numElems = 4
numThreads = 1

read(*,nml=align_data)

write(*,*) "Read namelist data: "
write(*,nml=align_data)

end subroutine read_namelist


