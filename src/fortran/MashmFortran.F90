module Mashm_enum_mod
  use, intrinsic :: iso_c_binding
end module


module Mashm_mod
    use, intrinsic :: iso_c_binding
#include "Mashmf.h"
implicit none

  !enum, bind(c) :: MashmSendReceive
  !  enumerator :: mashm_send, mashm_receive
  !end enum
  integer(c_int), parameter :: MASHM_SEND = 0, MASHM_RECEIVE = 1

  type MashmBufferPointer
    !real*8, pointer :: p(:)
    real(c_double), pointer, public :: p(:)
    type(c_ptr), public :: cPtr
  end type MashmBufferPointer

  type MashmPointer1d
    real*8, allocatable :: p(:)
  end type MashmPointer1d

interface

  subroutine MashmInit(in_mashm, in_comm)
    use, intrinsic :: iso_c_binding
    Mashm :: in_mashm
    integer, value :: in_comm
  end subroutine

  subroutine MashmDestroy(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm :: in_mashm
  end subroutine

  subroutine MashmPrintInfo(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine MashmPrintCommCollection(in_mashm);
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine MashmSetNumComms(in_mashm, numComms)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: numComms
  end subroutine 

  subroutine MashmAddSymComm(in_mashm, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: pairRank
    integer, value :: msgSize
  end subroutine 

  subroutine MashmSetCommC(in_mashm, commIndex, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer(c_int), value :: commIndex
    integer(c_int), value :: pairRank
    integer(c_int), value :: msgSize
  end subroutine 

  function MashmGetCommRank(in_mashm, commIndex) result(msgRank)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: commIndex
    integer :: msgRank
  end function

  function MashmGetCommSize(in_mashm, commIndex) result(msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer(c_int), value :: commIndex
    integer(c_int) :: msgSize
  end function

  subroutine MashmCommFinish(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine MashmInterNodeCommBegin(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine MashmInterNodeCommEnd(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine MashmIntraNodeCommBegin(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine MashmIntraNodeCommEnd(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  !subroutine MashmGetBufferPointerC(in_mashm, i, sendReceive, bufferPointer)
  !  use, intrinsic :: iso_c_binding
  !  use Mashm_enum_mod
  !  Mashm, value, intent(in) :: in_mashm
  !  integer(c_int), value, intent(in) :: i
  !  integer(c_int), value, intent(in) :: sendReceive
  !  type(c_ptr), intent(out) :: bufferPointer
  !end subroutine

  !function MashmGetBufferPointerC(in_mashm, i, sendReceive) result(bufferPointer)
  !  use, intrinsic :: iso_c_binding
  !  use Mashm_enum_mod
  !  Mashm, value, intent(in) :: in_mashm
  !  integer(c_int), value, intent(in) :: i
  !  integer(c_int), value, intent(in) :: sendReceive
  !  type(c_ptr) :: bufferPointer
  !end function

  subroutine MashmGetBufferPointer2C(in_mashm, i, sendReceive, bufferPointer)
    use, intrinsic :: iso_c_binding
    use Mashm_enum_mod
    Mashm, value, intent(in) :: in_mashm
    integer(c_int), value, intent(in) :: i
    integer(c_int), value, intent(in) :: sendReceive
    type(c_ptr), intent(out) :: bufferPointer
  end subroutine

  function MashmIsMsgOnNodeC(in_mashm, commIndex) result(isOnNode)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: commIndex
    integer :: isOnNode
   
  end function

end interface

contains

  subroutine MashmGetBufferPointer(in_mashm, i, sendReceive, ftnBufferPointer, inMsgSize)
    use, intrinsic :: iso_c_binding
    use Mashm_enum_mod
    implicit none
    Mashm, value :: in_mashm
    integer, intent(in) :: i
    integer(c_int), intent(in) :: sendReceive
    type(MashmBufferPointer), intent(inout) :: ftnBufferPointer
    integer(c_int) :: msgSize
    integer(c_int) :: iBaseZero
    integer :: inMsgSize
    iBaseZero = i - 1

    ! TODO: move to MashmGetBufferPointerC function call - not working for some
    !       reason
    call MashmGetBufferPointer2C(in_mashm, iBaseZero, sendReceive, ftnBufferPointer%cPtr)
    !ftnBufferPointer%cPtr = MashmGetBufferPointerC(in_mashm, iBaseZero, sendReceive)

    msgSize = MashmGetCommSize(in_mashm, iBaseZero)
    !print *, "MashmGetBufferPointer size, inMsgSize = ", msgSize, inMsgSize

    call c_f_pointer(cptr=ftnBufferPointer%cPtr,fptr=ftnBufferPointer%p,shape = (/ msgSize /))

  end subroutine

  subroutine MashmSetComm(in_mashm, commIndex, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: commIndex
    integer(c_int), value :: pairRank
    integer(c_int), value :: msgSize
    integer(c_int) :: commIndexBaseZero
    commIndexBaseZero = commIndex - 1
    call MashmSetCommC(in_mashm, commIndexBaseZero, pairRank, msgSize)
  end subroutine 

  function MashmIsMsgOnNode(in_mashm, commIndex) result(isOnNode)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: commIndex
    logical :: isOnNode
   
    isOnNode = ( MashmIsMsgOnNodeC(in_mashm, commIndex - 1) .eq. 1 )

  end function
end module Mashm_mod

