module Mashm_enum_mod
  use, intrinsic :: iso_c_binding
end module


module Mashm_mod
#include "Mashmf.h"
implicit none

  enum, bind(c)
    enumerator :: mashm_send, mashm_receive
  end enum

  type MashmBufferPointer
    real*8, pointer :: p(:)
    type(c_ptr) :: c_ptr
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
    integer, value :: commIndex
    integer, value :: pairRank
    integer, value :: msgSize
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
    integer, value :: commIndex
    integer :: msgSize
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

  subroutine MashmGetBufferPointerC(in_mashm, i, sendReceive, bufferPointer)
    use, intrinsic :: iso_c_binding
    use Mashm_enum_mod
    Mashm, value :: in_mashm
    integer :: i
    integer(kind(mashm_send)) :: sendReceive
    type(c_ptr), value :: bufferPointer

  end subroutine

end interface

contains

  subroutine MashmGetBufferPointer(in_mashm, i, sendReceive, ftnBufferPointer)
    use, intrinsic :: iso_c_binding
    use Mashm_enum_mod
    Mashm, value :: in_mashm
    integer, intent(in) :: i
    integer(kind(mashm_send)), intent(in) :: sendReceive
    type(MashmBufferPointer), intent(out) :: ftnBufferPointer
    integer :: msgSize
    integer  :: iBaseZero
    iBaseZero = i - 1
    call MashmGetBufferPointerC(in_mashm, iBaseZero, sendReceive, &
             ftnBufferPointer%c_ptr)
    ! get size?
    msgSize = MashmGetCommSize(in_mashm, iBaseZero)

    !call Mashm
    call c_f_pointer(ftnBufferPointer%c_ptr,ftnBufferPointer%p,(/msgSize/))

  end subroutine

  subroutine MashmSetComm(in_mashm, commIndex, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: commIndex
    integer, value :: pairRank
    integer, value :: msgSize
    integer :: commIndexBaseZero
    commIndexBaseZero = commIndex - 1
    call MashmSetCommC(in_mashm, commIndexBaseZero, pairRank, msgSize)
  end subroutine 

end module Mashm_mod

