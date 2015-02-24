module Mashm_enum_mod
  use, intrinsic :: iso_c_binding
end module

module Mashm_type
  use, intrinsic :: iso_c_binding
  implicit none
  type, bind(c) :: Mashm
    type(c_ptr) :: p
  end type
end module 

module Mashm_mod
    use, intrinsic :: iso_c_binding
#include "Mashmf.h"
implicit none
 
!  type, bind(c) :: Mashm
!    type(c_ptr) :: p
!  end type

  ! Set enump for MashmSendReceive
  enum, bind(c) 
    enumerator  :: MASHM_SEND = 0, MASHM_RECEIVE = 1
  end enum

  ! Set enump for MashmCommType
  enum, bind(c)
    enumerator :: MASHM_COMM_STANDARD, &
                  MASHM_COMM_INTRA_MSG, &
                  MASHM_COMM_INTRA_SHARED, &
                  MASHM_COMM_MIN_AGG
  end enum
#if 0
  ENUM, BIND(C) :: MashmCommType
    enumerator MASHM_COMM_STANDARD = 0, &
    enumerator MASHM_COMM_INTRA_MSG = 1, &
    enumerator MASHM_COMM_INTRA_SHARED = 2, &
    enumerator MASHM_COMM_MIN_AGG = 3
  END ENUM
#endif

  type MashmBufferPointer
    real(c_double), pointer, public :: p(:)
    type(c_ptr), public :: cPtr
  end type MashmBufferPointer

interface

  subroutine MashmInit(in_mashm, in_comm) &
    bind(c, name='MashmInit')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    !Mashm :: in_mashm
    type(Mashm) :: in_mashm
    integer(c_int), value :: in_comm
  end subroutine

  subroutine MashmDestroy(in_mashm) &
    bind(c, name='MashmDestroy')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm) :: in_mashm
  end subroutine

  subroutine MashmPrintInfo(in_mashm) &
    bind(c,name='MashmPrintInfo')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmPrintCommCollection(in_mashm) &
    bind(c, name='MashmPrintCommCollection')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmPrintMessageStats(in_mashm) &
    bind(c,name='MashmPrintMessageStats')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmSetNumComms(in_mashm, numComms) &
    bind(c,name='MashmSetNumComms')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, value :: numComms
  end subroutine 

  subroutine MashmSetCommMethod(in_mashm, commMethod) &
    bind(c, name='MashmSetCommMethod')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer(c_int), value :: commMethod
  end subroutine


  subroutine MashmAddSymComm(in_mashm, pairRank, msgSize) &
    bind(c, name='MashmAddSymCom')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, value :: pairRank
    integer, value :: msgSize
  end subroutine 

  subroutine MashmSetCommC(in_mashm, commIndex, pairRank, msgSize) &
    bind(c, name='MashmSetComm')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer(c_int), value :: commIndex
    integer(c_int), value :: pairRank
    integer(c_int), value :: msgSize
  end subroutine 

  function MashmGetCommRank(in_mashm, commIndex) result(msgRank) &
    bind(c, name='MashmGetCommRank')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, value :: commIndex
    integer :: msgRank
  end function

  function MashmGetCommSize(in_mashm, commIndex) result(msgSize) &
    bind(c, name='MashmGetCommSize')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    !type(Mashm), value :: in_mashm
    type(Mashm), value :: in_mashm
    integer(c_int), value :: commIndex
    integer(c_int) :: msgSize
  end function

  subroutine MashmCommFinish(in_mashm) &
    bind(c, name='MashmCommFinish')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    !Mashm, value :: in_mashm
    !type(Mashm), value :: in_mashm
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmInterNodeCommBegin(in_mashm) &
    bind(c, name='MashmInterNodeCommBegin')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmInterNodeCommEnd(in_mashm) &
    bind(c, name='MashmInterNodeCommEnd')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmIntraNodeCommBegin(in_mashm) &
    bind(c, name='MashmIntraNodeCommBegin')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
  end subroutine 

  subroutine MashmIntraNodeCommEnd(in_mashm) &
    bind(c, name='MashmIntraNodeCommEnd')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
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

  subroutine MashmGetBufferPointer2C(in_mashm, i, sendReceive, bufferPointer) &
    bind(c, name='MashmGetBufferPointer2')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value, intent(in) :: in_mashm
    integer(c_int), value, intent(in) :: i
    integer(c_int), value, intent(in) :: sendReceive
    type(c_ptr), intent(out) :: bufferPointer
  end subroutine

  function MashmIsMsgIntraNodalC(in_mashm, commIndex) result(isOnNode) &
    bind(c, name='MashmIsMsgIntraNodal')
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, value :: commIndex
    integer :: isOnNode
   
  end function

  subroutine MashmRetireBufferPointerC(in_mashm, bufferPointer)
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value, intent(in) :: in_mashm
    type(c_ptr), intent(out) :: bufferPointer

  end subroutine

end interface

contains

  subroutine MashmGetBufferPointer(in_mashm, i, sendReceive, ftnBufferPointer)
    use, intrinsic :: iso_c_binding
    use Mashm_enum_mod
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, intent(in) :: i
    integer(c_int), intent(in) :: sendReceive
    type(MashmBufferPointer), intent(inout) :: ftnBufferPointer
    integer(c_int) :: msgSize
    integer(c_int) :: iBaseZero
    integer :: fMsgSize
    iBaseZero = i - 1

    ! TODO: move to MashmGetBufferPointerC function call - not working for some
    !       reason
    call MashmGetBufferPointer2C(in_mashm, iBaseZero, sendReceive, ftnBufferPointer%cPtr)
    !ftnBufferPointer%cPtr = MashmGetBufferPointerC(in_mashm, iBaseZero, sendReceive)

    msgSize = MashmGetCommSize(in_mashm, iBaseZero)
    fMsgSize = msgSize

    call c_f_pointer(cptr=ftnBufferPointer%cPtr,fptr=ftnBufferPointer%p,shape = (/ fMsgSize /))

  end subroutine

  subroutine MashmRetireBufferPointer(in_mashm, ftnBufferPointer)
    use, intrinsic :: iso_c_binding
    use Mashm_enum_mod
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    type(MashmBufferPointer), intent(inout) :: ftnBufferPointer

    ! Nullify the Fortran pointer
    nullify(ftnBufferPointer%p)

    ! Retire (nullify) the C pointer
    call MashmRetireBufferPointerC(in_mashm, ftnBufferPointer%cPtr)
  end subroutine


  subroutine MashmSetComm(in_mashm, commIndex, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, value :: commIndex
    integer(c_int), value :: pairRank
    integer(c_int), value :: msgSize
    integer(c_int) :: commIndexBaseZero
    commIndexBaseZero = commIndex - 1
    call MashmSetCommC(in_mashm, commIndexBaseZero, pairRank, msgSize)
  end subroutine 

  function MashmIsMsgIntraNodal(in_mashm, commIndex) result(isOnNode) 
    use, intrinsic :: iso_c_binding
    use Mashm_type
    implicit none
    type(Mashm), value :: in_mashm
    integer, value :: commIndex
    logical :: isOnNode
    isOnNode = ( MashmIsMsgIntraNodalC(in_mashm, int(commIndex - 1,kind=c_int)) .eq. 1 )

  end function
end module Mashm_mod

