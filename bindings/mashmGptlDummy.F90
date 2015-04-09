module gptl
  integer :: gptlnanotime
contains 
  function gptlinitialize() result(gptlErrorCode)
  implicit none
  integer :: gptlErrorCode
  gptlErrorCode = 0
  end function 

  function gptlsetutr(gptlTimer) result(gptlErrorCode)
  implicit none
  integer, intent(in) :: gptlTimer
  integer :: gptlErrorCode
  gptlErrorCode = 0
  end function gptlsetutr

  function gptlstart(timerString) result(gptlErrorCode)
  implicit none
  character(*), intent(in) :: timerString
  integer :: gptlErrorCode
  gptlErrorCode = 0
  end function 

  function gptlstop(timerString) result(gptlErrorCode)
  implicit none
  character(*), intent(in) :: timerString
  integer :: gptlErrorCode
  gptlErrorCode = 0
  end function

  function gptlpr_summary_file(mpiComm,fileName) result(gptlErrorCode)
  implicit none
  integer, intent(in) :: mpiComm
  character(*), intent(in) :: fileName
  integer :: gptlErrorCode
  gptlErrorCode = 0
  end function 

  function gptlpr_file(fileName) result(gptlErrorCode)
  implicit none
  character(*), intent(in) :: fileName
  integer :: gptlErrorCode
  gptlErrorCode = 0
  end function

end module gptl

