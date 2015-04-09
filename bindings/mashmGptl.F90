module mashmGptl
use gptl
implicit none
contains
  function mashmGptlinitialize() result(gptlErrorCode)
  implicit none
  integer :: gptlErrorCode
  gptlErrorCode = gptlinitialize()
  end function 

  function mashmGptlsetutr(gptlTimer) result(gptlErrorCode)
  implicit none
  integer, intent(in) :: gptlTimer
  integer :: gptlErrorCode
  gptlErrorCode = gptlsetutr(gptlTimer)
  end function

  function mashmGptlstart(timerString) result(gptlErrorCode)
  implicit none
  character(*), intent(in) :: timerString
  integer :: gptlErrorCode
  gptlErrorCode = gptlstart(timerString)
  end function 

  function mashmGptlstop(timerString) result(gptlErrorCode)
  implicit none
  character(*), intent(in) :: timerString
  integer :: gptlErrorCode
  gptlErrorCode = gptlstop(timerString)
  end function

  function mashmGptlpr_summary_file(mpiComm,fileName) result(gptlErrorCode)
  implicit none
  integer, intent(in) :: mpiComm
  character(*), intent(in) :: fileName
  integer :: gptlErrorCode
  gptlErrorCode = gptlpr_summary_file(mpiComm,fileName)
  end function 

  function mashmGptlpr_file(fileName) result(gptlErrorCode)
  implicit none
  character(*), intent(in) :: fileName
  integer :: gptlErrorCode
  gptlErrorCode = gptlpr_file(fileName)
  end function

end module mashmGptl

