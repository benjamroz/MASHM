
module mashm_mod
#include "mashmf.h"
implicit none
interface

  subroutine mashmInit(in_mashm, in_comm)
    Mashm in_mashm
    integer in_comm
  end subroutine

  subroutine mashmDestroy(in_mashm)
    Mashm in_mashm
  end subroutine

end interface

end module mashm_mod
