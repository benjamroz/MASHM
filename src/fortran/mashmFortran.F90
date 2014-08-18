
module mashm_mod
#include "mashmf.h"
use, intrinsic :: iso_c_binding
implicit none
interface

  subroutine mashmInit(in_mashm, in_comm)
    use, intrinsic :: iso_c_binding
    Mashm :: in_mashm
    integer, value :: in_comm
  end subroutine

  subroutine mashmDestroy(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm :: in_mashm
  end subroutine

  subroutine mashmPrintInfo(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine mashmAddSymComm(in_mashm, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: pairRank
    integer, value :: msgSize
  end subroutine 

  subroutine mashmCommFinish(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine mashmInterNodeCommBegin(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine mashmInterNodeCommEnd(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine mashmIntraNodeCommBegin(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

  subroutine mashmIntraNodeCommEnd(in_mashm)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
  end subroutine 

end interface

end module mashm_mod
