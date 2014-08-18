
module Mashm_mod
#include "Mashmf.h"
implicit none
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

  subroutine MashmAddSymComm(in_mashm, pairRank, msgSize)
    use, intrinsic :: iso_c_binding
    Mashm, value :: in_mashm
    integer, value :: pairRank
    integer, value :: msgSize
  end subroutine 

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

end interface

end module Mashm_mod
