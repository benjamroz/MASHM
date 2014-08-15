

program nodalCommFtn
#include "mashmf.h"
use mashm_mod
use mpi
implicit none
Mashm :: myMashm
integer :: ierr
integer :: myComm

call MPI_Init(ierr)
print *, "ierr = ", ierr
myComm = MPI_COMM_WORLD
call mashmInit(myMashm, MPI_COMM_WORLD)

!call mashmDestroy(myMashm)

call MPI_Finalize(ierr)

end program nodalCommFtn
