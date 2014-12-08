module arrayOfPointers_mod
  type arrayOfPointers
    real*8, allocatable :: p(:)
  end type arrayOfPointers
end module 


program nodalCommFtn
#include "Mashmf.h"
use iso_c_binding
use Mashm_mod
use mpi
use arrayOfPointers_mod
use grid_data
implicit none

integer :: ierr
integer :: rank, numProcs
integer :: numElems
integer, allocatable :: gridIndices(:,:)

Mashm :: myMashm

call MPI_Init(ierr)
call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

! Read namelist
call read_grid_data_namelist(MPI_COMM_WORLD)


! Get the grid decomposition
call grid_3d_decomp_num_elements(rank, numProcs)

!allocate(gridIndices(3,numElems))

! 
!call grid_3d_decomp_get_elements(rank, numProcs)



!deallocate(gridIndices)

call MPI_Finalize(ierr)

end program nodalCommFtn


