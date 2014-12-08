
module grid_data
  use mpi
  private 
  integer :: nx = -1
  integer :: ny = -1
  integer :: nz = -1
  real*8 :: lx = -1
  real*8 :: ly = -1
  real*8 :: lz = -1
  integer :: gridComm
  integer :: gridRank
  integer :: gridNumProcs

  real*8 :: hx = -1
  real*8 :: hy = -1
  real*8 :: hz = -1

  integer :: numLocalElems
  integer :: numTotalElems
  real*8, allocatable :: coordinates(:,:)
  real*8, allocatable :: rhs(:), lhs(:)
  integer, allocatable :: cellIndices(:,:)

  public :: read_grid_data_namelist
  public :: grid_3d_decomp_num_elements
  public :: grid_3d_decomp_get_elements

contains

  function get_num_cells(iDir) result(numCells)
    integer, intent(in) :: iDir
    integer :: numCells
    select case (iDir)
      case (1)
        numCells = nx
      case (2)
        numCells = ny
      case (3)
        numCells = nz
    end select
  end function

  function get_grid_length(iDir) result(dirLength)
    integer, intent(in) :: iDir
    real*8 :: dirLength
    select case (iDir)
      case (1)
        dirLength = lx
      case (2)
        dirLength = ly
      case (3)
        dirLength = lz
    end select
  end function

  subroutine print_grid_data()
    

  end subroutine

  subroutine read_grid_data_namelist(comm)
    implicit none
    integer, intent(in) :: comm
    real*8, parameter :: PI = 3.1415926535897932385
    integer :: ierr
    namelist /grid_data/ nx, ny, nz, lx, ly, lz

    gridComm = comm
    call MPI_Comm_size(gridComm, gridNumProcs, ierr)
    call MPI_Comm_rank(gridComm, gridRank, ierr)

    ! Set default values
    nx = 16
    ny = 16
    nz = 16

    lx = 2.0*PI
    ly = 2.0*PI
    lz = 2.0*PI

    if (gridRank .eq. 0) then 
      read(*,nml=grid_data)

      write(*,*) "Read grid_data namelist: "
      write(*,nml=grid_data)
    endif


    ! Broadcast the data to the other ranks
    call MPI_Bcast(nx, 1, MPI_INTEGER, 0, gridComm, ierr)
    call MPI_Bcast(ny, 1, MPI_INTEGER, 0, gridComm, ierr)
    call MPI_Bcast(nz, 1, MPI_INTEGER, 0, gridComm, ierr)
    call MPI_Bcast(lx, 1, MPI_REAL8, 0, gridComm, ierr)
    call MPI_Bcast(ly, 1, MPI_REAL8, 0, gridComm, ierr)
    call MPI_Bcast(lz, 1, MPI_REAL8, 0, gridComm, ierr)

  end subroutine read_grid_data_namelist

  subroutine setup_grid_data(rank, numProcs)
    implicit none
    integer, intent(in) :: rank, numProcs

    integer :: modIndex
    integer :: numTotalCells, numLocalCells
    hx = lx / nx
    hy = ly / ny
    hz = lz / nz
   
    numTotalCells = nx*ny*nz
    ! Determine the number of cells each MPI task owns
   
    modIndex = mod(numTotalCells, numProcs)
    if (rank < modIndex) then
      numLocalCells = numTotalCells/numProcs + 1
    else
      numLocalCells = numTotalCells/numProcs
    endif 

    ! allocate the ... data
    allocate(rhs(numLocalCells))
    allocate(lhs(numLocalCells))
    rhs = 0.0
    lhs = 0.0

    ! Determine connectivity

  end subroutine 

  subroutine grid_3d_decomp_num_elements(rank, numProcs)
  implicit none
  integer, intent(in) :: rank, numProcs

  integer :: factors(3), remainder
  integer :: cellsInDir(3), localCellsInDir
  integer :: iDir, iRank
  integer :: startIndex(3), endIndex(3)
  integer :: ierr


  ! Determine the domain size 
  numTotalElems = nx*ny*nz

  print *, "nx, ny, nz, numTotalElems = ", nx, ny, nz, numTotalElems
  if (mod(numTotalElems, numProcs) .eq. 0) then
    numLocalElems = numTotalElems/numProcs
  else
    print *, "Error: nx * ny * nz = ", numTotalElems, " not divisible by the number of processes ", numProcs
  endif

  allocate(coordinates(3,numLocalElems))
  allocate(rhs(numLocalElems))
  allocate(lhs(numLocalElems))
  allocate(cellIndices(3,numLocalElems))

  ! Create rectilinear blocks of cells.
  cellsInDir(1) = nx
  cellsInDir(2) = ny
  cellsInDir(3) = nz

  remainder = numLocalElems
  factors(:) = 1
  do while (.true.)
    do iDir = 1,3
      if (mod(remainder, 2) .eq. 0) then
        factors(iDir) = factors(iDir)*2
        remainder = remainder/2
      endif
    enddo
    if (remainder .eq. 1) then
      exit 
    endif
  enddo

  if (rank .eq. 0) then
    print *, "factors = ", factors
  endif

  startIndex(1) = 1
  startIndex(2) = 1
  startIndex(3) = 1

  startIndex(2) = startIndex(2) + (startIndex(1) + gridRank*factors(1))/nx*factors(2)

  startIndex(1) = mod(startIndex(1) + gridRank*factors(1),nx)
  startIndex(3) = startIndex(3) + (startIndex(2))/ny*factors(3)
  startIndex(2) = mod(startIndex(2), ny)

  endIndex(1) = startIndex(1) + factors(1) - 1
  endIndex(2) = startIndex(2) + factors(2) - 1
  endIndex(3) = startIndex(3) + factors(3) - 1

  call MPI_Barrier(gridComm, ierr)
  call MPI_Barrier(gridComm, ierr)
  call flush(6)
  do iRank = 0, gridNumProcs - 1
    if (iRank .eq. gridRank) then
      print *, "rank: ", gridRank
      print *, "  start: ", startIndex
      print *, "  end: ", endIndex
    endif
    call MPI_Barrier(gridComm, ierr)
    call MPI_Barrier(gridComm, ierr)
    call flush(6)
    call MPI_Barrier(gridComm, ierr)
    call MPI_Barrier(gridComm, ierr)
  enddo

  end subroutine grid_3d_decomp_num_elements


  subroutine grid_3d_decomp_get_elements(rank, numProcs, gridIndices)
  implicit none
  integer, intent(in) :: rank, numProcs
  integer, intent(out) :: gridIndices(:,:)


  end subroutine grid_3d_decomp_get_elements

end module 

