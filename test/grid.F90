
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
  integer, target :: startIndex(3), endIndex(3)
  integer :: factors(3)

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
  public :: grid_3d_get_indices
  public :: getElemRank
  public :: determineCommSchedule
  public :: get_num_cells

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

  integer :: remainder
  integer :: cellsInDir(3), localCellsInDir
  integer :: iDir, iRank
  integer :: ierr
  integer :: counter


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
  counter = 0
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
    counter = counter + 1
    if (counter .gt. 100) then
      print *, "Error: nx,ny,nz must be factors of 2"
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


  subroutine grid_3d_get_indices(numElems, gridIndicesStart, gridIndicesEnd)
  implicit none
  integer, intent(out) :: numElems
  integer, pointer, intent(out) :: gridIndicesStart(:), gridIndicesEnd(:)

  numElems = numLocalElems
  gridIndicesStart => startIndex
  gridIndicesEnd => endIndex

  end subroutine 

  function getElemRank(ix, iy, iz) result(iRank)
  implicit none
  integer, intent(in) :: ix, iy, iz
  integer :: iRank

  integer :: facX, facY, facZ

  if (ix < 1 .or. ix > nx .or. &
      iy < 1 .or. iy > ny .or. &
      iz < 1 .or. iz > nz) then
    iRank = -1
  else
    facX = (ix-1)/factors(1)
    facY = (iy-1)/factors(2)
    facZ = (iz-1)/factors(3)

    iRank = facX + nx/factors(1)*facY + nx/factors(1)*ny/factors(2)*facZ
  endif

  end function

  subroutine determineCommSchedule(rank, numProcs, gridIndicesStart, gridIndicesEnd, &
                                   numMessages, msgIndices, msgSizes, msgOffsets, neighborRanks)
  implicit none
  integer, intent(in) :: rank, numProcs
  integer, intent(in) :: gridIndicesStart(3), gridIndicesEnd(3)

  integer, intent(out) :: numMessages
  integer, intent(out) :: msgIndices(-1:1,-1:1,-1:1)
  integer, allocatable, intent(out) :: msgSizes(:), msgOffsets(:)
  integer, allocatable, intent(out) :: neighborRanks(:)

  integer :: iRank, ix, iy, iz, elemRank, rankCounter
  integer :: ierr
  integer :: elemRankCounter(0:numProcs-1)
  integer :: tmpDir(-1:1,-1:1,-1:1)
  integer :: tmpDirX, tmpDirY, tmpDirZ
  integer :: totalMessageSize, i
  logical :: found

  do iRank = 0, numProcs - 1
    if (iRank == rank) then
      print *, "Rank: start ", gridIndicesStart, ", end ", gridIndicesEnd
      do iz = gridIndicesStart(3), gridIndicesEnd(3)
        do iy = gridIndicesStart(2), gridIndicesEnd(2)
          do ix = gridIndicesStart(1), gridIndicesEnd(1)
            elemRank = getElemRank(ix,iy,iz)
            print *, "  Element (", ix, ",", iy, ",", iz, ") belongs to rank ", elemRank
          enddo
        enddo
      enddo
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call flush(6)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  enddo

  elemRankCounter = 0
  ! Determine number of neighbors
  do iz = gridIndicesStart(3) - 1, gridIndicesEnd(3) + 1
    do iy = gridIndicesStart(2) - 1, gridIndicesEnd(2) + 1
      do ix = gridIndicesStart(1) - 1, gridIndicesEnd(1) + 1
        ! If a corner - then cycle
        if (    (ix < gridIndicesStart(1) .or. ix > gridIndicesEnd(1)) &
          .and. (iy < gridIndicesStart(2) .or. iy > gridIndicesEnd(2)) &
          .and. (iz < gridIndicesStart(3) .or. iz > gridIndicesEnd(3)) ) then
          cycle
        endif
        elemRank = getElemRank(ix,iy,iz)
        if (elemRank .ne. -1) then
          tmpDir = 0
          elemRankCounter(elemRank) = 1

        endif
      enddo
    enddo
  enddo

  rankCounter = 0
  do iRank = 0, numProcs - 1
    if (elemRankCounter(iRank) .gt. 0) then
      rankCounter = rankCounter + 1
    endif
  enddo

  ! Subtract one for self
  rankCounter = rankCounter - 1
  numMessages = rankCounter

  allocate(neighborRanks(numMessages))
  ! Determine neighbor ranks
  msgIndices = 0
  rankCounter = 0
  do iz = gridIndicesStart(3) - 1, gridIndicesEnd(3) + 1
    do iy = gridIndicesStart(2) - 1, gridIndicesEnd(2) + 1
      do ix = gridIndicesStart(1) - 1, gridIndicesEnd(1) + 1
        ! If a corner - then cycle
        if (    (ix < gridIndicesStart(1) .or. ix > gridIndicesEnd(1)) &
          .and. (iy < gridIndicesStart(2) .or. iy > gridIndicesEnd(2)) &
          .and. (iz < gridIndicesStart(3) .or. iz > gridIndicesEnd(3)) ) then
          cycle
        endif
        elemRank = getElemRank(ix,iy,iz)
        if (elemRank .ne. -1 .and. elemRank .ne. rank) then
          found = .false.
          do iRank = 1, rankCounter
            if (neighborRanks(iRank) .eq. elemRank) then
              found = .true.
              exit
            endif
          enddo
          if (.not. found) then
            rankCounter = rankCounter + 1
            neighborRanks(rankCounter) = elemRank 
            tmpDirZ = 0
            tmpDirY = 0
            tmpDirX = 0
            if (iz .gt. gridIndicesEnd(3)) then
              tmpDirZ = 1
            else if (iz .lt. gridIndicesStart(3)) then
              tmpDirZ = -1
            endif
            if (iy .gt. gridIndicesEnd(2)) then
              tmpDirY = 1
            else if (iy .lt. gridIndicesStart(2)) then
              tmpDirY = -1
            endif
            if (ix .gt. gridIndicesEnd(1)) then
              tmpDirX = 1
            else if (ix .lt. gridIndicesStart(1)) then
              tmpDirX = -1
            endif
            msgIndices(tmpDirX,tmpDirY,tmpDirZ) = rankCounter
            

          endif
        endif
      enddo
    enddo
  enddo

  numMessages = numMessages
  allocate(msgSizes(numMessages))
  allocate(msgOffsets(numMessages+1))
  msgSizes = 0
  ! Determine message sizes
  do iz = gridIndicesStart(3) - 1, gridIndicesEnd(3) + 1
    do iy = gridIndicesStart(2) - 1, gridIndicesEnd(2) + 1
      do ix = gridIndicesStart(1) - 1, gridIndicesEnd(1) + 1
        ! If a corner - then cycle
        if (    (ix < gridIndicesStart(1) .or. ix > gridIndicesEnd(1)) &
          .and. (iy < gridIndicesStart(2) .or. iy > gridIndicesEnd(2)) &
          .and. (iz < gridIndicesStart(3) .or. iz > gridIndicesEnd(3)) ) then
          cycle
        endif
        elemRank = getElemRank(ix,iy,iz)
        if (elemRank .ne. -1 .and. elemRank .ne. rank) then
          do iRank = 1, numMessages
            if (neighborRanks(iRank) .eq. elemRank) then
              msgSizes(iRank) = msgSizes(iRank) + 1
              exit
            endif
          enddo
        endif
      enddo
    enddo
  enddo

  totalMessageSize = 0
  do i = 1, numMessages
    msgOffsets(i) = totalMessageSize
    totalMessageSize = totalMessageSize + msgSizes(i)
  enddo
  msgOffsets(numMessages+1) = totalMessageSize

  do iRank = 0, numProcs - 1
    if (iRank == rank) then
      print *, "Rank ", rank, " has ", numMessages, " neighbors"
      do i = 1, numMessages
        print *, "  pair Rank ", neighborRanks(i), ", size ", msgSizes(i)
      enddo
      print *, "  totalSize ", totalMessageSize
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call flush(6)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  enddo

  end subroutine

end module 

