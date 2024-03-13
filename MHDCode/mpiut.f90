  module mpiut
! -----------------------------------------------------------
!                    MODULE mpiut
! -----------------------------------------------------------
  use common
  use mpi
  implicit none

  integer :: ierr
  integer :: dims(2), coords(2)
  logical :: periods(2), reorder  
  public
  private :: exchange_x1, &
             exchange_x3, &
             exchange_y1, &
             exchange_y3

  interface mpiut__exchange_x
  module procedure exchange_x1, &
                    exchange_x3
  end interface

  interface mpiut__exchange_y
  module procedure exchange_y1, &
                   exchange_y3
  end interface


  contains

! ======================================================
  subroutine mpiut__init
! ======================================================
  call mpi_init(ierr)

! ---------------------
! Original communicator 
! ---------------------
  call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)

! --------------- -------
! Cartesion communicator 
! ----------------------
  dims(1) = nproc_x
  dims(2) = nproc_y

  periods(1) = .true.
  periods(2) = .true.

  reorder = .true.

  call mpi_cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,MPI_COMM_CART,ierr)

  call mpi_comm_size(MPI_COMM_CART,nproc,ierr)
  call mpi_comm_rank(MPI_COMM_CART,myrank,ierr)
!
  if(nproc .ne. nproc_x*nproc_y) then
     write(*,*) '::ERR:mpiut__init, nproc not equal to nproc_x*nproc_y', &
                nproc, nproc_x, nproc_y
     call mpiut__finalize
     stop
  end if

! -----------------------
! My Coordinate 
! -----------------------
  call MPI_CART_COORDS(MPI_COMM_CART, myrank, 2, coords, ierr)
  index_x = coords(1)
  index_y = coords(2)

! -----------------------
! Coordinate of neighbors 
! -----------------------
  call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, rank_left, rank_right, ierr)
  call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, rank_down, rank_up, ierr)

  call mpiut__barrier 
  end subroutine mpiut__init

! =============================================================
  subroutine mpiut__finalize
! =============================================================
  call mpi_finalize(ierr)

  end subroutine mpiut__finalize

! =============================================================
  subroutine mpiut__barrier
! =============================================================
  call MPI_BARRIER(MPI_COMM_CART,ierr)

  end subroutine mpiut__barrier

! =====================================================================
  subroutine mpiut__max(d)
! =====================================================================
  real(DP) :: sendbuf, recvbuf, d

  sendbuf = d
  call MPI_REDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                  root, MPI_COMM_CART, ierr)
  call MPI_BCAST(recvbuf, 1, MPI_DOUBLE_PRECISION, &
                 root, MPI_COMM_CART, ierr)
  d = recvbuf

  end subroutine mpiut__max

! ==============================================================
  subroutine mpiut__min(d)
! ==============================================================
  real(DP) :: sendbuf, recvbuf, d

  sendbuf = d
  call MPI_REDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                  root, MPI_COMM_CART, ierr)

  call MPI_BCAST(recvbuf, 1, MPI_DOUBLE_PRECISION, &
                 root, MPI_COMM_CART, ierr)
  d = recvbuf

  end subroutine mpiut__min

! ===============================================================
  subroutine mpiut__sum(d)
! ===============================================================
  real(DP), intent(inout) :: d
  real(DP) :: sendbuf, recvbuf

  sendbuf = d
  call MPI_REDUCE(sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                  root, MPI_COMM_CART, ierr)

  call MPI_BCAST(recvbuf, 1, MPI_DOUBLE_PRECISION, &
                 root, MPI_COMM_CART, ierr)
  d = recvbuf

  end subroutine mpiut__sum

! ================================================================
  subroutine exchange_x1(ax)
! ================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(inout)  :: ax
  real(DP), dimension(0:NZ,-1:NY,1) :: sendbuf, recvbuf
  integer :: sendcount = (NY+2)*(NZ+1)*1
  integer :: recvcount = (NY+2)*(NZ+1)*1
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------------
! from left to right
! ------------------
  sendbuf(:,:,1) = ax(:,:,NX-1)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_right, 0, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_left, 0,  &
                    MPI_COMM_CART, status, ierr)
  ax(:,:,-1) = recvbuf(:,:,1)

! ------------------ 
! from right to left
! ------------------
  sendbuf(:,:,1) = ax(:,:,0)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_left, 1,  & 
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_right, 1, &
                    MPI_COMM_CART, status, ierr)
  ax(:,:,NX) = recvbuf(:,:,1)

  end subroutine exchange_x1

! ===================================================================
  subroutine exchange_x3(ax,ay,az)
! ===================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(inout)  :: ax, ay, az
  real(DP), dimension(0:NZ,-1:NY,3) :: sendbuf, recvbuf
  integer :: sendcount = (NY+2)*(NZ+1)*3
  integer :: recvcount = (NY+2)*(NZ+1)*3
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------------
! from left to right
! ------------------
  sendbuf(:,:,1) = ax(:,:,NX-1)
  sendbuf(:,:,2) = ay(:,:,NX-1)
  sendbuf(:,:,3) = az(:,:,NX-1)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_right, 0, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_left, 0,  &
                    MPI_COMM_CART, status, ierr)
  ax(:,:,-1) = recvbuf(:,:,1)
  ay(:,:,-1) = recvbuf(:,:,2)
  az(:,:,-1) = recvbuf(:,:,3)

! ------------------
! from right to left
! ------------------
  sendbuf(:,:,1) = ax(:,:,0)
  sendbuf(:,:,2) = ay(:,:,0)
  sendbuf(:,:,3) = az(:,:,0)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_left, 1,  &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_right, 1, &
                    MPI_COMM_CART, status, ierr)
  ax(:,:,NX) = recvbuf(:,:,1)
  ay(:,:,NX) = recvbuf(:,:,2)
  az(:,:,NX) = recvbuf(:,:,3)

  return
  end subroutine exchange_x3

! =====================================================================
  subroutine exchange_y1(ax)
! =====================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(inout)  :: ax
  real(DP), dimension(0:NZ,-1:NX,1) :: sendbuf, recvbuf
  integer :: sendcount = (NX+2)*(NZ+1)*1
  integer :: recvcount = (NX+2)*(NZ+1)*1
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ---------------
! from down to up
! ---------------
  sendbuf(:,:,1) = ax(:,NY-1,:)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_up, 0,   &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_down, 0, &
                    MPI_COMM_CART, status, ierr)
  ax(:,-1,:) = recvbuf(:,:,1)

! ---------------
! from up to down
! ---------------
  sendbuf(:,:,1) = ax(:,0,:)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_down, 1, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_up, 1,   &
                    MPI_COMM_CART, status, ierr)
   ax(:,NY,:) = recvbuf(:,:,1)

   end subroutine exchange_y1

! =====================================================================
  subroutine exchange_y3(ax,ay,az)
! =====================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(inout)  :: ax, ay, az
  real(DP), dimension(0:NZ,-1:NX,3) :: sendbuf, recvbuf
  integer :: sendcount = (NX+2)*(NZ+1)*3
  integer :: recvcount = (NX+2)*(NZ+1)*3
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ---------------
! from down to up
! ---------------
  sendbuf(:,:,1) = ax(:,NY-1,:)
  sendbuf(:,:,2) = ay(:,NY-1,:)
  sendbuf(:,:,3) = az(:,NY-1,:)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_up, 0,   &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_down, 0, &
                    MPI_COMM_CART, status, ierr)

  ax(:,-1,:) = recvbuf(:,:,1)
  ay(:,-1,:) = recvbuf(:,:,2)
  az(:,-1,:) = recvbuf(:,:,3)

! ---------------
! from up to down
! ---------------
  sendbuf(:,:,1) = ax(:,0,:)
  sendbuf(:,:,2) = ay(:,0,:)
  sendbuf(:,:,3) = az(:,0,:)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_down, 1, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_up, 1,   &
                    MPI_COMM_CART, status, ierr)
  ax(:,NY,:) = recvbuf(:,:,1)
  ay(:,NY,:) = recvbuf(:,:,2)
  az(:,NY,:) = recvbuf(:,:,3)

  end subroutine exchange_y3

!======================= public =======================================

! =====================================================================
  subroutine mpiut__exchange_2d_x(ax)
! =====================================================================
  real(DP), dimension(-1:NY,-1:NX), intent(inout)  :: ax
  real(DP), dimension(-1:NY,1) :: sendbuf, recvbuf
  integer :: sendcount = (NY+2)*1
  integer :: recvcount = (NY+2)*1
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ------------------
! from left to right
! ------------------
  sendbuf(:,1) = ax(:,NX-1)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_right, 0, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_left, 0,  &
                    MPI_COMM_CART, status, ierr)
  ax(:,-1) = recvbuf(:,1)

! ------------------ 
! from right to left
! ------------------
  sendbuf(:,1) = ax(:,0)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_left, 1,  &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_right, 1, &
                    MPI_COMM_CART, status, ierr)
  ax(:,NX) = recvbuf(:,1)

  end subroutine mpiut__exchange_2d_x

! =====================================================================
  subroutine mpiut__exchange_2d_y(ax)
! =====================================================================
  real(DP), dimension(-1:NY,-1:NX), intent(inout)  :: ax
  real(DP), dimension(-1:NX,1) :: sendbuf, recvbuf
  integer :: sendcount = (NX+2)*1
  integer :: recvcount = (NX+2)*1
  integer, dimension(MPI_STATUS_SIZE) :: status 

! ---------------
! from down to up
! ---------------
  sendbuf(:,1) = ax(NY-1,:)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_up, 0,   &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_down, 0, &
                    MPI_COMM_CART, status, ierr)
  ax(-1,:) = recvbuf(:,1)

! ---------------
! from up to down
! ---------------
  sendbuf(:,1) = ax(0,:)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_down, 1, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_up, 1,   &
                    MPI_COMM_CART, status, ierr)
   ax(NY,:) = recvbuf(:,1)

  end subroutine mpiut__exchange_2d_y

! =====================================================================
  subroutine mpiut__exchange_2dxl_y(ax)
! =====================================================================
  real(DP), dimension( 0:NZ,-1:NY), intent(inout)  :: ax
  real(DP), dimension( 0:NZ,1)                     :: sendbuf, recvbuf
  integer                                          :: sendcount = (NY+2)*1
  integer                                          :: recvcount = (NY+2)*1
  integer, dimension(MPI_STATUS_SIZE)              :: status 

! ---------------
! from down to up
! ---------------
  sendbuf(:,1) = ax(:,NY-1)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_up, 0,   &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_down, 0, &
                    MPI_COMM_CART, status, ierr)
  ax(:,-1) = recvbuf(:,1)

! ---------------
! from down to up
! ---------------
  sendbuf(:,1) = ax(:,0)
  call MPI_SENDRECV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, rank_down, 1, &
                    recvbuf, recvcount, MPI_DOUBLE_PRECISION, rank_up, 1,   &
                    MPI_COMM_CART, status, ierr)
   ax(:,NY) = recvbuf(:,1)

  end subroutine mpiut__exchange_2dxl_y
  end module mpiut










