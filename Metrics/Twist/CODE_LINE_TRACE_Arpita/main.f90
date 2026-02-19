! **************************************************************************
  program bline
! **************************************************************************

! ----------------------------------------------------------------------------
! Initial set up
! ----------------------------------------------------------------------------
  include "common.f90"
  integer :: i, j, k,l,numbl, isol=1, idir, line, dim_lines, line_each_dim
! isol = 1: Euler,  2: RKG
!
! filed line starting points are determined as
! a point for dim_lines=0
! equi-interval poins on a line for dim_lines = 1
! grid points on a plane for dim_lines        = 2
! grid points in a cube for dim_lines         = 3

  character clp*3,filen*11
  real*8,dimension(0:NX,0:NY,0:NZ) :: a3d
  real(8) :: r_start(3), x0, y0, z0
  real(8) :: dxl, dyl, dzl
  real(8) :: r_0(3),r_a(3),r_b(3),r_c(3)
  real(8) :: r_0a(3),r_0b(3),r_0c(3)
  real(8) :: alpha, beta, gamma
  real(8) :: del_alpha, del_beta, del_gamma
  
  character(3) :: cloop  
  integer      :: loop

! -----------------------------------------------------------------------------
! Read coordinate data input from AVS field data 
! -----------------------------------------------------------------------------
  open(8,file=dir_r//'coord.xgc',status='OLD')
  do i = 0, NX
     read(8,*) xc0(i)
  end do
  close(8)
!
  open(8,file=dir_r//'coord.ygc',status='OLD')
  do j = 0, NY
     read(8,*) yc0(j)
  end do
  close(8)
!
  open(8,file=dir_r//'coord.zgc',status='OLD')
  do k = 0, NZ
     read(8,*) zc0(k)
  end do
!
  do i = 0, nxl
     xc(i) = xc0(nx0+i*leapx)
  end do

  do j = 0, nyl
     yc(j) = yc0(ny0+j*leapy)
  end do
!
  do k = 0, nzl
     zc(k) = zc0(nz0+k*leapz)
  end do

! -----------------------------------------------------------------------------
! Read 3D field from AVS field data
! -----------------------------------------------------------------------------
  open(10,file=dir_r//'B3D.001.BX.R.014',form='UNFORMATTED',status='OLD')
  read(10) a3d
  close(10)
  call reduced_array(a3d,bx)
!
  open(10,file=dir_r//'B3D.001.BY.R.014',form='UNFORMATTED',status='OLD')
  read(10) a3d
  close(10)
  call reduced_array(a3d,by)
!
  open(10,file=dir_r//'B3D.001.BZ.R.014',form='UNFORMATTED',status='OLD')
  read(10) a3d
  close(10)
  call reduced_array(a3d,bz)
!
  open(10,file=dir_r//'B3D.001.CX.R.014',form='UNFORMATTED',status='OLD')
  read(10) a3d
  close(10)
  call reduced_array(a3d,cx)
!
  open(10,file=dir_r//'B3D.001.CY.R.014',form='UNFORMATTED',status='OLD')
  read(10) a3d
  close(10)
  call reduced_array(a3d,cy)
!
  open(10,file=dir_r//'B3D.001.CZ.R.014',form='UNFORMATTED',status='OLD')
  read(10) a3d
  close(10)
  call reduced_array(a3d,cz)
!
!
! -----------------------------------------------------------------------------
! Initial Point
! -----------------------------------------------------------------------------
      open(10,file='Param_bline_2D',status='OLD')
 500  read(10,*,end=999) dim_lines,numbl,rbl,gbl,bbl,wbl

! -----------------------------------------------------------------------------
! filed line tracing
! -----------------------------------------------------------------------------
  call preprg

! -----------------------------------------------------------------------------
! POINT
! -----------------------------------------------------------------------------
  if(dim_lines.eq.0) then
     read(10,*,end=999) r_start(1:3)
  l = 1
  idir = 1
  call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
  idir = -1
  call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
 
! -----------------------------------------------------------------------------
! LINE  
! -----------------------------------------------------------------------------
  else if(dim_lines.eq.1) then
  read(10,*,end=999) r_0(1:3), r_0a(1:3)

  do l = 1, numbl
     r_start(:) = r_0(:) + (l-1)*r_0a(:)

  idir = 1
  call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
  idir = -1
  call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
  end do

! ------------------------------------------------------------------------------
! PLANE
! -------------------------------------------------------------------------------
  else if(dim_lines.eq.2) then
  read(10,*,end=999) r_0(1:3), r_0a(1:3), r_0b(1:3)
  line_each_dim = max(2, int(real(numbl)**0.5))
  del_alpha = 1.0d0
  del_beta  = 1.0d0

  l = 0
  do i = 1, line_each_dim
     alpha = del_alpha*(i-1)
  do j = 1, line_each_dim
     beta  = del_beta *(j-1)
     r_start(:) = r_0(:) + alpha*r_0a(:) + beta*r_0b(:)
!
     l = l + 1
     idir = 1
     call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
     idir = -1
     call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
   end do
   end do

! ------------------------------------------------------------------------------
! CUBE 
! ------------------------------------------------------------------------------
  else if(dim_lines.eq.3) then
  read(10,*,end=999) r_0(1:3), r_0a(1:3), r_0b(1:3), r_0c(1:3)
  line_each_dim = max(2, int(real(numbl)**(1.0d0/3.0d0)))

  del_alpha = 1.0d0
  del_beta  = 1.0d0
  del_gamma = 1.0d0

  do i = 1, line_each_dim
     alpha = del_alpha*(i-1)
  do j = 1, line_each_dim
     beta  = del_beta*(j-1)
  do k = 1, line_each_dim
     gamma  = del_gamma*(k-1)

     r_start(:) = r_0(:) + alpha*r_0a(:) + beta*r_0b(:) + gamma*r_0c(:)
     idir = 1
     call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)
     idir = -1
     call calbl(r_start(1),r_start(2),r_start(3),isol,idir,l)

   end do
   end do
   end do

   end if

      go to 500

 999  continue
!
      close(10)
      close(20)
      stop
 99   write(6,*) ' Input ERR '
      close(20)
      stop
  contains
!
! ============================================================================
  subroutine reduced_array(a,b)
! ============================================================================
  real*8,dimension(0:NX,0:NY,0:NZ)    :: a
  real*8,dimension(0:nxl,0:nyl,0:nzl) :: b
  integer :: i, j, k
!
  do k = 0, nzl
  do j = 0, nyl
  do i = 0, nxl
   b(i,j,k) = a(nx0+i*leapx,ny0+j*leapy,nz0+k*leapz)
  end do
  end do
  end do
!
  end subroutine reduced_array
  end program bline
