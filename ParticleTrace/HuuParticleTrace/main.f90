! ****************************************************************************
  program bline
! Particle trace porgram written by S.Inoue on Jan.1,2014
! ****************************************************************************

! ----------------------------------------------------------------------------
! Initial set up
! ----------------------------------------------------------------------------
  include "common.f90"
  integer :: i, j, k,l,numbl, line, dim_lines, line_each_dim

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! filed line starting points are determined as
! a point for dim_lines=0
! equi-interval poins on a line for dim_lines = 1
! grid points on a plane for dim_lines        = 2
! grid points in a cube for dim_lines         = 3
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  character clp*3,filen*11
  real(8) :: r_start(3), x0, y0, z0, atime, atime_a
  real(8) :: dxl, dyl, dzl
  real(8) :: r_0(3),r_a(3),r_b(3),r_c(3)
  real(8) :: r_0a(3),r_0b(3),r_0c(3)
  
  character(3) :: time_loop
  character(3) :: time_loop_a
  integer      :: loop


! =============================================================================
! Read coordinate data input from AVS field data 
! =============================================================================
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

! =============================================================================
! Read PARAM_PLOT
! =============================================================================
  open(10,file='Param_bline_0D',status='OLD')
  read(10,*)dim_lines,numbl,rbl,gbl,bbl,wbl

  if(dim_lines.eq.0) then          !! Point
     read(10,*) r_start(1:3)

  else if(dim_lines.eq.1) then     !! Line
  read(10,*) r_0(1:3), r_0a(1:3)

  else if(dim_lines.eq.2) then     !! plane
  read(10,*) r_0(1:3), r_0a(1:3), r_0b(1:3)

  else if(dim_lines.eq.3) then     !! Cube
  read(10,*) r_0(1:3), r_0a(1:3), r_0b(1:3), r_0c(1:3)

  end if
  close (10)

! ==========================================================================
! Time loop
! ==========================================================================
  do loop = loop0,40 !changed from 34 to 41 as that is number of timesteps
  write(time_loop,'(i3.3)')   loop
  write(time_loop_a,'(i3.3)') loop+1   
  write(6,*) time_loop

! --------------------------------------------------------------------------
! Read 3D field from AVS field data
! --------------------------------------------------------------------------
!  open(11,file=dir_r//'B3D.001.BX.R.'//trim(time_loop),  &
!       form='UNFORMATTED',status='OLD')
!  read(11) bx
!  close(11)
  open(10,file=dir_r//'B3D.001.BX.R.'//trim(time_loop),  &
       form='UNFORMATTED',status='OLD')
  read(10) bx
  close(10)
!          Changed above lines from open(11), read(11), close(10) to (11) instead
  open(10,file=dir_r//'B3D.001.BY.R.'//trim(time_loop),  &
       form='UNFORMATTED',status='OLD')
  read(10) by
  close(10)
!                                                                              
  open(10,file=dir_r//'B3D.001.BZ.R.'//trim(time_loop),  &
       form='UNFORMATTED',status='OLD')
  read(10) bz
  close(10)
!                 
  open(10,file=dir_r//'B3D.001.VX.R.'//trim(time_loop),  &
       form='UNFORMATTED',status='OLD')
  read(10) vx
  close(10)
!
  open(10,file=dir_r//'B3D.001.VY.R.'//trim(time_loop),  &
       form='UNFORMATTED',status='OLD')
  read(10) vy
  close(10)
!
  open(10,file=dir_r//'B3D.001.VZ.R.'//trim(time_loop),  &
       form='UNFORMATTED',status='OLD')
  read(10) vz
  close(10)
!
!  open(10,file=dir_r//'B3D.001.CX.R.'//trim(time_loop),  &
!       form='UNFORMATTED',status='OLD')
!  read(10) cx
!  close(10)
!
!  open(10,file=dir_r//'B3D.001.CY.R.'//trim(time_loop),  &
!       form='UNFORMATTED',status='OLD')
!  read(10) cy
!  close(10)
!
!  open(10,file=dir_r//'B3D.001.CZ.R.'//trim(time_loop),  &
!       form='UNFORMATTED',status='OLD')
!  read(10) cz
!  close(10)
!
  open(10,file=dir_r//'TIME.'//trim(time_loop),  &
       form='FORMATTED',status='OLD')
  read(10,*) atime
  close(10)
  write(10,*) atime

  open(10,file=dir_r//'TIME.'//trim(time_loop_a),  &
       form='FORMATTED',status='OLD')
  read(10,*) atime_a
  close(10)
  
! -----------------------------------------------------------------------------
! Debug ... check the original data
! -----------------------------------------------------------------------------
!  open(11,file = dir//'BOTTOM_V.'//trim(time_loop),form='formatted')
!       do i = 0, nxl
!       do j = 0, nyl
!  write(11,*) xc0(i),yc0(j),vz(i,j,3)
!       end do
!  write(11,*)' '
!       end do
!  close(11)  

! =============================================================================
! Tracing a Particle
! =============================================================================
  call preprg
! -----
! POINT
! -----
  if(dim_lines.eq.0) then
  l = 1
  call calbl(r_start(1),r_start(2),r_start(3),loop,l,atime,atime_a)
! ----
! LINE  
! ----
  else if(dim_lines.eq.1) then
!  do line = 1, numbl
  do l = 1, numbl
     r_start(:) = r_0(:) + (l-1)*r_0a(:)
     call calbl(r_start(1),r_start(2),r_start(3),loop,l,atime,atime_a)
  end do

! -----
! PLANE
! -----
  else if(dim_lines.eq.2) then
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
     call calbl(r_start(1),r_start(2),r_start(3),loop,l,atime,atime_a)
   end do
   end do

! ----
! CUBE 
! ----
  else if(dim_lines.eq.3) then
  line_each_dim = max(2, int(real(numbl)**(1.0d0/3.0d0)))

  del_alpha = 1.0d0
  del_beta  = 1.0d0
  del_gamma = 1.0d0
!
  do i = 1, line_each_dim
     alpha = del_alpha*(i-1)
  do j = 1, line_each_dim
     beta  = del_beta*(j-1)
  do k = 1, line_each_dim
     gamma  = del_gamma*(k-1)
!
     r_start(:) = r_0(:) + alpha*r_0a(:) + beta*r_0b(:) + gamma*r_0c(:)
     call calbl(r_start(1),r_start(2),r_start(3),loop,l,atime,atime_a)
!
  end do
  end do
  end do
!
   end if

   end do !! Time loop

  stop
  end program bline
