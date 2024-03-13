!  Files:
!      RUN_NUMBER: if not exist, the calculation starts from the
!                  initial condition.
!                  if exits, read the files of the last run_number 
!                  written in the file RUN_NUMBER
  module pset
! -----------------------------------------------------------
!                    MODULE pset
! -----------------------------------------------------------
  use common
  use mpiut
  implicit none

  integer :: index_x_check, index_y_check
  real(DP) :: char_speed=1.0

  public :: pset__init, &
            pset__dt, &
            pset__namelist, &
            pset__integrate

  private :: pset__x_grid, &
             pset__y_grid, &
             pset__z_grid, &
             pset__mode
  contains

! ===========================================================
  subroutine pset__init
! ===========================================================
  integer :: isend, ierr

! ------------------------------<< sysout file >>
  open(FILE_SYSOUT,file = dir//trim(cfile_sysout)//cmyrank, &
       form='formatted')
! ------------------------------<< run number >>
  
  if(myrank == root) then

    open(FILE_RUN_NUMBER,file=dir//trim(cfile_run_number),form='formatted')
       run_number = 0
       do
          read(FILE_RUN_NUMBER,*,end=101) run_number
       end do
101    run_number = run_number + 1
!
!       write(FILE_RUN_NUMBER,*) run_number   ! wrtie the current run_number 
    close(FILE_RUN_NUMBER)

  end if

! ------------------------------
! send run_number to other nodes
! ------------------------------
  isend = run_number
  call mpi_bcast(isend,1,MPI_INTEGER, &
                 root,MPI_COMM_CART,ierr)
  run_number = isend
  crun_number = '_'//chari2(run_number)

! ================= 
! start from initial 
! ==================
  if(run_number == 1) then
     start_from_initial = .true.
  else
     start_from_initial = .false.
  end if

! --------------------
! open coordinate file 
! --------------------
 open(FILE_COORDINATE_X,file = dir//trim(cfile_coordinate_x)//cmyrank, &
      form='formatted')
 open(FILE_COORDINATE_Y,file = dir//trim(cfile_coordinate_y)//cmyrank, &
      form='formatted')
 open(FILE_COORDINATE_Z,file = dir//trim(cfile_coordinate_z)//cmyrank, &
      form='formatted')

! -----------
! output file 
! -----------
  if(myrank == root) then
  open(FILE_TIME_LIST, file = dir//trim(cfile_time_list)//crun_number, &
       form='formatted')
  end if

 open(FILE_OUTPUT_LIST, &
      file = dir//trim(cfile_output_list)//crun_number//cmyrank, &
      form='formatted')

  write(FILE_OUTPUT_LIST,*) 
  write(FILE_OUTPUT_LIST,*) ' ++ 3D zero-beta nonlinear MHD MODEL ++'
  write(FILE_OUTPUT_LIST,*) ' NX, NY, NZ = ',nx,ny,nz
  write(FILE_OUTPUT_LIST,*) 

  write(FILE_OUTPUT_LIST,nml=nlist00)
  write(FILE_OUTPUT_LIST,nml=nlist01)
  write(FILE_OUTPUT_LIST,nml=nlist02)
  write(FILE_OUTPUT_LIST,nml=nlist03)
  write(FILE_OUTPUT_LIST,nml=nlist04)
  write(FILE_OUTPUT_LIST,nml=nlist05)

! -------------------
! set grid coordinate 
! -------------------
  call pset__x_grid
  call pset__y_grid
  call pset__z_grid

  delta_min = min(minval(dx),minval(dy),minval(dz))
  call mpiut__min(delta_min)

  dt_max = delta_min/char_speed


! ---------------
! set mode number 
! ---------------
  call pset__mode

  end subroutine pset__init

! ===========================================================
  subroutine pset__namelist
! ===========================================================
! ++++ input parameter ++++
!
  open(FILE_NAMELIST,file=dir//'NAMELIST')
      read(FILE_NAMELIST,nml=nlist00)
      read(FILE_NAMELIST,nml=nlist01)
      read(FILE_NAMELIST,nml=nlist02)
      read(FILE_NAMELIST,nml=nlist03)
      read(FILE_NAMELIST,nml=nlist04)
      read(FILE_NAMELIST,nml=nlist05)
  close(FILE_NAMELIST)

! ----------------------
! WARNING for time steps 
! ----------------------
  if(nloop_output .gt. nloop_incmax) then
     write(FILE_SYSOUT,*) ':: WARNING : nloop_output > nloop_incmax ', &
     nloop_output, nloop_incmax
  end if

  end subroutine pset__namelist

! ===========================================================
  subroutine pset__x_grid
! ===========================================================
  real(DP) :: delx, xi0, dxmin, dxmax, dxchg, bunbo
  real(DP) :: xleft ! uniform x coordinate at the left end
  real(DP) :: xright ! uniform x coordinate at the right end
  real(DP) :: xlength ! uniform x length included in this rank
  real(DP), dimension(-1:NX_W) :: xc_w
  real(DP), dimension(-1:NY_W) :: yc_w
  integer :: i
!

  xlength = xl / nproc_x  
  delx = xlength / nx
  xleft = xlength * index_x
  xright= xleft + xlength

  do i = -1, nx
     xi0 = xleft + delx*i
     xc(i) = xi0 + epsx*xl/(PI2)*sin(PI2*xi0/xl)
  end do

! ----------------------------------------------------------   
!  open(10,file = dir_r//"COORDINATE_128.dat",form='unformatted')
!  read(10) xc_w,yc_w
!  close(10)!
!
!  do i = -1,nx
!     xc(i) = xc_w(i+index_x*nx)
!  end do
! ----------------------------------------------------------
!
  do i =-1, nxm1
     dx(i) = xc(i+1) - xc(i)
    ddx(i) = 1.0d0/dx(i)

  if(i.eq.-1) then
     dxmin = dx(i)
     dxmax = dx(i)
  else
     dxmin = min(dx(i), dxmin)
     dxmax = max(dx(i), dxmax)
  end if

  if(i.eq.0) then
     dxchg = abs(dx(i)-dx(i-1))/dx(i)
  else if(i.gt.0) then
     dxchg = max(dxchg, abs(dx(i)-dx(i-1))/dx(i))
  end if

  end do
!
   dx(nx) = dx(nxm1)
  ddx(nx)= ddx(nxm1)
!
  do i = 0, nx
     dx2(i) = dx(i-1) + dx(i)
     ddx2(i) = 1.0d0/dx2(i)
!
     bunbo = 1.0d0/(dx(i-1)*dx(i)*dx2(i))
     d1x(i,-1) = -dx(i)**2*bunbo
     d1x(i, 0) = (dx(i)**2-dx(i-1)**2)*bunbo
     d1x(i, 1) =  dx(i-1)**2*bunbo
!
     d2x(i,-1) = 2*dx (i  )*bunbo
     d2x(i, 0) =-2*dx2(i  )*bunbo
     d2x(i, 1) = 2*dx (i-1)*bunbo
  end do


! --------------------------------------------
! Lateral Boundary for x 
! --------------------------------------------
  do i = -1,nx
  if (xc(i).eq.0) then
! 
     dx2(i) = dx(i)
    ddx2(i) = 1.0d0/dx2(i)
!
     bunbo = 1./(dx(i)*dx(i+1)*dx2(i+1))
     d1x(i,-1) = -dx (i+1)*(dx(i+1)+2*dx(i))*bunbo
     d1x(i, 0) =  dx2(i+1)**2*bunbo
     d1x(i, 1) = -dx (i)**2*bunbo
!
     d2x(i,-1) = 2*dx (i+1)*bunbo
     d2x(i, 0) =-2*dx2(i+1)*bunbo
     d2x(i, 1) = 2*dx (i  )*bunbo
!
   end if
   end do
!
   do i = -1,nx
   if(xc(i).eq.xl) then
!
     dx2(i) = dx(i)
    ddx2(i) = 1.0d0/dx2(i)
!
     bunbo = 1./(dx(i-1)*dx(i-2)*dx2(i-1))
     d1x(i,-1) =  dx (i-1)**2*bunbo
     d1x(i, 0) = -dx2(i-1)**2*bunbo
     d1x(i, 1) =  dx (i-2)*(dx(i-2)+2*dx(i-1))*bunbo
!
     d2x(i,-1) = 2*dx (i-1)*bunbo
     d2x(i, 0) =-2*dx2(i-1)*bunbo
     d2x(i, 1) = 2*dx (i-2)*bunbo 
!!
   end if
   end do


! ------------- print minimum and maximum grid size ---------
      write(FILE_OUTPUT_LIST,*) ':: set_x_grid : ',dxmin, dxmax, dxchg
! -----------------------------------------------------------

! ------------- write coordinate file X ---------------------
  write(FILE_COORDINATE_X,'(e25.16)') xc(-1:NX)
! -----------------------------------------------------------
  end subroutine pset__x_grid

! ===========================================================
  subroutine pset__y_grid
! ===========================================================
  implicit none
  real(DP) :: dely, yj0, dymin, dymax, dychg, bunbo
  real(DP) :: ybottom ! uniform y coordinate at the bottom
  real(DP) :: ytop ! uniform y coordinate at the top
  real(DP) :: ylength ! uniform y length included in this rank
  real(DP), dimension(-1:NX_W) :: xc_w
  real(DP), dimension(-1:NY_W) :: yc_w
  integer :: j
!
  ylength = yl / nproc_y
  dely = ylength / ny
  ybottom = ylength * index_y
  ytop    = ybottom + ylength

  do j = -1, ny
     yj0 = ybottom + dely*j
     yc(j) = yj0 + epsy*yl/(PI2)*sin(PI2*yj0/yl)
  end do
! ----------------------------------------------------------                   
!  open(10,file = dir_r//"COORDINATE_128.dat",form='unformatted')
!  read(10) xc_w,yc_w
!  close(10)
!
!  do j = -1,ny
!     yc(j) = yc_w(j+index_y*ny)
!  end do
! ----------------------------------------------------------
!
  do j =-1, nym1
     dy(j) = yc(j+1) - yc(j)
     ddy(j) = 1.0d0/dy(j)

   if(j.eq.-1) then
      dymin = dy(j)
      dymax = dy(j)
   else
      dymin = min(dy(j), dymin)
      dymax = max(dy(j), dymax)
   end if

   if(j.eq.0) then
      dychg = abs(dy(j)-dy(j-1))/dy(j)
   else if(j.gt.0) then
      dychg = max(dychg, abs(dy(j)-dy(j-1))/dy(j))
   end if

   end do
!
   dy(ny) = dy(nym1)
   ddy(ny)= ddy(nym1)
!
   do j = 0, ny
      dy2(j) = dy(j-1) + dy(j)
     ddy2(j) = 1.0d0/dy2(j)
!
      bunbo = 1.0d0/(dy(j-1)*dy(j)*dy2(j))
      d1y(j,-1) = -dy(j)**2*bunbo
      d1y(j, 0) = (dy(j)**2-dy(j-1)**2)*bunbo
      d1y(j, 1) =  dy(j-1)**2*bunbo
!
      d2y(j,-1) = 2*dy (j  )*bunbo
      d2y(j, 0) =-2*dy2(j  )*bunbo
      d2y(j, 1) = 2*dy (j-1)*bunbo
   end do

! ----------------------------------------------
! Latteral Boundary for y
! ----------------------------------------------

  do j = -1,ny
  if(yc(j).eq.0.0) then
!
     dy2(j) = dy(j)
    ddy2(j) = 1.0d0/dy2(j)
!
     bunbo = 1./(dy(j)*dy(j+1)*dy2(j+1))
     d1y(j,-1) = -dy (j+1)*(dy(j+1)+2*dy(j))*bunbo
     d1y(j, 0) =  dy2(j+1)**2*bunbo
     d1y(j, 1) = -dy (j)**2*bunbo
!
     d2y(j,-1) = 2*dy (j+1)*bunbo
     d2y(j, 0) =-2*dy2(j+1)*bunbo
     d2y(j, 1) = 2*dy (j  )*bunbo
!
   end if
   end do
!
   do j = -1,ny
   if(yc(j).eq.yl) then
!
      dy2(j) = dy(j)
      ddy2(j) = 1.0d0/dy2(j)
!
      bunbo = 1./(dy(j-1)*dy(j-2)*dy2(j-1))
      d1y(j,-1) =  dy (j-1)**2*bunbo
      d1y(j, 0) = -dy2(j-1)**2*bunbo
      d1y(j, 1) =  dy (j-2)*(dy(j-2)+2*dy(j-1))*bunbo
!
      d2y(j,-1) = 2*dy (j-1)*bunbo
      d2y(j, 0) =-2*dy2(j-1)*bunbo
      d2y(j, 1) = 2*dy (j-2)*bunbo
!
    end if
    end do
! ------------- print minimum and maximum grid size ---------
      write(FILE_OUTPUT_LIST,*) ':: set_y_grid : ',dymin, dymax, dychg
! -----------------------------------------------------------

! ------------- write coordinate file Y ---------------------
      write(FILE_COORDINATE_Y,'(e25.16)') yc(-1:NY)
! -----------------------------------------------------------
  end subroutine pset__y_grid

! ===========================================================
  subroutine pset__z_grid
! ===========================================================
  implicit none
  real(DP) :: z1, delz, zca, zcb, zcc, zcd, &
              zk, dzmin, dzmax, dzchg, &
              bunbo
  integer :: k

! ------------------------------------------------
! inhomogeneous grid distribution for Z coordinate
! ------------------------------------------------
!  0.0 < fraction_fine_grid < 1.0
!  0.0 < z_fine_grid < zl
!
  if(   z_fine_grid > 0.0 &
  .and. z_fine_grid < zl  &
  .and. fraction_fine_grid > 0.0 &
  .and. fraction_fine_grid < 1.0 &
    ) then

      z1=z_fine_grid/zl
      delz = z1/fraction_fine_grid
!
      zca = (1.-z1-delz*(1.0-fraction_fine_grid)) &     
           /(1.0 &
            -3*fraction_fine_grid &
            +3*fraction_fine_grid**2 &
            -  fraction_fine_grid**3)
      zcb = -3*zca*fraction_fine_grid
      zcc = delz-3*zca*fraction_fine_grid**2 &
                -2*zcb*fraction_fine_grid
      zcd = 1.0-(zca+zcb+zcc)

      do k = 0, nz
         zk = float(k)/float(nz)
         if(k.le.nz*fraction_fine_grid) then
            zc(k) = delz*zk*zl
         else
            zc(k) =(zca*zk**3+zcb*zk**2+zcc*zk+zcd)*zl
         end if
      end do

   else

! ------------------------------------------------
! homogeneous grid distribution for Z coordinate
! ------------------------------------------------
  delz = zl/nz
  do k = 0, nz
     zc(k) = delz*k
  end do

  end if
!
  zc(0) = 0.0
  zc(nz) = zl
!
!
  do k = 0, nzm1
     dz(k) = zc(k+1) - zc(k)
    ddz(k) = 1.0d0/dz(k)

   if(k.eq.0) then
      dzmin = dz(k)
      dzmax = dz(k)
    else
      dzmin = min(dz(k), dzmin)
      dzmax = max(dz(k), dzmax)
    end if

    if(k.eq.1) then
       dzchg = abs(dz(k)-dz(k-1))/dz(k)
    else if(k.gt.1) then
       dzchg = max(dzchg, abs(dz(k)-dz(k-1))/dz(k))
    end if

   end do

      dz(nz) = dz(nzm1)
      ddz(nz) = ddz(nzm1)
!
      do k = 1, nzm1
         dz2(k) = dz(k-1) + dz(k)
        ddz2(k) = 1.0d0/dz2(k)
!
         bunbo = 1./(dz(k)*dz(k-1)*dz2(k))
         d1z(k,-1) = -dz(k)**2*bunbo
         d1z(k, 0) = (dz(k)**2-dz(k-1)**2)*bunbo
         d1z(k, 1) =  dz(k-1)**2*bunbo
!
         d2z(k,-1) = 2*dz (k  )*bunbo
         d2z(k, 0) =-2*dz2(k  )*bunbo
         d2z(k, 1) = 2*dz (k-1)*bunbo
      end do
!
!
    k = 0
      dz2(k) = dz(k)
     ddz2(k) = 1.0d0/dz2(k)

      bunbo = 1./(dz(k)*dz(k+1)*dz2(k+1))
     d1z(k,-1) = -dz (k+1)*(dz(k+1)+2*dz(k))*bunbo
     d1z(k, 0) =  dz2(k+1)**2*bunbo
     d1z(k, 1) = -dz (k)**2*bunbo
!
     d2z(k,-1) = 2*dz (k+1)*bunbo
     d2z(k, 0) =-2*dz2(k+1)*bunbo
     d2z(k, 1) = 2*dz (k  )*bunbo
!
    k = nz
     dz2(k) = dz(k)
    ddz2(k) = 1.0d0/dz2(k)

     bunbo = 1./(dz(k-1)*dz(k-2)*dz2(k-1))
     d1z(k,-1) =  dz (k-1)**2*bunbo
     d1z(k, 0) = -dz2(k-1)**2*bunbo
     d1z(k, 1) =  dz (k-2)*(dz(k-2)+2*dz(k-1))*bunbo
!
     d2z(k,-1) = 2*dz (k-1)*bunbo
     d2z(k, 0) =-2*dz2(k-1)*bunbo
     d2z(k, 1) = 2*dz (k-2)*bunbo         

! ------------- print minimum and maximum grid size ---------
      write(FILE_OUTPUT_LIST,*) ':: set_z_grid : ',dzmin, dzmax, dzchg
! -----------------------------------------------------------

! ------------- write coordinate file Z ---------------------
      write(FILE_COORDINATE_Z,'(e25.16)') zc(0:NZ)
! -----------------------------------------------------------

  end subroutine pset__z_grid

! ============================================================================
  subroutine pset__dt
! ============================================================================ 
  implicit none
  real(DP) :: dt_diff, dt
  real(DP) :: alfven_w, sound_w, ttx, tty, ttz
  integer  :: i,j,k

  dt = dt_max
  do i = 0, NX-1
  do j = 0, NY-1
  do k = 1, NZ-1

! ---------------------------------------------------------------------------
! Determine the time step by fast mode
! ---------------------------------------------------------------------------
  alfven_w =  sqrt(bx(k,j,i)**2 + by(k,j,i)**2 + bz(k,j,i)**2)/sqrt(ro(k,j,i))
  sound_w  =  sqrt(gamma*pr(k,j,i)/ro(k,j,i))
  ttx = dx(i) / ( abs(vx(k,j,i)) + sqrt(alfven_w + sound_w)**2 )
  tty = dy(j) / ( abs(vy(k,j,i)) + sqrt(alfven_w + sound_w)**2 )
  ttz = dz(k) / ( abs(vz(k,j,i)) + sqrt(alfven_w + sound_w)**2 ) 
!
  if(eta(k,j,i).ne.0) then
     dt_diff = min(dx(i),dy(j),dz(k))**2/eta(k,j,i)
  else
     dt_diff = dt_max
  end if
  dt = min(dt,ttx,tty,ttz,dt_diff)
 

  end do
  end do
  end do

  dtstep = cfl*dt
  call mpiut__min(dtstep) 
!
! ---------------------------------------------------------------------------
! write summary data  
! ---------------------------------------------------------------------------
  if(myrank == root) then
     write(FILE_TIME_LIST,*) atime,ave_ene_mag,ave_ene_kin
  end if
!
  if(myrank == root) then
  open(89, file = dir//"TIME_DIVB"//crun_number,form='formatted')
      write(89,*) atime, ave_divb 
  end if
!
 if(myrank == root) then
  open(90, file = dir//"TIME_STEP"//crun_number,form='formatted')
      write(90,*) atime, nloop, dt, dtstep 
  end if

! ---------------------------------------------------------------------------
! time progresses
! ---------------------------------------------------------------------------
  atime = atime + dtstep    
  end subroutine pset__dt

! ===========================================================================
  subroutine pset__integrate(a,f)
! ===========================================================================
  implicit none
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: a
  real(DP) :: f
  integer :: i, j, k

  do i = 0,nx
  if(xc(i).eq.0.or.xc(i).eq.xl) then
   a(:,:,i) = 0.0d0
  end if
  end do

  do j = 0,ny
  if(yc(j).eq.0.or.yc(j).eq.yl) then
   a(:,j,:) = 0.0d0
  end if
  end do

  do k = 0,nz
  if(zc(k).eq.0.or.zc(k).eq.zl) then
   a(k,:,:) = 0.0d0
  end if
  end do

  f = 0.0d0
  do i =  1, NX-1
  do j =  1, NY-1
  do k =  1, NZ-1
     f = f + (a(k  ,j  ,i)+a(k  ,j  ,i+1) &
           +a(k  ,j+1,i)+a(k  ,j+1,i+1))*0.25D0*dx(i)*dy(j)*dz(k)
  end do
  end do
  end do

  call mpiut__sum(f)

  end subroutine pset__integrate

! ============================================================================
  subroutine pset__mode
! ============================================================================
  implicit none
  integer :: i, j

      do i = 1, NX_W
        if(i.le.NX_W/2+1) then
         md(i) = i-1
        else
         md(i) = i-1 - NX_W
        end if

        akx(i) = PI2*md(i)/xl
      end do
! 
      do j = 1, NY_W
        if(j.le.NY_W/2+1) then
         nd(j) = j-1
        else
         nd(j) = j-1 - NY_W
        end if

        aky(j) = PI2*nd(j)/yl
      end do
!
      do i = 1, NX_W
      do j = 1, NY_W
         akk(i,j) = akx(i)**2 + aky(j)**2
         ak (i,j) = sqrt(akk(i,j))
      end do
      end do
!
  end subroutine pset__mode
  end module pset






