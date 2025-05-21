  module iset
! -----------------------------------------------------------------------------
!                    MODULE iset
! -----------------------------------------------------------------------------
  use common
  use job
  use mhd
  implicit none
  public
  contains

! =============================================================================
  subroutine iset__initial(iset_err)
! =============================================================================
  integer,intent(out) :: iset_err
  integer :: TEST_FILE1 = 10

   if(start_from_initial) then
       call iset__read
       nloop_end = nloop_incmax
   else
       call iset__restart
       nloop_end = nloop + nloop_incmax
   end if

  iset_err = 0
  return

  end subroutine iset__initial
!
!
! ============================================================================
  subroutine iset__read
! ============================================================================
  use common
  real(DP), dimension(-1:NX_W,-1:NY_W,0:NZ_W) :: bx_w, by_w, bz_w
  integer :: iw

! ----------------------------------------------------------------------------
! Magnetic Field
! ----------------------------------------------------------------------------
!  open(11,file = dir_r//"B3D_POTE"//cmyrank, form='unformatted')
!  read(11,err=801) bx, by, bz
!  close(11)

! Non-subdomain version                                                                   
   open(10, file = dir_r//"B3D_POTE_360_360_360_M1",form='unformatted')
   read(10) bx_w, by_w, bz_w
   close(10)

    do i = -1,nx
    do j = -1,ny
    do k =  0,nz
      bx(k,j,i) = bx_w(i+index_x*nx,j+index_y*ny,k)
      by(k,j,i) = by_w(i+index_x*nx,j+index_y*ny,k)
      bz(k,j,i) = bz_w(i+index_x*nx,j+index_y*ny,k)
   end do
   end do
   end do

  iw    = 0
  nloop = 0
  atime = 0
  iwrite = iw

! ----------------------------------------------------------------------------
! Dedner Potential
! ----------------------------------------------------------------------------
  ph(:,:,:) = 0.0
  
! ----------------------------------------------------------------------------
! Initial Density and pressure profile
! ----------------------------------------------------------------------------
!  ro    (:,:,:) = 1.0d0
  ro(:,:,:) = sqrt(bx(:,:,:)**2 + by(:,:,:)**2 + bz(:,:,:)**2) 
  ro_org(:,:,:) = ro(:,:,:)

  pr(:,:,:) = 0.0d0

! ----------------------------------------------------------------------------
! Velocity innrt and Bottom Boundary
! ----------------------------------------------------------------------------
  vx(:,:,:) = 0.0d0 
  vy(:,:,:) = 0.0d0
  vz(:,:,:) = 0.0d0  

!  bz2t_0(:,:) = bz(0,:,:)**2 !! Bz^2 at initial condition
  call iset__velocity

  return

801 write(FILE_SYSOUT,*) ' ### ERR: File Not Found in iset__restart '
    call job__finalize

  end subroutine iset__read
!
! ====================================================================
  subroutine iset__restart
! ====================================================================
  use common
  integer :: iw

! --------------------------------------------------------------------
! read 3d data     
! --------------------------------------------------------------------        
  open(11,file = dir//"B3D_RESTART_030"//cmyrank, form='unformatted')
  read(11,err=801) iw, nloop, atime, dtstep, bx, by, bz, ph
  close(11)

!  open(11,file = dir//"B3D_RESTART2_045"//cmyrank, form='unformatted')   
!  read(11) vx,vy,vz,ro,ro_org
!  close(11)                                                           

  iwrite = iw
  atime0 = atime

! Reset the density
  !  ro    (:,:,:) = 1.0d0
  ro    (:,:,:) = sqrt(bx(:,:,:)**2 + by(:,:,:)**2 + bz(:,:,:)**2)
  ro_org(:,:,:) = ro(:,:,:)

  pr(:,:,:) = 0.0d0

  vx(:,:,:) = 0.0d0
  vy(:,:,:) = 0.0d0
  vz(:,:,:) = 0.0d0
  
  return

801 write(FILE_SYSOUT,*) ' ### ERR: File Not Found in iset__restart '
    call job__finalize

  end subroutine iset__restart

!
! ==============================================================
  subroutine iset__velocity
! ==============================================================
  real(DP), dimension(-1:NY, -1:NX)  :: psi
  real(DP), dimension(-1:NY, -1:NX)  :: vx2d,   vy2d,     &
                                        vx2d_t, vy2d_t,   &
                                        vx2d_c, vy2d_c
  real(DP)                           :: grow_1, bzmax, vmax

! Twist parameter                    
  real(DP), parameter                :: psic       = 0.0003125
  real(DP), parameter                :: gzai       = 1.000
  real(DP), parameter                :: timea      = 1.0
  real(DP), parameter                :: timew      = 0.5
  real(DP), parameter                :: tw         = 0.1d0
! Changed this back to 1 for M1
  real(DP), parameter                :: v0         = 1.0e-02

! Divergence & Convergence
  real(DP), parameter                :: vel_conf_1 = 3.e-04
  real(DP), parameter                :: vel_conf_2 = 4.0e-02
  real(DP), parameter                :: time_cri_c = 3.0
  real(DP)                           :: grow_2

! -------------------------------------------------------------    
! Maximum Value                                               
! -------------------------------------------------------------   
  bzmax = 0.0d0
  do i = -1,nx
  do j = -1,ny
     if(abs(bz(0,j,i)).gt.bzmax) then
        bzmax = abs(bz(0,j,i))
     end if
  end do
  end do
  call mpiut__max(bzmax)

! -------------------------------------------------------------
! Lump Function for Twist and Converging Motion
! -------------------------------------------------------------  
!    grow_1 =  1.0d0
   grow_1 = -0.5*tanh(2.0*((atime - time_cri_t) - timea)/timew) + 0.5
!   grow_2 = -0.5*tanh(2.0*((atime - time_cri_c) - timea)/timew) + 0.5 

! -------------------------------------------------------------
! Velocity Field on Twisting Modtion
! -------------------------------------------------------------
!  psi(:,:)  =  psic*grow_1*bz(0,:,:)**2*exp( ( bz(0,:,:)**2 &
!             - bzmax**2 )/(gzai*bzmax)**2 )

! Amari et al. 2003-2011
  psi(:,:)  =  grow_1*(bz(0,:,:)**2)*exp(-(bz(0,:,:)**2 - bzmax**2)/bzmax**2) 

!  Dale...2015
!  psi(:,:) = (bz(0,:,:)**2/bzmax**2)*exp((bz(0,:,:)**2 -bzmax**2)/(tw*bzmax**2))

  do i = -1, nxm1 ! Restricting region to set twisting motion
     ip = i+1
     im = i-1

  do j = -1, nym1
     jp = j+1
     jm = j-1

!
     if(0.43<=xc(i).and.xc(i)<=0.57.and.0.43<=yc(j).and.yc(j)<=0.57) then
            vx2d_t(j,i) =   -( psi(jp,i)*d1y(j,+1)    &
                            + psi(j, i)*d1y(j, 0)    &
                            + psi(jm,i)*d1y(j,-1) )
                                                   
            vy2d_t(j,i) = ( psi(j,ip)*d1x(i,+1)    &
                            + psi(j,i )*d1x(i, 0)    &
                            + psi(j,im)*d1x(i,-1) )
     else
            vx2d_t(j,i) = 0.0d0
            vy2d_t(j,i) = 0.0d0
     endif
  enddo
  enddo

! ------------------------------------------------------------
! Velocity Field on Divergence Motion
! ------------------------------------------------------------
!  if(atime.lt.time_cri_t) then
!     vx2d_c(:,:) = 0.0d0
!     vy2d_c(:,:) = 0.0d0

!  else 

!  do i = 0, nxm1
!     ip = i+1
!     im = i-1
!
!  do j = 0, nym1
!     jp = j+1
!     jm = j-1
!
! Divergence Motion                                       
!     vx2d_c(j,i) =  - grow_2*vel_conf_1*( bz2t_0(j,ip)*d1x(i,+1)    &
!                              +  bz2t_0(j, i)*d1x(i, 0)             &
!                              +  bz2t_0(j,im)*d1x(i,-1) )
!                                                                
!     vy2d_c(j,i) =  - grow_2*vel_conf_1*( bz2t_0(jp,i)*d1y(j,+1)    &
!                              +  bz2t_0(j, i)*d1y(j, 0)             &
!                              +  bz2t_0(jm,i)*d1y(j,-1) )

! Convergence Motion
!      vx2d_c(j,i) = -grow_2*vel_conf_2*bz(0,j,i)                  
!      vy2d_c(j,i) =  grow_2*vel_conf_2*bz(0,j,i) 
!
!  end do
!  end do

!  end if

! ------------------------------------------------------------
! Total Velocity
! ------------------------------------------------------------

! THIS IS WHAT TO EDIT TO CHANGE STRENGTH OF TWISTING MOTION

  vx2d(:,:) = vx2d_t(:,:)
  vy2d(:,:) = vy2d_t(:,:)

! ------------------------------------------------------------
! Data Exchange for X and Y direction 
! ------------------------------------------------------------
  call mpiut__exchange_2d_x(vx2d)
  call mpiut__exchange_2d_x(vy2d)

  call mpiut__exchange_2d_y(vx2d)
  call mpiut__exchange_2d_y(vy2d)

  vx(0,:,:) = vx2d(:,:)
  vy(0,:,:) = vy2d(:,:)
  vz(0,:,:) = 0.0d0

! ------------------------------------------------------------
! Set the initial amplitude
! ------------------------------------------------------------
  vmax = 0.0d0
  do i = -1,nx
  do j = -1,ny
     if(abs(vy(0,j,i)).gt.vmax) then
        vmax = abs(vy(0,j,i))
     end if
  end do
  end do
  call mpiut__max(vmax)

  vx(0,:,:) = v0*vx(0,:,:)/vmax
  vy(0,:,:) = v0*vy(0,:,:)/vmax


!  if(myrank == root) then
!    open(31,file = dir//'GROW',form = 'formatted')
!    write(31,*) atime,time_cri_t,grow_1
!  endif 

  end subroutine iset__velocity

! --------------------------------------------------------------------------
  subroutine iset__fluxemergence
! --------------------------------------------------------------------------
  real(DP) :: x_r,y_r, z_r, r_a, r_b, r_c, dt, r_top, ts, rxy, xx_r, yy_r, zz_r
  real(DP), parameter :: x_center = 0.5d0                ! Location of X
  real(DP), parameter :: y_center = 0.5d0                ! Location of Y
  real(DP), parameter :: r_cri    = 1.5e-02              ! Radius of EF 
  real(DP), parameter :: B_0      = 0.1d0                ! Strength of EF
  real(DP), parameter :: VE       = 1.0e-02              ! Speed of EF
  real(DP), parameter :: theta    = 180.0*pi/180.0  ! Angle of EF to PIL
  logical             :: bool_emerging

  if (R_CRI - VE*(atime-atime0) > 0.d0) then
       bool_emerging = .true.
  else
       bool_emerging = .false.
  endif

  do j =  0, ny
  do i =  0, nx
     x_r = xc(i) - X_CENTER
     y_r = yc(j) - Y_CENTER
     z_r = max(R_CRI - VE*(atime-atime0), 0.0)

  ! Rotation of coordinate to theat                                                       
     xx_r = x_r*cos(theta) - y_r*sin(theta)
     yy_r = x_r*sin(theta) + y_r*cos(theta)
     zz_r = z_r

    rxy = sqrt(R_CRI**2 - zz_r**2)
         r_b = sqrt(xx_r**2 + zz_r**2)
         r_c = sqrt(xx_r**2 + yy_r**2)
         if (r_c <= rxy) then
         bx(0,j,i)   =  - zz_r/r_b*B_0*cos(theta)
         by(0,j,i)   =    zz_r/r_b*B_0*sin(theta)
         bz(0,j,i)   =  + xx_r/r_b*B_0
 !        if (bool_emerging) then
 !            vz(0,j,i)   =  VE !! <== ?? 
 !        else
 !            vz(0,j,i)   =  0.0d0
 !       endif
     end if
     enddo
     enddo

  end subroutine iset__fluxemergence

! ---------------------------------------------------------------------------------
  end module iset






