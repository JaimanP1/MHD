! ===========================================================
!      3D zero-beta nonlinear MHD Model 
!      NL3DpwD_f90_00 by Kanya Kusano (kusano@jamstec.go.jp)
! ===========================================================
  module common
! -----------------------------------------------------------
!                    MODULE COMMON
! -----------------------------------------------------------
  use constants
  implicit none

!!## This case is that raw datas were rearranged (x,y,z) ##
! === define of basic strucutre ===
  type vector_3d_
   real(DP), dimension(0:NZ,-1:NY,-1:NX)      :: x
   real(DP), dimension(0:NZ,-1:NY,-1:NX)      :: y
   real(DP), dimension(0:NZ,-1:NY,-1:NX)      :: z
  end type vector_3d_

! === 2D vector field
  real(DP), dimension( 0:NX, 0:NY)        :: bz2t
  real(DP), dimension( 0:NX, 0:NZ)        :: by2l, by2r
  real(DP), dimension( 0:NY, 0:NZ)        :: bx2b, bx2f

! === 3D vector fields in a sub-domain ===
  real(DP), dimension( 0:NZ, -1:NY, -1:NX)     ::  vx,  vy,  vz
  real(DP), dimension( 0:NZ, -1:NY, -1:NX)     ::  bx,  by,  bz
  real(DP), dimension( 0:NZ, -1:NY, -1:NX)     ::  cx,  cy,  cz
  real(DP), dimension( 0:NZ, -1:NY, -1:NX)     ::  ro,  ph,  cb2
  real(DP), dimension( 0:NZ, -1:NY, -1:NX)     ::  divb,divb2,divv

! === 3D vector field from -1 to NX/NY/NZ ===========
  real(DP), dimension(-1:NX_R, -1:NY_R, 0:NZ_R)::  bx_p1   ,by_p1    ,bz_p1, &
                                                   vx_p1   ,vy_p1    ,vz_p1, &
                                                   ro_p1,   ph_p1

! === 3D vector sparse fields in the whole-domain ===
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  vx_r   ,vy_r    ,vz_r
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  bx_r   ,by_r    ,bz_r
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  cx_r   ,cy_r    ,cz_r
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  ro_r   ,ph_r
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  cb2_r  ,bt,  ct ,ctf
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  divb_r ,divv_r
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  alf    ,alf_ave ,beta, &
                                                   cbt ,c_xrt_r
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  vt
  real,     dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  bx_sp, by_sp, bz_sp
  real,     dimension(0:NX_R, 0:NY_R, 0:NZ_R)  ::  vx_sp, vy_sp, vz_sp

! == 3D reverse
  real(DP), dimension(0:NZ_R, 0:NY_R, 0:NZ_R)  ::  bx_rs, by_rs, bz_rs  

! == 2D
  real(DP), dimension(0:NX_R,0:NY_R)           ::  c_xrt

! === Modefied 3D magnetic field ====================
  real(DP), dimension(0:NX_R, 0:NY_R, 0:NZ_R-10) ::  bx_rr, by_rr, bz_rr

! === 3D rearrange field for AVS ====================
 real(DP),dimension(0:NX_ARR,0:NY_ARR,0:NZ_ARR) :: bx_arr, by_arr, bz_arr
 real(DP),dimension(0:NX_ARR,0:NY_ARR,0:NZ_ARR) :: vx_arr, vy_arr, vz_arr
 real(DP),dimension(0:NX_ARR,0:NY_ARR,0:NZ_ARR) :: ro_arr
 real(DP),dimension(0:NX_ARR,0:NY_ARR,0:NZ_ARR) :: cx_arr, cy_arr, cz_arr
 real(DP),dimension(0:NX_ARR,0:NY_ARR,0:NZ_ARR) :: ct_arr, bt_arr, vt_arr, &
                                                   ct_bt_arr, divv_arr

! === 2D Vector Magnetogram
  real(DP), dimension( 0:NX_R, 0:NY_R) ::  bx2d,    by2d,    bz2d
  real(DP), dimension(-1:NX_R,-1:NY_R) ::  bx2d_p1, by2d_p1, bz2d_p1, &
                                           vx2d_p1, vy2d_p1, vz2d_p1

! === geometrical parameters ===
  real(DP) :: xl                     ! size of the x coordinate
  real(DP) :: yl                     ! size of the y coordinate
  real(DP) :: zl                     ! size of the z coordinate
  real(DP) :: epsx                   ! grid packing parameter for x
  real(DP) :: epsy                   ! grid packing parameter for y
  real(DP) :: z_fine_grid            ! z at the top of zine grid region
  real(DP) :: fraction_fine_grid     ! the rate of the fine grid number for 
                                   ! z < z_fine_grid to the all grid.
  real(DP) :: cfl                    ! CFL coefficient
  real(DP) :: cb2_min                ! Limitter for magnetofrictional method
  real(DP), dimension(-1:nx_r)      :: xc_r
  real(DP), dimension(-1:ny_r)      :: yc_r
  real(DP), dimension( 0:nz_r)      :: zc_r
  real(DP), dimension(-1:NX)        :: xc, dx, ddx, dx2, ddx2
  real(DP), dimension(-1:NY)        :: yc, dy, ddy, dy2, ddy2
  real(DP), dimension( 0:NZ)        :: zc, dz, ddz, dz2, ddz2
  real(DP), dimension(-1:NX,-1:1)   :: d1x, d2x
  real(DP), dimension(-1:NY,-1:1)   :: d1y, d2y
  real(DP), dimension(-1:NZ,-1:1)   :: d1z, d2z
  real(DP), dimension(-1:1)         :: d1y_0, d2y_0
  real(DP), dimension(-1:1)         :: d1y_ny,d2y_ny

! === coordinate in the whole domain ===
  real(DP), dimension(-1:NX_W)      :: xc_w, dx_w, dx2_w
  real(DP), dimension(-1:NY_W)      :: yc_w, dy_w, dy2_w
  real(DP), dimension( 0:NZ_W)      :: zc_w, dz_w, dz2_w

! === cooeficient for the differential operators ===
  real(DP), dimension(-1:NX_W,-1:1) :: d1x_w, d2x_w
  real(DP), dimension(-1:NY_W,-1:1) :: d1y_w, d2y_w
  real(DP), dimension(-1:1)         :: d1y_0_w, d2y_0_w
  real(DP), dimension(-1:1)         :: d1y_ny_w, d2y_ny_w
  real(DP), dimension(-1:NZ_W,-1:1) :: d1z_w, d2z_w

! === diffusion parameters ===
  real(DP) :: eta, eta1, cc0, visc
  real(DP) :: vis, vis_dif
  real(DP) :: grav, gamma
  real(DP) :: c_h,  c_p
  real(DP) :: time_cri_t

! === paramter for time and timestep ===
  real(DP)     :: atime, dtstep, cfr 
  real(DP)     :: max_speed
  integer      :: run_number     ! run number
  character(3) :: crun_number    ! character for run number
  integer      :: nloop          ! current loop number
  integer      :: nloop_output   ! output interval step
  integer      :: nloop_incmax   ! maximum calculation loop in the current task
  integer      :: nloop_end      ! last loop number
  integer      :: iwrite         !  index of the data on disk
  character(3) :: cloop
  character(4) :: cwrite

! === parameter for initial condition ===
integer  :: mmode1, mmode2 ! Fourier mode number (min & max) for the perturbation 
logical  :: start_from_initial 
logical  :: reset_velocity 
logical  :: add_perturbation
real(DP) :: ampp ! amplitude perturbation

! === parameter for the boundary condition
real(DP) :: polarity_change_bt

! === parameter for mpi ===
integer  :: nproc, nproc_x, nproc_y, MPI_COMM_CART
integer  :: myrank, root=0
integer  :: rank_right, rank_left    ! neighbor region for x-coordinate
integer  :: rank_up, rank_down       ! neighbor region for y-coordinate 
integer  :: index_x, index_y         ! coordinate index for each region
character(5) :: cmyrank

! === parameter for the optimized method ===
real(DP) :: dt_initial, dt_inc, dt_dec, mu, nu
real(DP) :: ffl_min ! if ffl < ffl_min then bc_para = bc_para + dbc_para
real(DP) :: bc_para, dbc_para

! === file name ===
character(100) :: cfile_sysout
character(100) :: cfile_run_number
character(100) :: cfile_coordinate_x
character(100) :: cfile_coordinate_y
character(100) :: cfile_coordinate_z
character(100) :: cfile_output_list
character(100) :: cfile_time_list
character(100) :: cfile_2d_field
character(100) :: cfile_3d_field
character(100) :: cfile_restart
character(100) :: cfile_2d_boundary='B2D_boundary_condition'

namelist /nlist00/ xl,   yl, zl
namelist /nlist01/ vis,  vis_dif,      &   
                   grav, gamma,        &
                   c_h,  c_p,          &  
                   cfl,                &
                   time_cri_t         
namelist /nlist02/ nloop_incmax,       &
                   nloop_output
namelist /nlist03/ epsx, epsy,         &
                   z_fine_grid,        &
                   fraction_fine_grid        
namelist /nlist04/ cfile_sysout,       &
                   cfile_run_number,   &
                   cfile_coordinate_x, &
                   cfile_coordinate_y, &
                   cfile_coordinate_z, &
                   cfile_output_list,  &
                   cfile_time_list,    &
                   cfile_2d_field,     &
                   cfile_3d_field,     &
                   cfile_restart
namelist /nlist05/ nproc_x, nproc_y

end module common


