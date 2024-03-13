  module common
! -----------------------------------------------------------------------------------
!                                 MODULE COMMON
! -----------------------------------------------------------------------------------
  use constants
  implicit none

! === define of basic strucutre ============
  type vector_3d_
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: x
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: y
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: z
  end type vector_3d_

! === 2D vector fields =====================
  real(DP), dimension(-1:NX_W,-1:NY_W)  :: bx2d, by2d, bz2d 
  real(DP), dimension(0:NZ,-1:NY)       :: dbx_b, dby_b, dbz_b,  &
                                           dbx_f, dby_f, dbz_f,  &
                                           dvx_b, dvy_b, dvz_b,  &
                                           dvx_f, dvy_f, dvz_f
  real(DP), dimension(0:NZ,-1:NX)      ::  dbx_r, dby_r, dbz_r,  &
                                           dbx_l, dby_l, dbz_l,  &
                                           dvx_r, dvy_r, dvz_r,  &
                                           dvx_l, dvy_l, dvz_l  
  real(DP), dimension(-1:NY,-1:NX)     ::  bz2t_0


! === Fourier mode number and wave vector ===
  real(DP), dimension(NX_W,NY_W)        :: akk, ak
  real(DP), dimension(NX_W)             :: akx
  real(DP), dimension(NY_W)             :: aky
  integer                               :: md(NX_W), nd(NY_W)

! === 3D vector fields ======================
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: bx,   by,     bz  ! vector B
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: vx,   vy,     vz  ! velocity 
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: ro,   ro_org, pr  ! density, pressure
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: ph                ! Dedner Potential
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: cx,   cy,     cz  ! rotB 
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: ex,   ey,     ez  ! electric field   
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dbx,  dby,    dbz ! delta_B
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dvx,  dvy,    dvz ! delta_V
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dpr,  dro,    dph ! delta_RO
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: cbx,  cby,    cbz ! rotBxB
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: divb, divb2       ! divB
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: eta               ! resistivity

! === geometrical parameters ================
  real(DP) :: xl                      ! box size for x
  real(DP) :: yl                      ! box size for y
  real(DP) :: zl                      ! box size for z
  real(DP) :: epsx                    ! grid packing parameter for x
  real(DP) :: epsy                    ! grid packing parameter for y
  real(DP) :: z_fine_grid             ! z at the top of zine grid region
  real(DP) :: fraction_fine_grid      ! the rate of the fine grid number for 
                                      ! z < z_fine_grid to the all grid.
  real(DP) :: cfl                     ! CFL coefficient
  real(DP), dimension(-1:NX)      :: xc, dx, ddx, dx2, ddx2 
  real(DP), dimension(-1:NY)      :: yc, dy, ddy, dy2, ddy2
  real(DP), dimension( 0:NZ)      :: zc, dz, ddz, dz2, ddz2
  real(DP), dimension(-1:NX,-1:1) :: d1x, d2x
  real(DP), dimension(-1:NY,-1:1) :: d1y, d2y
  real(DP), dimension(-1:NZ,-1:1) :: d1z, d2z
  real(DP) :: delta_min

! === paramter for time and timestep ========
  real(DP) :: atime, atime0, ffl
  real(DP) :: dtstep, dt_max
  real(DP) :: vis, vis_dif, grav, gamma 
  real(DP) :: c_h, c_p
  real(DP) :: time_cri_t

  real(DP) :: ave_ene_mag, ave_ene_kin,ave_divb
  integer :: run_number       ! run number
  character(3) :: crun_number ! character for run number
  integer :: nloop            ! current loop number
  integer :: nloop_output     ! output interval step
  integer :: nloop_incmax     ! maximum calculation loop in the current task
  integer :: nloop_end        ! last loop number
  integer :: iwrite           ! index of the data on disk
  character(4) :: cwrite
  logical :: start_from_initial 

! === parameter for mpi =====================
  integer :: nproc, nproc_x, nproc_y, MPI_COMM_CART
  integer :: myrank, root=0
  integer :: rank_right, rank_left ! neighbor region for x-coordinate
  integer :: rank_up, rank_down    ! neighbor region for y-coordinate 
  integer :: index_x, index_y      ! coordinate index for each region
  character(5) :: cmyrank

! === file name =============================
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


  namelist /nlist00/ xl, yl, zl
  namelist /nlist01/ vis, vis_dif,  & ! viscosity, viscosity for density
                     grav,          & ! gravity
                     gamma,         & ! ratio for heat conduction
                     c_h, c_p,      & ! Coefficient of Dedner's equation 
                     cfl,           &
                     time_cri_t       ! critical time for the twisting motion
  namelist /nlist02/ nloop_incmax,  &
                     nloop_output
  namelist /nlist03/ epsx, epsy,    &
                     z_fine_grid,   &
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
  contains
! ------------------------------------------------------------------
  function chari2(in)
! ------------------------------------------------------------------
  character(3) chari2
  integer, INTENT(IN) :: in
! to make character of length 2 corresponding to 
! the integer 'in' with pedding.
!
  if(in.gt.99) then
     write(*,*) ' ++ warning: chari3, in > 999:',in
  end if
  write(chari2,'(I2.2)') in
  return
  end function chari2

! -------------------------------------------------------------------
  function chari3(in)
! -------------------------------------------------------------------
  character(3) chari3
  integer, INTENT(IN) :: in
! to make character of length 3 corresponding to 
! the integer 'in' with pedding.
!
  if(in.gt.999) then
     write(*,*) ' ++ warning: chari3, in > 999:',in
  end if
  write(chari3,'(I3.3)') in
  return
  end function chari3

! --------------------------------------------------------------------
  function chari4(in)
! --------------------------------------------------------------------
  character(4) chari4
  integer, INTENT(IN) :: in
! to make character of length 4 corresponding to 
! the integer 'in' with pedding.
!
  if(in.gt.9999) then
     write(*,*) ' ++ warning: chari3, in > 999:',in
  end if
  write(chari4,'(I4.4)') in
  return
  end function chari4

  end module common
