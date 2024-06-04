! ===========================================================
!              Merging code for divided data by MPI
! ===========================================================
  module constants
! -----------------------------------------------------------
!                      MODULE CONSTANTS
! -----------------------------------------------------------
  implicit none

! === processor number ===
  integer, parameter :: N_PROC  = 144  ! number of the processors >= 1
  integer, parameter :: N_PROC_X= 12
  integer, parameter :: N_PROC_Y= 12

! === grid number in a sub-domain ===
  integer, parameter :: NX = 30
  integer, parameter :: NY = 30
  integer, parameter :: NZ = 360

! === grid number in a whole-domain ===
  integer, parameter :: NX_W=NX*N_PROC_X
  integer, parameter :: NY_W=NY*N_PROC_Y
  integer, parameter :: NZ_W=NZ

! === leaped grid to make the sparse data
  integer, parameter :: LEAPX = 1
  integer, parameter :: LEAPY = 1
  integer, parameter :: LEAPZ = 1

! === grid number for the sparse data ===
  integer, parameter :: NX_R=NX_W/LEAPX
  integer, parameter :: NY_R=NY_W/LEAPY
  integer, parameter :: NZ_R=NZ_W/LEAPZ

! === AVS rearrange grid ===
  integer, parameter :: NX_ARR = 360 !200  
  integer, parameter :: NY_ARR = 360 !200  
  integer, parameter :: NZ_ARR = 360 !200  

  integer, parameter :: NX_avs = 0 !50   !! Start point
  integer, parameter :: NY_avs = 0 !50   !! Start Point

! === sub-grid number ===
  integer, parameter :: nxm1=NX-1
  integer, parameter :: nxp1=NX+1
  integer, parameter :: nym1=NY-1
  integer, parameter :: nyp1=NY+1
  integer, parameter :: nzm1=NZ-1
  integer, parameter :: nzp1=NZ+1

! === precision type ===
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

! === mathematical constants ===
  real(DP), parameter :: PI = 3.1415926535897932384626433_DP
  real(DP), parameter :: PI2 = 2*PI
  complex(DPC), parameter :: IUNIT = (0.0_DP, 1.0_DP)

! === READ & OUTPUT FOLDER ============
  character(*) ,parameter :: dir_r = '/project/wangj/node819/Jaiman/Wulver/sp24/ram_4_24/test4/DATA/'
  character(*) ,parameter :: dir   = '/project/wangj/node819/Jaiman/Wulver/sp24/ram_4_24/test4/VAPOR/Merge/'
  character(*) ,parameter :: dir_a = '/project/wangj/node819/Jaiman/Wulver/sp24/ram_4_24/test4/VAPOR/'

! === I/O file number
  integer, parameter :: FILE_SYSOUT       = 06 ! sysout file
  integer, parameter :: FILE_READN        = 08 ! readN
  integer, parameter :: FILE_RUN_NUMBER   = 09 ! number list calculated so far
  integer, parameter :: FILE_NAMELIST     = 10 ! namelist
  integer, parameter :: FILE_COORDINATE_X = 11 ! coordinate x with '(6e15.5)'
  integer, parameter :: FILE_COORDINATE_Y = 12 ! coordinate y with '(6e15.5)'
  integer, parameter :: FILE_COORDINATE_Z = 13 ! coordinate z with '(6e15.5)'
  integer, parameter :: FILE_COORD_R_X    = 15 ! coordinate x for sparse data
  integer, parameter :: FILE_COORD_R_Y    = 16 ! coordinate y for sparse data
  integer, parameter :: FILE_COORD_R_Z    = 17 ! coordinate z for sparse data
  integer, parameter :: FILE_EIGENMODE    = 21 ! eigenmode perturbation
  integer, parameter :: FILE_OUTPUT_LIST  = 50 ! output list
  integer, parameter :: FILE_TIME_LIST    = 60 ! sequential list (time,loop,etc.)
  integer, parameter :: FILE_2D_FIELD     = 65 ! initial 2d equilibrium
  integer, parameter :: FILE_3D_FIELD     = 70 ! 3d results 
  integer, parameter :: FILE_3D_SPARSE    = 75 ! 3d sparse data
  integer, parameter :: FILE_AVS_FIELD    = 77 ! 3d sparse data
  integer, parameter :: FILE_RESTART      = 80 ! 3d last result for restart
  integer, parameter :: FILE_SLICE_XZ     = 90 ! 2d slice data
  integer, parameter :: FILE_SLICE_XY     = 91 ! 2d slice data

  end module constants

