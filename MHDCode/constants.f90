  module constants
! -----------------------------------------------------------------------------
!                      MODULE CONSTANTS
! -----------------------------------------------------------------------------
  implicit none

! === grid number ===
  integer, parameter :: NX=30
  integer, parameter :: NY=30
  integer, parameter :: NZ=360

! === grid number for the whole domain ===
  integer, parameter :: NX_W=360
  integer, parameter :: NY_W=360
  integer, parameter :: NZ_W=NZ

! === sub-grid number ===
  integer, parameter :: nxm1=NX-1
  integer, parameter :: nxp1=NX+1
  integer, parameter :: nym1=NY-1
  integer, parameter :: nyp1=NY+1
  integer, parameter :: nzm1=NZ-1
  integer, parameter :: nzp1=NZ+1

! === precision type ===
  integer, parameter :: DP  = kind(1.0d0)
  integer, parameter :: SP  = kind(1.0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

! === mathematical constants ===
  real(DP), parameter :: PI = 3.1415926535897932384626433_DP
  real(DP), parameter :: PI2 = 2*PI
  complex(DPC), parameter :: IUNIT = (0.0_DP, 1.0_DP)

! -------
! I/O
! -------
 character(*),parameter :: dir   = '/project/si22/jdp46/sp25/Test1/MHD/'
 character(*),parameter :: dir_r = '/project/wangj/node819/Jaiman/Wulver/sp24/inDATA/'

! character(*),parameter :: dir_r  = '/research/solarlab_nobackup/inosato/MHD/FLUX_EME/POTE_DATA/320/'
! character(*),parameter :: dir    = '/research/solarlab_nobackup/inosato/MHD/FLUX_EME/MHD/DATA/EMF/'
 

! === I/O file number
  integer, parameter :: FILE_SYSOUT       = 06 ! sysout file
  integer, parameter :: FILE_RUN_NUMBER   = 09 ! number list calculated so far
  integer, parameter :: FILE_NAMELIST     = 10 ! namelist
  integer, parameter :: FILE_COORDINATE_X = 11 ! coordinate x with '(6e15.5)'
  integer, parameter :: FILE_COORDINATE_Y = 12 ! coordinate y with '(6e15.5)'
  integer, parameter :: FILE_COORDINATE_Z = 13 ! coordinate z with '(6e15.5)'
  integer, parameter :: FILE_OUTPUT_LIST  = 50 ! output list
  integer, parameter :: FILE_TIME_LIST    = 60 ! sequential list (time,loop,etc.)
  integer, parameter :: FILE_2D_FIELD     = 65 ! 2d magnetogram data
  integer, parameter :: FILE_2D_BOUNDARY  = 66 ! boundary condition (potential&SFT)
  integer, parameter :: FILE_3D_FIELD     = 70 ! 3d results 
  integer, parameter :: FILE_3D_INITIAL   = 75 ! 3d initial field
  integer, parameter :: FILE_RESTART      = 80 ! 3d last result for restart
  integer, parameter :: FILE_AVS_FIELD    = 85 ! AVS field file
  integer, parameter :: FILE_SLICE_XZ     = 90 ! 2d slice data

  end module constants








