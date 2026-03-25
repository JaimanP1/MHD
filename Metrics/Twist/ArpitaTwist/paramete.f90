! x    0    1    2    3    4    5    6    7    8    9   10
! 2**x 1    2    4    8   16   32   64  128  256  512 1024
! equivalent of constants.f90 ----------------------------
integer,parameter :: NX  = 360        ! NX_W
integer,parameter :: NY  = 360        ! NY_W
integer,parameter :: NZ  = 360        ! NZ_W
integer,parameter :: nx0 = NX*0, nx1 = NX*1.0
integer,parameter :: ny0 = NY*0, ny1 = NY*1.0
integer,parameter :: nz0 = 0,    nz1 = NZ
integer,parameter :: leapx=1
integer,parameter :: leapy=1
integer,parameter :: leapz=1
!.............................................................
!integer,parameter ::                      &
!                     nxl=(nx1-nx0)/leapx, &
!                     nyl=(ny1-ny0)/leapy, &
!                     nzl=(nz1-nz0)/leapz
integer,parameter  :: nxl  = 360
integer,parameter  :: nyl  = 360
integer,parameter  :: nzl  = 360
 
integer,parameter :: nxyz=(nxl+1)*(nyl+1)*(nzl+1)    ! In this case 361*361*361
!
! grind index 0 means origin for each dimension.
!
real(8),parameter :: pi=3.14159265358979329d0,pi2=6.2831853071795864d0

! Adding z_cri to calculate Twist above the bottom boundary
real(8),parameter :: z_cri = 1.0d0 / 360.0d0   ! 1 grid above the bottom boundary
!real(8),parameter :: z_cri = 2.0d0 / 360.0d0    ! 2 grids above the bottom boundary

character(*) ,parameter ::dir_r  = '/project/cstr/Jaiman/sp26/Project1/Test3/VAPOR_Twist/'
character(*) ,parameter ::dir    = '/project/cstr/Jaiman/sp26/Project1/Test3/Metrics/Twist2/'
