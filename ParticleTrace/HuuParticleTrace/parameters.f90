! x    0    1    2    3    4    5    6    7    8    9   10
! 2**x 1    2    4    8   16   32   64  128  256  512 1024
integer,parameter :: NX    = 360
integer,parameter :: NY    = 360
integer,parameter :: NZ    = 360
integer,parameter :: loop0 = 18
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
 
integer,parameter :: nxyz=(nxl+1)*(nyl+1)*(nzl+1)
!
! grind index 0 means origin for each dimension.
!
real(8),parameter :: pi=3.14159265358979329d0,pi2=6.2831853071795864d0

character(*) ,parameter :: dir   = '/project/cstr/Jaiman/sp25/DecayIndex/'
character(*) ,parameter :: dir_r = '/project/si22/jdp46/sp25/Test1/VAPOR/'
