! =================================================================================
  subroutine out(xi,yi,zi,dl,svx,svy,svz,sbx,sby,sbz,scx,scy,scz,run_number,atime)
! =================================================================================
  implicit none
  include "common.f90"
  real*8,intent(inout) :: xi,yi,zi,svx,svy,svz,sbx,sby,sbz,scx,scy,scz,dl,atime
  real*8               :: angle, abs_b0, abs_b, cbz
  integer              :: i,j,k, run_number
  character*2          :: I_NUMBER

! --------------------------------------------------------------------------------
! Deviation of angle from initial state
! --------------------------------------------------------------------------------
!  abs_b0 = sqrt(sbx0**2 + sby0**2 + sbz0**2)
!  abs_b  = sqrt(sbx**2  + sby**2  + sbz**2)

!  angle = (sbx0*sbx + sby0*sby + sbz0*sbz)/abs_b0/abs_b

! --------------------------------------------------------------------------------
!　Lorentz Force
! --------------------------------------------------------------------------------  
  cbz = scx*sby-scy*sbx
  
! --------------------------------------------------------------------------------
! Output
! --------------------------------------------------------------------------------
  write(I_NUMBER,'(i2.2)') run_number
  open(97,file = dir//'TRACE',form='formatted')
  write(97,*)xi,yi,zi

  open(98,file = dir//'TIME_V',form='formatted')
  write(98,*)atime,svz
  
  open(99,file = dir//'TIME_H',form='formatted')
  write(99,*)atime,zi

  open(96,file = dir//'TIME_CBZ',form='formatted')
  write(96,*)atime,cbz
  
!  open(96,file = dir//'TIME_ANGLE',form='formatted')
!  write(96,*)atime,(180/pi)*acos(angle)

30  end subroutine out
