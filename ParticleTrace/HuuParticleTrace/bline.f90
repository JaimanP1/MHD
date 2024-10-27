! ============================================================================
  subroutine preprg
! ============================================================================
  implicit none
  include "common.f90"

  common /blxy/ blx(nxl),bly(nyl),blz(nzl) !array values were previously hardcoded with 100, e.g. blx(100)
  real(8)               :: blx, bly, blz, blx0, bly0, blz0
  real(8)               :: sbx, sby, sbz, svx, svy, svz, scx, scy, scz
  real(8)               :: dlx, dly, dlz
  real(8)               :: dt, atime, atime_a
  real(8),dimension(nzl):: x,y,z,dl !was previously hardcoded to an array of size 100
  integer               :: iout, idir, ipcell, isol, nop, l0, l1
  integer               :: icell=1,isci=2
  integer               :: iperix=0, iperiy=0
                 !          >>> iperio = 1, periodic line plotted
  integer               :: i, j, k, loop, run_number
  character*3           :: I_NUMBER

  return

! ===========================================================================
  entry calbl(blx0,bly0,blz0,loop,run_number,atime,atime_a)
! ===========================================================================

  i = loop

  if(i == loop0) then
  blx(i) = blx0
  bly(i) = bly0
  blz(i) = blz0
  dl (i) = 0.0d0
  end if
 
! ---------------------------------------------------------------------------------
! Interpolation
! ---------------------------------------------------------------------------------
  call sint(blx(i),bly(i),blz(i),sbx, sby, sbz, svx, svy, svz, scx, scy, scz, iout)

  if(i == loop0) then
  sbx0 = sbx
  sby0 = sby
  sbz0 = sbz
  end if
  write(6,*) sby0
  
! --------------------------------------------------------------------------
! output 
! --------------------------------------------------------------------------
  call out(blx(i),bly(i),blz(i),dl(i),svx,svy,svz,sbx,sby,sbz,scx,scy,scz, &
           run_number,atime)
! 
! --------------------------------------------------------------------------
! Integration
! --------------------------------------------------------------------------
  dt = atime_a - atime

  i = i+1                  
  dlx = dt*svx
  dly = dt*svy 
  dlz = dt*svz
!
  blx(i) = blx(i-1) + dlx
  bly(i) = bly(i-1) + dly
  blz(i) = blz(i-1) + dlz

  dl(i) = sqrt(dlx**2 + dly**2 + dlz**2)

  return
  end subroutine preprg
