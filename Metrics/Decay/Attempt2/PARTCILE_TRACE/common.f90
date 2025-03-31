include "paramete.f90"
!
      common /comcq/ &
       c1, c2, c3, c4, cq1, cq2, cs2, cq3, cs3
      real*8 :: c1, c2, c3, c4, cq1, cq2, cs2, cq3, cs3
!
      common /cmbfld/ bx(0:nxl,0:nyl,0:nzl), &
                      by(0:nxl,0:nyl,0:nzl), &
                      bz(0:nxl,0:nyl,0:nzl), &
                      vx(0:nxl,0:nyl,0:nzl), &
                      vy(0:nxl,0:nyl,0:nzl), &
                      vz(0:nxl,0:nyl,0:nzl), &
                      cx(0:nxl,0:nyl,0:nzl), &
                      cy(0:nxl,0:nyl,0:nzl), &
                      cz(0:nxl,0:nyl,0:nzl)

      real*8 bx, by, bz, vx, vy, vz, cx, cy, cz

!
      common /combl/ width,rbl,gbl,bbl,wbl
      real*8 width,rbl,gbl,bbl,wbl
!
      common /cmxyz/ xl,yl,zl,                      &
                     xc(0:nxl),yc(0:nyl),zc(0:nzl), &
                     xc0(0:nx),yc0(0:ny),zc0(0:nz)
      real(8) xl,yl,zl,xc,yc,zc,xc0,yc0,zc0
      real(8) x0_dammy
!
      common /comtotal/ sbx0,sby0,sbz0
      real*8 sbx0,sby0,sbz0
!
      common /cometch/ in_dir, out_dir
      character*12,dimension(10) :: in_dir
      character*2 :: out_dir
!
      namelist /nlist00/ xl, yl, zl




