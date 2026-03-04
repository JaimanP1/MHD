include "paramete.f90"
!
      common /comcq/ &
       c1, c2, c3, c4, cq1, cq2, cs2, cq3, cs3
      real*8 :: c1, c2, c3, c4, cq1, cq2, cs2, cq3, cs3
!
      common /cmbfld/ bx(0:nxl,0:nyl,0:nzl), &
                      by(0:nxl,0:nyl,0:nzl), &
                      bz(0:nxl,0:nyl,0:nzl), &
                      cx(0:nxl,0:nyl,0:nzl), &
                      cy(0:nxl,0:nyl,0:nzl), &
                      cz(0:nxl,0:nyl,0:nzl), &
                     alf(0:nxl,0:nyl,0:nzl), &
                     ff_alpha_flux_total(5000)

      real*8 bx, by, bz, cx, cy, cz, alf, ff_alpha_flux_total

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
      common /comtotal/ alpha_total, bl_total, bz_start, cz_start, ct_start
      real*8 alpha_total, bl_total, dbl_total, bz_start, cz_start, ct_start
!
      common /cometch/ in_dir, out_dir
      character*12,dimension(10) :: in_dir
      character*2 :: out_dir
!
      namelist /nlist00/ xl, yl, zl




