! =====================================================================
  subroutine sint(xi,yi,zi,sbx,sby,sbz,svx,svy,svz,scx,scy,scz,iout)
!! Variable name
!  xi,   yi , z1         :: coordinate field line
!  svx, svy, svz         :: velocity field at selected point.
!  sbx, sby, sbz         :: Magnetic field at selected point.
!  scx, scy, scz    
! =====================================================================
  implicit none
  include "common.f90"
  real(8),intent(in)    :: xi,yi,zi
  real(8),intent(out)   :: sbx,sby,sbz,svx,svy,svz,scx,scy,scz
  integer,intent(inout) :: iout
!
  integer               :: ii, ix0, ix1, iy0, iy1, iz0, iz1,i,j,k
  real(8)               :: x, y, z, dx, dy, dz
  real(8)               :: w000, w001, w010, w100, w011, w101, w110, w111
  real(8)               :: wx0, wy0, wz0, wx1, wy1, wz1
  real(8)               :: celldx, celldy, celldz

  iout = 0
  do i = 0,nxl-1
  if(xc(i).le.xi.and.xc(i+1).ge.xi) then
     celldx = xc(i+1) - xc(i)
  end if
  end do
!
  do j = 0,nyl-1
  if(yc(j).le.yi.and.yc(j+1).ge.yi) then
     celldy = yc(j+1) - yc(j)
  end if
  end do
!
  do k = 0,nzl-1
  if(zc(k).le.zi.and.zc(k+1).ge.zi) then
     celldz = zc(k+1) - zc(k)
  end if
  end do

! --------------------------------------
! check the location inside the domain
! --------------------------------------
  if(xi.lt.xc(0)) then
     iout = -100
     return
  else if(xi.gt.xc(nxl)) then
     iout = 100
     return
  end if
!
  if(yi.lt.yc(0)) then
     iout = -101
     return
  else if(yi.gt.yc(nyl)) then
     iout = 101
     return
  end if
!
  if(zi.lt.zc(0)) then
     iout = -102
     return
  else if(zi.gt.zc(nzl)) then
     iout =  102
     return
  end if
!
  x = xi
  y = yi
  z = zi

! ------------------------------------
! Find the cell in which the point is
! ------------------------------------
! [X]
     do i = 0,nxl-1
     if((xc(i)-x)*(xc(i+1)-x).le.0.0) then
         ix0 = i
         go to 10
     end if
     end do
10   ix1 = ix0+1
     wx1 = (x-xc(ix0))/celldx
     wx0 = (xc(ix1)-x)/celldx


! [Y]
     do j = 0,nyl-1
     if((yc(j)-y)*(yc(j+1)-y).le.0.0) then
         iy0 = j
         go to 11
     end if
     end do
11   iy1 = iy0+1
     wy1 = (y-yc(iy0))/celldy
     wy0 = (yc(iy1)-y)/celldy


! [Z]
     do k = 0,nzl-1
     if((zc(k)-z)*(zc(k+1)-z).le.0.0) then
         iz0 = k
         go to 12
     end if
     end do
12   iz1 = iz0+1
     wz1 = (z-zc(iz0))/celldz
     wz0 = (zc(iz1)-z)/celldz

! -------------------------------------
! interpolation (linear)
! -------------------------------------
  w000 = wx0*wy0*wz0
  w001 = wx0*wy0*wz1
  w010 = wx0*wy1*wz0
  w011 = wx0*wy1*wz1
  w100 = wx1*wy0*wz0
  w101 = wx1*wy0*wz1
  w110 = wx1*wy1*wz0
  w111 = wx1*wy1*wz1

    sbx = w000*bx(ix0,iy0,iz0) &
        +w001*bx(ix0,iy0,iz1)  &
        +w010*bx(ix0,iy1,iz0)  &
        +w011*bx(ix0,iy1,iz1)  &
        +w100*bx(ix1,iy0,iz0)  &
        +w101*bx(ix1,iy0,iz1)  &
        +w110*bx(ix1,iy1,iz0)  &
        +w111*bx(ix1,iy1,iz1)
 !                                                                                                                                                                                  
    sby = w000*by(ix0,iy0,iz0) &
        +w001*by(ix0,iy0,iz1)  &
        +w010*by(ix0,iy1,iz0)  &
        +w011*by(ix0,iy1,iz1)  &
        +w100*by(ix1,iy0,iz0)  &
        +w101*by(ix1,iy0,iz1)  &
        +w110*by(ix1,iy1,iz0)  &
        +w111*by(ix1,iy1,iz1)
!                                                                                                                                                                                  
    sbz = w000*bz(ix0,iy0,iz0) &
        +w001*bz(ix0,iy0,iz1)  &
        +w010*bz(ix0,iy1,iz0)  &
        +w011*bz(ix0,iy1,iz1)  &
        +w100*bz(ix1,iy0,iz0)  &
        +w101*bz(ix1,iy0,iz1)  &
        +w110*bz(ix1,iy1,iz0)  &
        +w111*bz(ix1,iy1,iz1)   
!
    svx = w000*vx(ix0,iy0,iz0)  &
        +w001*vx(ix0,iy0,iz1)  &
        +w010*vx(ix0,iy1,iz0)  &
        +w011*vx(ix0,iy1,iz1)  &
        +w100*vx(ix1,iy0,iz0)  &
        +w101*vx(ix1,iy0,iz1)  &
        +w110*vx(ix1,iy1,iz0)  &
        +w111*vx(ix1,iy1,iz1)
 !
    svy = w000*vy(ix0,iy0,iz0)  &
        +w001*vy(ix0,iy0,iz1)  &
        +w010*vy(ix0,iy1,iz0)  &
        +w011*vy(ix0,iy1,iz1)  &
        +w100*vy(ix1,iy0,iz0)  &
        +w101*vy(ix1,iy0,iz1)  &
        +w110*vy(ix1,iy1,iz0)  &
        +w111*vy(ix1,iy1,iz1)
!
    svz = w000*vz(ix0,iy0,iz0)  &
        +w001*vz(ix0,iy0,iz1)  &
        +w010*vz(ix0,iy1,iz0)  &
        +w011*vz(ix0,iy1,iz1)  &
        +w100*vz(ix1,iy0,iz0)  &
        +w101*vz(ix1,iy0,iz1)  &
        +w110*vz(ix1,iy1,iz0)  &
        +w111*vz(ix1,iy1,iz1)
!
   scx = w000*cx(ix0,iy0,iz0)  &
        +w001*cx(ix0,iy0,iz1)  &
        +w010*cx(ix0,iy1,iz0)  &
        +w011*cx(ix0,iy1,iz1)  &
        +w100*cx(ix1,iy0,iz0)  &
        +w101*cx(ix1,iy0,iz1)  &
        +w110*cx(ix1,iy1,iz0)  &
        +w111*cx(ix1,iy1,iz1)
!
   scy = w000*cy(ix0,iy0,iz0)  &
        +w001*cy(ix0,iy0,iz1)  &
        +w010*cy(ix0,iy1,iz0)  &
        +w011*cy(ix0,iy1,iz1)  &
        +w100*cy(ix1,iy0,iz0)  &
        +w101*cy(ix1,iy0,iz1)  &
        +w110*cy(ix1,iy1,iz0)  &
        +w111*cy(ix1,iy1,iz1)
!
   scz = w000*cz(ix0,iy0,iz0)  &
        +w001*cz(ix0,iy0,iz1)  &
        +w010*cz(ix0,iy1,iz0)  &
        +w011*cz(ix0,iy1,iz1)  &
        +w100*cz(ix1,iy0,iz0)  &
        +w101*cz(ix1,iy0,iz1)  &
        +w110*cz(ix1,iy1,iz0)  &
        +w111*cz(ix1,iy1,iz1)

   
  return
  end subroutine sint



