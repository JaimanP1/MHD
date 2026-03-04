! =====================================================================
  subroutine preprg
! =====================================================================
  implicit none
  include "common.f90"
  integer,parameter     :: mbl=40000

  common /blxy/ blx(mbl),bly(mbl),blz(mbl)
  real(8)               :: blx, bly, blz
  real(8)               :: blx0,bly0,blz0,blx_start,bly_start,blz_start
  real*8                :: bt,ct,ft
  real(8)               :: qbx, qby, qbz, dbx, dby, dbz, bl
  real*8                :: ff_alpha0
  real(8)               :: sbx, sby, sbz, scz, dcell, bb, dbl,ff_alpha
  real(8),dimension(mbl):: x,y,z
  integer               :: iout, idir, ipcell, isol, nop, l0, l1
  integer               :: icell=1,isci=2
  integer               :: iperix=0, iperiy=0
                 !          >>> iperio = 1, periodic line plotted
  integer               :: i, j, k, run_number
  character*3           :: I_NUMBER

! ---------------------------------------------------------------------
! Cofficent of Runge-Kutta-Gill 
! ---------------------------------------------------------------------
  c1  = 0.5d0
  c2  = 1.0d0 - sqrt(0.5d0)
  c3  = 1.0d0 + sqrt(0.5d0)
  c4  = 1.0d0 / 6.0d0
  cq1 = -2.0d0
  cq2 =  1.0d0 - 3.0d0 * c2
  cs2 =  2.0d0 * c2
  cq3 =  1.0d0 - 3.0d0 * cs3
  cs3 =  2.0d0 * c3

  return

! ====================================================================
  entry calbl(blx0,bly0,blz0,isol,idir,run_number)
! Valiable name
! blx0, bly0, blz0 :: the coordinate of an initial condition
! isol = 1   --> eular, otherwise  --> RKG
! ====================================================================
  if(abs(idir).ne.1) then
     write(6,*) '>>WARNING calbl: idir', idir
     return
  end if
  ipcell = idir*icell

! --------------------------------------------------------------------
! start point(i=1) on the filed line
! --------------------------------------------------------------------
  i = 1
  blx(i) = blx0
  bly(i) = bly0
  blz(i) = blz0

! ---------------------------------------------------------------------
! Initial set along the field line
! ---------------------------------------------------------------------
  alpha_total = 0.0
  bl_total    = 0.0

 1  continue                      !! tracing field line
 
! =====================================================================
  if(isol.eq.1) then !! Eular Method
! =====================================================================
  call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,scz,bt,ct,     &
            ff_alpha,dcell,iout)

! ---------------------------------------------------------------------
! check the boundary and integral number 
! "iout.ne.0" indicate the out of boudndary. See the sint.f90
! mbl is an integral number.
! ---------------------------------------------------------------------
  if(iout.ne.0.or.i.gt.mbl-1) then  
     go to 10
  end if  

! ---------------------------------------------------------------------
! check the polarity
! i=1 indicates "on the surface"
! sbz*idl < 0 indicate "the direction of an line integral"
! ---------------------------------------------------------------------
  if(sbz.eq.0.0) then  !! avoid null point
     sbz=1.e-08
  end if

  if(i.eq.1.and.sbz*idir.lt.0) then
     go to 20
  end if

! ---------------------------------------------------------------------
! Initial set physical value
! ---------------------------------------------------------------------
  if(i.eq.1) then
     ff_alpha0 = ff_alpha
     blx_start = blx(i)
     bly_start = bly(i)
     blz_start = blz(i)
     bz_start  = sbz
     cz_start  = scz
     ct_start  = ct 
  end if
! 
! ---------------------------------------------------------------------
! Integration
! ---------------------------------------------------------------------
  bb  = sqrt(sbx**2+sby**2+sbz**2)
!  if(bb.ne.0) then
     sbx = sbx/bb
     sby = sby/bb
     sbz = sbz/bb
!  else
!     sbx = 1.0e-05 
!     sby = 1.0e-05
!     sbz = 1.0e-05
!  end if
  dbl = dcell/ipcell

!
  i = i+1                         !! Up the integral number
  if(i.gt.mbl) go to 10           !! Trap for the loop end
  qbx = dbl*sbx
  qby = dbl*sby 
  qbz = dbl*sbz
!
  blx(i) = blx(i-1) + qbx
  bly(i) = bly(i-1) + qby
  blz(i) = blz(i-1) + qbz

! Line Length
  bl       = sqrt(qbx**2 + qby**2 + qbz**2)
  bl_total = bl_total + bl
!  write(6,*) qbx,qby,qbz

! Calculating Force Free alpha along the field line --------
! Case-1 Above the bottom surface after a critical height:z_cri
  if(blz(i) >= z_cri) then
     alpha_total = alpha_total + ff_alpha*bl
  end if

! Case-2 At the bottom surface
!  alpha_total = alpha_total + ff_alpha*bl

! =========================================================
  else !! Runge Kutaa Gill
! =========================================================

! ------
! step 1
! ------
  call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,scz,bt,ct, &
            ff_alpha,dcell,iout)

  bb  = sqrt(sbx**2+sby**2+sbz**2)
  sbx = sbx/bb
  sby = sby/bb
  sbz = sbz/bb

!  dbl = dcell/ipcell
  dbl = idir*xl/50
  if(blz(i).lt.0.5) dbl=idir*xl/200
  i = i+1
  if(iout.ne.0.or.i.gt.mbl) go to 10
  qbx = dbl*sbx
  qby = dbl*sby 
  qbz = dbl*sbz
!
    blx(i) = blx(i-1) + c1*qbx
    bly(i) = bly(i-1) + c1*qby
    blz(i) = blz(i-1) + c1*qbz
!
! ------
! step 2
! ------
  call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,scz,bt,ct, &
            ff_alpha,dcell,iout)

  bb  = sqrt(sbx**2+sby**2+sbz**2)
  sbx = sbx/bb
  sby = sby/bb
  sbz = sbz/bb

  dbl = dcell/ipcell
  if(iout.ne.0) go to 10
   dbx = dbl*sbx
   dby = dbl*sby
   dbz = dbl*sbz
!
   blx(i) = blx(i) + c2 * (dbx -qbx)
   bly(i) = bly(i) + c2 * (dby -qby)
   blz(i) = blz(i) + c2 * (dbz -qbz)
!
      qbx = cq2 * qbx + cs2 * dbx
      qby = cq2 * qby + cs2 * dby
      qbz = cq2 * qbz + cs2 * dbz
!
! ------
! step 3
! ------
  call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,scz,bt,ct,  &
            ff_alpha,dcell,iout)

  bb  = sqrt(sbx**2+sby**2+sbz**2)
  sbx = sbx/bb
  sby = sby/bb
  sbz = sbz/bb

  dbl = dcell/ipcell
  if(iout.ne.0) go to 10
!
  dbx = dbl*sbx
  dby = dbl*sby
  dbz = dbl*sbz
!
  blx(i) = blx(i) + c3 * (dbx -qbx)
  bly(i) = bly(i) + c3 * (dby -qby)
  blz(i) = blz(i) + c3 * (dbz -qbz)         
!
     qbx = cq3 * qbx + cs3 * dbx
     qby = cq3 * qby + cs3 * dby
     qbz = cq3 * qbz + cs3 * dbz
!
! ------
! step 4
! ------
  call sint(blx(i),bly(i),blz(i),sbx,sby,sbz,scz,bt,ct, &
            ff_alpha,dcell,iout)

  bb  = sqrt(sbx**2+sby**2+sbz**2)
  sbx = sbx/bb
  sby = sby/bb
  sbz = sbz/bb
!
  dbl = dcell/ipcell
  if(iout.ne.0) go to 10
!
! ---------------------------------------------------------
! Output
! ---------------------------------------------------------
!  write(I_NUMBER,'(i2.2)') run_number
!  open(99,file = dir//'ALPHA.'//I_NUMBER,form='formatted')
!  write(99,*) blz(i),ff_alpha
!
!  open(98,file = dir//'Blines.'//I_NUMBER,form='formatted')
!  write(98,*) blx(i),bly(i),blz(i)

     dbx = dbl*sbx
     dby = dbl*sby
     dbz = dbl*sbz
!
  blx(i) = blx(i) + c4 * (dbx - 2.*qbx)
  bly(i) = bly(i) + c4 * (dby - 2.*qby)
  blz(i) = blz(i) + c4 * (dbz - 2.*qbz)
!
  end if
  go to 1


! ==================================================================
! Analysis 
! using final value after tracing "a" field line.
! ==================================================================
 10  nop = i - 1
     call ana(bz_start,cz_start,ct_start,blx_start,bly_start,blz_start,  &
              blx(i),bly(i),blz(i),sbx,sby,sbz,bt,ct,ff_alpha0,ff_alpha, &
              run_number)

 20  return
     end subroutine preprg
