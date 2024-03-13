  module ana
  use common
  implicit none
  contains
! ===============================================================
  subroutine ana_1
! ===============================================================
  real(DP),dimension(0:nx_r,0:ny_r,0:nz_r) :: fx, fy, fz, cb, roi
  real(DP),dimension(0:nx_r,0:ny_r)        :: bh2d
  real(DP),dimension(nz)                   :: decay_index   
  real(DP)                                 :: ene_mag, ene_ken, bz_flux
  real(DP)                                 :: dx_r, dy_r
  real(DP)                                 :: bh2d_total,bh2d_loop, bh2d_max,&
                                              vt_max, zc_mfr
  integer :: i,j,k,index_x,index_y, x_max, y_max, z_max 
  integer :: index_dx, index_dy !! Decay Index

! ===========================================================================
! Geometry set
! ===========================================================================
  dx_r = 1./nx_r
  dy_r = 1./ny_r

! ===========================================================================
! Set the Physical Value
! ===========================================================================
  bt(:,:,:) = sqrt( bx_r(:,:,:)**2 + by_r(:,:,:)**2 + bz_r(:,:,:)**2 )
  vt(:,:,:) = sqrt( vx_r(:,:,:)**2 + vy_r(:,:,:)**2 + vz_r(:,:,:)**2 )
  ct(:,:,:) = sqrt( cx_r(:,:,:)**2 + cy_r(:,:,:)**2 + cz_r(:,:,:)**2 )
  cb(:,:,:) =     ( cx_r(:,:,:)*bx_r(:,:,:)    &
                  + cy_r(:,:,:)*by_r(:,:,:)    &
                  + cz_r(:,:,:)*bz_r(:,:,:))

! ----------------------------------------------------------------------------
! Lorentz Force
! ----------------------------------------------------------------------------
   roi(:,:,:) = 1.0d0/ro_r(:,:,:)
   fx (:,:,:) = roi(:,:,:)*(cy_r(:,:,:)*bz_r(:,:,:) - cz_r(:,:,:)*by_r(:,:,:))
   fy (:,:,:) = roi(:,:,:)*(cz_r(:,:,:)*bx_r(:,:,:) - cx_r(:,:,:)*bz_r(:,:,:))
   fz (:,:,:) = roi(:,:,:)*(cx_r(:,:,:)*by_r(:,:,:) - cy_r(:,:,:)*bx_r(:,:,:))
   cbt(:,:,:) = roi(:,:,:)*(fx(:,:,:)**2 + fy(:,:,:)**2) ! + fz(i,j,k)**2 

   cbt(:,:,0) = 0
   ct (:,:,0) = 0

! -------------------------------------------------------------------------------  
! Height of MFR                  
! ------------------------------------------------------------------------------- 
!  do k = 0,nz_r/2
!    if(bx_r(nx_r/2,ny_r/2,k+1)*bx_r(nx_r/2,ny_r/2,k)<0) then
!      zc_mfr = zc_r(k)
!    endif
!  enddo

!  open(11,file = dir//'BX_1D.'//trim(cloop),form='formatted')
!  do k=  0, nz_r
!     write(11,*) k,bx_r(nx_r/2,ny_r/2,k)
!  end do
!  close(11)

!  open(91,file = dir//'TIME_MFR',form='formatted')
!  write(91,*) atime, zc_mfr

! ------------------------------------------------------------------------------  
! OUTPUT                                                                          
! ------------------------------------------------------------------------------
!  open(11,file = dir//'BXBY.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k = 0, nz_r
!     write(11,*) i,k,bx_r(i,ny_r/2,k),by_r(i,ny_r/2,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'CXCY.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k = 0, nz_r
!     write(11,*) i,k,cx_r(i,ny_r/2,k),cy_r(i,ny_r/2,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'FZ.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k = 0, nz_r
!     write(11,*) i,k,fz(i,ny_r/2,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'DIVV_2D.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,divv_r(i,150,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'CT_BT_1D.'//trim(cloop),form='formatted')
!  do k=  0, nz_r
!     write(11,*) k,ct(120,165,k)/bt(120,165,k),vz_r(120,165,k)
!  end do
!  close(11)

!  open(11,file = dir//'divv_1D.'//trim(cloop),form='formatted')
!  do k=  0, nz_r
!     write(11,*) k,divv_r(150,150,k)
!  end do
!  close(11)

 
! ====================================================================
! Bh
! ====================================================================
!  bh2d(:,:) = sqrt(bx_r(:,:,0)**2 + by_r(:,:,0)**2)
!  bh2d_total = 0.0
!  bh2d_loop  = 0.0
!  do i = 105,195 !0.35<x<0.65
!  do j = 120,210 !0.4 <y<0.70
!     if(bh2d(i,j).ge.0.4) then
!        bh2d_total = bh2d_total + bh2d(i,j)
!        bh2d_loop  = bh2d_loop + 1 
!     end if
!  end do
!  end do
!  bh2d_total = bh2d_total/bh2d_loop


! Maximum Value
!  bh2d_max=0.0
!  do i = 105,195 !0.35<x<0.65
!  do j = 120,210 !0.4 <y<0.70
!     if(bh2d(i,j).ge.bh2d_max) then
!     bh2d_max = bh2d(i,j)
!     end if
!  end do
!  end do

!  open(23,file = dir//'TIME_BH2D',form='formatted')
!  write(23,*)bh2d_total, bh2d_max


! ===================================================================
! Decay Index
! ===================================================================
!  index_dx = 61 ! 59
!  index_dy = 62 ! 60
!  do i =  43,78
!  do j =  50,77

!  if(abs(bz_r(i,j,0)).le.0.01) then
!     decay_index(k) = 0.0d0
!  do k = 1,nz_r-1
!     decay_index(k) = decay_index(k)         &
!                      -                      &
!                     (zc_r(k)/bt(i,j,k))* &
!                      (                      &
!                      ( (bt(i, j, k+1)  & 
!                       - bt(i, j, k-1)) &
!                        *0.5/dz_w(k))/35/27    &
!                      )
!  end do
!  end if

!  end do
!  end do

!  open(12,file = dir//'DECAY_INDEX',form='formatted')
!  do k = 1,nz_r   
!     write(12,*)zc_r(k), decay_index(k)
!  end do   
!  close(12)

! ===================================================================
! Magnetic, Kinetic, Energy 
! ===================================================================
!  ene_mag = 0.0d0
!  do i = nx_r/4.,3.*nx_r/4.
!  do j = ny_r/4.,3.*ny_r/4
!  do k = 0, nz_r*0.2
!     ene_mag = ene_mag + (                      &
!                          bx_r(i,j,k)**2        &
!                       +  by_r(i,j,k)**2        &
!                       +  bz_r(i,j,k)**2        &
!                               )*dx_w(i)*dy_w(j)*dz_w(k)
!!
!  end do
!  end do
!  end do
!  ene_mag = 0.5d0 * ene_mag

!  do i = nx_r/4.,3.*nx_r/4.
!  do j = ny_r/4.,3.*ny_r/4.
!  do k = 1, nz_r*0.2
!     ene_mag = ene_mag + (                                     &
!                       +  bx_r(i,j,k)**2                        &
!                       +  by_r(i,j,k)**2                        &
!                       +  bz_r(i,j,k)**2                        &
!                         )*dx_w(i)*dy_w(j)*dz_w(k)
!  end do
!  end do
!  end do

!     ene_mag = 0.5d0 * ene_mag 

!  open(21,file = dir//'TIME_ENE',form='formatted')
!  write(21,*)ene_mag


! --------------------------------------------------------------------
! Net Flux
! --------------------------------------------------------------------
!  bz_flux = 0.0d0
!  do i = 0,nx_r-1
!  do j = 0,ny_r-1
!         bz_flux = bz_flux + bz_r(i,j,0)*8.0/nx_r/ny_r
!  end do 
!  end do
!  bz_flux = bz_flux/8.0

!  open(11,file=dir//'B_ENE',form='formatted')
!  write(11,*)ene_mag
!  close(11)




! ---------------------------------------------------
! Max of Velocity
! ---------------------------------------------------
!  vt_max = 0.0d0 
!  do i = nx_r/4, 3.*nx_r/4.
!  do j = ny_r/4.,3.*ny_r/4.
!  do k = 1, nz_r*0.2
!     if(vt(i,j,k).gt.vt_max) then
!        vt_max = vt(i,j,k)
!     end if
!  end do
!  end do
!  end do

!  open(23,file = dir//'TIME_VTMAX',form='formatted')
!  write(23,*)vt_max

  end subroutine ana_1
    

  end module ana










