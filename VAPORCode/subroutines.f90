  module subroutines
  use common
  contains
! =================================================================
  subroutine boundary(a)
! =================================================================
  implicit none
  real(DP),dimension(0:NX_R, 0:NY_R, 0:NZ_R),intent(inout) :: a
!
   if(mod(NX_W,LEAPX).eq.0) then
      a(NX_R,:,:) = a(0,:,:)
   end if
!
  if(mod(NY_W,LEAPY).eq.0) then
      a(:,NY_R,:) = a(:,0,:)
   end if
  end subroutine boundary
! =================================================================
  subroutine copy_to_sparseglobal(a_s,a_r,p_x,p_y)
! =================================================================
  implicit none
!  node coordinate
     integer, intent(in) :: p_x,p_y
!
!  local data defined in the subdomain
     real(DP),dimension(0:NZ,-1:NY,-1:NX),intent(inout)       :: a_s 
!
!  global data in the sparse array
     real(DP),dimension(0:NX_R, 0:NY_R, 0:NZ_R),intent(inout) :: a_r
!
     integer :: i_min, i_max, j_min, j_max, &
                i, i_g, i_r, j, j_g, j_r, k, k_r

   i_min = LEAPX - mod(NX*(p_x-1), LEAPX)
   if(i_min.eq.LEAPX) i_min =  0

   i_max = NX - mod(NX*(p_x-1)+NX-1, LEAPX)
   if(i_max.eq.NX) i_max = i_max

   j_min = LEAPY - mod(NY*(p_y-1), LEAPY)
   if(j_min.eq.LEAPY) j_min = 0

   j_max = NY - mod(NY*(p_y-1)+NY-1, LEAPY)
   if(j_max.eq.NY) j_max = j_max

     do i = i_min, i_max, LEAPX
        i_g = NX*(p_x-1) + i
        i_r = i_g/LEAPX
        if(mod(i_g,LEAPX).ne.0) write(*,*) '##WARNING i, i_g',i,i_g

     do j = j_min, j_max, LEAPY
        j_g = NY*(p_y-1) + j
        j_r = j_g/LEAPY
        if(mod(j_g,LEAPY).ne.0) write(*,*) '##WARNING j, j_g',j,j_g

     do k = 0, NZ, LEAPZ
        k_r = k/LEAPZ

        a_r(i_r, j_r, k_r) = a_s(k,j,i)

     end do
     end do
     end do

  end subroutine copy_to_sparseglobal
! =================================================================
  subroutine domain_coord(p,p_x,p_y)
! =================================================================
  implicit none
     integer, intent(in) :: p
     integer, intent(out) :: p_x, p_y
!....... p is domain number beginning from 1 not 0!
!....... p_x & p_y are the domain coordinates beginning form 1
!
     p_x = p/N_PROC_Y+1
     p_y = mod(p, N_PROC_Y)+1

  return
  end subroutine domain_coord

! =================================================================
  subroutine pset
! =================================================================
  implicit none
  integer :: i,j,k
  real(DP) :: xi0, delx
  real(DP) :: yj0, dely
  real(DP) :: z1, zca, zcb, zcc, zcd, zk
  real(DP) :: bunbo
  real(DP) :: delz
  real(DP) :: dx_r, dy_r, dz_r

! -----------------------------
! read NAMELIST
! -----------------------------
  open(FILE_NAMELIST,file=dir_r//'NAMELIST')
  read(FILE_NAMELIST,nml=nlist00)
  read(FILE_NAMELIST,nml=nlist01)
  read(FILE_NAMELIST,nml=nlist02)
  read(FILE_NAMELIST,nml=nlist03)
  read(FILE_NAMELIST,nml=nlist04)
  read(FILE_NAMELIST,nml=nlist05)
  close(FILE_NAMELIST)

!### DEBUG
!      write(*,nml=nlist00)
!      write(*,nml=nlist01)
!      write(*,nml=nlist02)
!      write(*,nml=nlist03)
!      write(*,nml=nlist04)
!      write(*,nml=nlist05)

! ------------------------------
! read nonuniform cord
! ------------------------------
!  open(20,file=dir_r//"COORD256_NONUNI_XYZ", form='unformatted')
!  read(20)xc_w, yc_w, zc_w
!  close(20)

! ------------------------------
! set x-coordinate
! ------------------------------
  delx = xl / NX_W
  do i = -1, NX_W
     xi0 = delx*i
     xc_w(i) = xi0 + epsx*xl/(PI2)*sin(PI2*xi0/xl)
  end do

  do i =-1, NX_W-1
     dx_w(i) = xc_w(i+1) - xc_w(i)
  end do
     dx_w(NX_W) = dx_w(NX_W-1)

  do i = 0, NX_W
     dx2_w(i) = dx_w(i-1) + dx_w(i)
!
     bunbo = 1.0d0/(dx_w(i-1)*dx_w(i)*dx2_w(i))
     d1x_w(i,-1) = -dx_w(i)**2*bunbo
     d1x_w(i, 0) = (dx_w(i)**2-dx_w(i-1)**2)*bunbo
     d1x_w(i, 1) =  dx_w(i-1)**2*bunbo
!
     d2x_w(i,-1) = 2*dx_w (i  )*bunbo
     d2x_w(i, 0) =-2*dx2_w(i  )*bunbo
     d2x_w(i, 1) = 2*dx_w (i-1)*bunbo

  end do

!  ------------------------
!  set y-coordinate
!  ------------------------
   dely = yl / NY_W
   do j = -1, NY_W
      yj0 = dely*j
      yc_w(j) =yj0 + epsy*yl/(PI2)*sin(PI2*yj0/yl)
   end do

   do j =-1, NY_W-1
      dy_w(j) = yc_w(j+1) - yc_w(j)
   end do
      dy_w(NY_W) = dy_w(NY_W-1)
!
   do j = 0, NY_W
      dy2_w(j) = dy_w(j-1) + dy_w(j)
!
      bunbo = 1.0d0/(dy_w(j-1)*dy_w(j)*dy2_w(j))
      d1y_w(j,-1) = -dy_w(j)**2*bunbo
      d1y_w(j, 0) = (dy_w(j)**2-dy_w(j-1)**2)*bunbo
      d1y_w(j, 1) =  dy_w(j-1)**2*bunbo
!
      d2y_w(j,-1) = 2*dy_w (j  )*bunbo
      d2y_w(j, 0) =-2*dy2_w(j  )*bunbo
      d2y_w(j, 1) = 2*dy_w (j-1)*bunbo
   end do
!
!
!!      j = 0
!!         dy2_w(j) = dy_w(j)
!!$
!!         bunbo = 1./(dy_w(j)*dy_w(j+1)*dy2_w(j+1))
!!         d1y_0_w(-1) = -dy_w (j+1)*(dy_w(j+1)+2*dy_w(j))*bunbo
!!         d1y_0_w( 0) =  dy2_w(j+1)**2*bunbo
!!         d1y_0_w( 1) = -dy_w (j)**2*bunbo
!!$
!!         d2y_0_w(-1) = 2*dy_w (j+1)*bunbo
!!         d2y_0_w( 0) =-2*dy2_w(j+1)*bunbo
!!         d2y_0_w( 1) = 2*dy_w (j  )*bunbo
!!$
!!     j = NY_W
!!$
!!         dy2_w(j) = dy_w(j)
!!$
!!         bunbo = 1./(dy_w(j-1)*dy_w(j-2)*dy2_w(j-1))
!!         d1y_ny_w(-1) =  dy_w (j-1)**2*bunbo
!!         d1y_ny_w( 0) = -dy2_w(j-1)**2*bunbo
!!         d1y_ny_w( 1) =  dy_w (j-2)*(dy_w(j-2)+2*dy_w(j-1))*bunbo
!!$
!!        d2y_ny_w(-1) = 2*dy_w (j-1)*bunbo
!!         d2y_ny_w( 0) =-2*dy2_w(j-1)*bunbo
!!         d2y_ny_w( 1) = 2*dy_w (j-2)*bunbo         

! ----------------------------------
! set z-coordinate
! ----------------------------------
!! inhomogeneous grid distribution for Z coordinate
!  0.0 < fraction_fine_grid < 1.0
!  0.0 < z_fine_grid < zl
!
! DEBUG===================
  write(*,*) '##DEBUG## z_fine_grid, fraction_fine_grid',    &
                        z_fine_grid,fraction_fine_grid
! DEBUG===================

  if(   z_fine_grid > 0.0 &
  .and. z_fine_grid < zl  &
  .and. fraction_fine_grid > 0.0 &
  .and. fraction_fine_grid < 1.0 &
    ) then

  z1=z_fine_grid/zl
  delz = z1/fraction_fine_grid
  zca = (1.-z1-delz*(1.0-fraction_fine_grid)) &     
       /(1.0 &
        -3*fraction_fine_grid &
        +3*fraction_fine_grid**2 &
        -  fraction_fine_grid**3)
   zcb = -3*zca*fraction_fine_grid
   zcc = delz-3*zca*fraction_fine_grid**2 &
            -2*zcb*fraction_fine_grid
   zcd = 1.0-(zca+zcb+zcc)

   do k = 0, nz
      zk = float(k)/float(nz)
      if(k.le.nz*fraction_fine_grid) then
         zc_w(k) = delz*zk*zl
      else
         zc_w(k) =(zca*zk**3+zcb*zk**2+zcc*zk+zcd)*zl
      end if
   end do

   else


!! homogeneous grid distribution for Z coordinate

   delz = zl/nz
   do k = 0, nz
      zc_w(k) = delz*k
   end do
  end if
!
   zc_w(0) = 0.0
   zc_w(nz) = zl

! --------------------------------------------------
! finite difference
! --------------------------------------------------
      do k = 0, NZ_W-1
         dz_w(k) = zc_w(k+1) - zc_w(k)
      end do

      dz_w(NZ_W) = dz_w(NZ_W-1)
!
      do k = 1, NZ_W-1
         dz2_w(k) = dz_w(k-1) + dz_w(k)
!
         bunbo = 1./(dz_w(k)*dz_w(k-1)*dz2_w(k))
         d1z_w(k,-1) = -dz_w(k)**2*bunbo
         d1z_w(k, 0) = (dz_w(k)**2-dz_w(k-1)**2)*bunbo
         d1z_w(k, 1) =  dz_w(k-1)**2*bunbo
!
         d2z_w(k,-1) = 2*dz_w (k  )*bunbo
         d2z_w(k, 0) =-2*dz2_w(k  )*bunbo
         d2z_w(k, 1) = 2*dz_w (k-1)*bunbo
      end do
!
      k = 0
!
         dz2_w(k) = dz_w(k)

         bunbo = 1./(dz_w(k)*dz_w(k+1)*dz2_w(k+1))
         d1z_w(k,-1) = -dz_w (k+1)*(dz_w(k+1)+2*dz_w(k))*bunbo
         d1z_w(k, 0) =  dz2_w(k+1)**2*bunbo
         d1z_w(k, 1) = -dz_w (k)**2*bunbo
!
         d2z_w(k,-1) = 2*dz_w (k+1)*bunbo
         d2z_w(k, 0) =-2*dz2_w(k+1)*bunbo
         d2z_w(k, 1) = 2*dz_w (k  )*bunbo
!
      k = NZ_W
!
         dz2_w(k) = dz_w(k)

         bunbo = 1./(dz_w(k-1)*dz_w(k-2)*dz2_w(k-1))
         d1z_w(k,-1) =  dz_w (k-1)**2*bunbo
         d1z_w(k, 0) = -dz2_w(k-1)**2*bunbo
         d1z_w(k, 1) =  dz_w (k-2)*(dz_w(k-2)+2*dz_w(k-1))*bunbo
!
         d2z_w(k,-1) = 2*dz_w (k-1)*bunbo
         d2z_w(k, 0) =-2*dz2_w(k-1)*bunbo
         d2z_w(k, 1) = 2*dz_w (k-2)*bunbo         

! -----------------------------------------------------
!         output coordinate file for sparse data
! -----------------------------------------------------
! 100  format(1e14.6)
 100  format(1e25.16)
         open(FILE_COORD_R_X,file=dir//'coord.xgc',form='formatted')
         do i = 0, NX_R
            write(FILE_COORD_R_X,100) xc_w(i*LEAPX) 
         end do
         close(FILE_COORD_R_X)
!
         open(FILE_COORD_R_Y,file=dir//'coord.ygc',form='formatted')
         do i = 0, NY_R
            write(FILE_COORD_R_Y,100) yc_w(i*LEAPY) 
         end do
         close(FILE_COORD_R_Y)
!
         open(FILE_COORD_R_Z,file=dir//'coord.zgc',form='formatted')
         do i = 0, NZ_R
            write(FILE_COORD_R_Z,100) zc_w(i*LEAPZ)
         end do
         close(FILE_COORD_R_Z)


!         open(11,file=dir//'coord.xgc',form='formatted')
!         do i = 0, NX_R
!            write(11,*) xc_w(i*LEAPX) 
!         end do
!         close(11)
!
!         open(11,file=dir//'coord.ygc',form='formatted')
!         do i = 0, NY_R
!            write(11,*) yc_w(i*LEAPY) 
!         end do
!         close(11)
!
!         open(11,file=dir//'coord.zgc',form='formatted')
!         do i = 0, NZ_R
!            write(11,*) zc_w(i*LEAPZ)
!         end do
!         close(11)

! ----------------------------------------
! New coordinate
! ----------------------------------------
  do i = -1,nx_r,LEAPX 
     xc_r(i/LEAPX) = xc_w(i)
  end do

  do j = -1,ny_r,LEAPY
     yc_r(j/LEAPY) = yc_w(j)
  end do

  do k = 0,nz_r, LEAPZ
     zc_r(k/LEAPZ) = zc_w(k)
  end do


  return
  end subroutine pset
!
!
! ================================================================
  subroutine calJ(p_x,p_y)
! ================================================================
  implicit none
  integer, intent(in) :: p_x, p_y
  integer :: ip,im,jp,jm,kp,kpp,km,kmm,ig,jg,kg,i,j,k
!
      do i = 0, NX-1
         ip = i+1
         im = i-1
         ig  = NX*(p_x-1) + i ! GLOBAL INDEX
      do j = 0, NY-1
         jp = j+1
         jm = j-1
         jg  = NY*(p_y-1) + j ! GLOBAL INDEX
      do k = 1, NZ-1
         kp = k+1
         km = k-1
         kg  = k              ! GLOBAL INDEX
!
           cx(k,j,i) = (bz(k ,jp,i )*d1y_w(jg,+1)  &
                       +bz(k, j ,i )*d1y_w(jg, 0)  &
                       +bz(k ,jm,i )*d1y_w(jg,-1)) &
                      -(by(kp,j ,i )*d1z_w(kg,+1)  &
                       +by(k ,j ,i )*d1z_w(kg, 0)  &
                       +by(km,j ,i )*d1z_w(kg,-1))  

          cy(k,j,i) =  (bx(kp,j ,i )*d1z_w(kg,+1)  &
                       +bx(k ,j ,i )*d1z_w(kg, 0)  &
                       +bx(km,j ,i )*d1z_w(kg,-1)) &
                      -(bz(k ,j ,ip)*d1x_w(ig,+1)  &
                       +bz(k ,j ,i )*d1x_w(ig, 0)  &
                       +bz(k ,j ,im)*d1x_w(ig,-1))  

          cz(k,j,i) =  (by(k ,j ,ip)*d1x_w(ig,+1)  &
                       +by(k ,j ,i )*d1x_w(ig, 0)  &
                       +by(k ,j ,im)*d1x_w(ig,-1)) &
                      -(bx(k ,jp,i )*d1y_w(jg,+1)  &
                       +bx(k ,j ,i )*d1y_w(jg, 0)  &
                       +bx(k ,jm,i )*d1y_w(jg,-1))  
!
      end do
      end do
      end do

! ---------------------------------------------------------
! electric current at the top & bottom boundary 
! ---------------------------------------------------------
      do i = 0, NX-1
         ip = i+1
         im = i-1
         ig  = NX*(p_x-1) + i ! GLOBAL INDEX
      do j = 0, NY-1
         jp = j+1
         jm = j-1
         jg  = NY*(p_y-1) + j ! GLOBAL INDEX
!
         k  = 0
         kp = k+1
         kpp= k+2
         kg = k
!
         cx(k,j,i) =   (bz(k ,jp,i )*d1y_w(jg,+1)  &
                       +bz(k , j,i )*d1y_w(jg, 0)  &
                       +bz(k ,jm,i )*d1y_w(jg,-1)) &
                      -(by(kpp,j,i )*d1z_w(kg,+1)  &
                       +by(kp, j,i )*d1z_w(kg, 0)  &
                       +by(k , j,i )*d1z_w(kg,-1))  

         cy(k,j,i) =   (bx(kpp,j,i )*d1z_w(kg,+1)  &
                       +bx(kp,j ,i )*d1z_w(kg, 0)  &
                       +bx(k ,j ,i )*d1z_w(kg,-1)) &
                      -(bz(k ,j ,ip)*d1x_w(ig,+1)  &
                       +bz(k ,j ,i )*d1x_w(ig, 0)  &
                       +bz(k ,j ,im)*d1x_w(ig,-1))   

         cz(k,j,i) =   (by(k ,j ,ip)*d1x_w(ig,+1)  &
                       +by(k ,j ,i )*d1x_w(ig, 0)  &
                       +by(k ,j ,im)*d1x_w(ig,-1)) &
                      -(bx(k ,jp,i )*d1y_w(jg,+1)  &
                       +bx(k ,j ,i )*d1y_w(jg, 0)  &
                       +bx(k ,jm,i )*d1y_w(jg,-1))   
!
         k = NZ
         km = k-1
         kmm= k-2
         kg = k
!
         cx(k,j,i) =   (bz(k ,jp,i )*d1y_w(jg,+1)  &
                       +bz(k, j, i )*d1y_w(jg, 0)  &
                       +bz(k ,jm,i )*d1y_w(jg,-1)) &
                      -(by(k ,j ,i )*d1z_w(kg,+1)  &
                       +by(km,j ,i )*d1z_w(kg, 0)  &
                       +by(kmm,j,i )*d1z_w(kg,-1))  

         cy(k,j,i) =   (bx(k ,j ,i )*d1z_w(kg,+1)  &
                       +bx(km,j ,i )*d1z_w(kg, 0)  &
                       +bx(kmm,j,i )*d1z_w(kg,-1)) &
                      -(bz(k ,j ,ip)*d1x_w(ig,+1)  &
                       +bz(k ,j ,i )*d1x_w(ig, 0)  &
                       +bz(k ,j ,im)*d1x_w(ig,-1))  

         cz(k,j,i) =   (by(k ,j ,ip)*d1x_w(ig,+1)  &
                       +by(k ,j ,i )*d1x_w(ig, 0)  &
                       +by(k ,j ,im)*d1x_w(ig,-1)) &
                      -(bx(k ,jp,i )*d1y_w(jg,+1)  &
                       +bx(k, j ,i )*d1y_w(jg, 0)  &
                       +bx(k ,jm,i )*d1y_w(jg,-1))
!
      end do
      end do
!
  end subroutine calJ
!
!
! =================================================================
  subroutine div__b(p_x, p_y)   !! calculation div_b
! =================================================================
  implicit none
  integer, intent(in) :: p_x, p_y
  integer             :: ip,im,jp,jm,kp,kpp,km,kmm,ig,jg,kg,i,j,k
  real(DP)            :: div2

      do i = 0, NX-1
         ip = i+1
         im = i-1
         ig  = NX*(p_x-1) + i ! GLOBAL INDEX
      do j = 0, NY-1
         jp = j+1
         jm = j-1
         jg  = NY*(p_y-1) + j ! GLOBAL INDEX
      do k = 1, NZ-1
         kp = k+1
         km = k-1
         kg  = k              ! GLOBAL INDEX
!
           divb(k,j,i) =  ( bx(k ,j,ip )*d1x_w(ig,+1)   &
                         +  bx(k, j ,i )*d1x_w(ig, 0)   &
                         +  bx(k ,j,im )*d1x_w(ig,-1) ) &
                         +( by(k ,jp,i )*d1y_w(jg,+1)   &
                         +  by(k ,j ,i )*d1y_w(jg, 0)   &
                         +  by(k ,jm,i )*d1y_w(jg,-1) ) &
                         +( bz(kp,j ,i )*d1z_w(kg,+1)   &
                         +  bz(k ,j ,i )*d1z_w(kg, 0)   &
                         +  bz(km,j ,i )*d1z_w(kg,-1) ) 

           divv(k,j,i) =                               &
!                          ( vx(k ,j,ip )*d1x_w(ig,+1)   &
!                         +  vx(k, j ,i )*d1x_w(ig, 0)   &
!                         +  vx(k ,j,im )*d1x_w(ig,-1) ) &
!                         +( vy(k ,jp,i )*d1y_w(jg,+1)   &
!                         +  vy(k ,j ,i )*d1y_w(jg, 0)   &
!                         +  vy(k ,jm,i )*d1y_w(jg,-1) ) &
                         +( vz(kp,j ,i )*d1z_w(kg,+1)   &
                         +  vz(k ,j ,i )*d1z_w(kg, 0)   &
                         +  vz(km,j ,i )*d1z_w(kg,-1) )

!
       end do
       end do
       end do

       divb(0, :,:) = 0.0d0
       divb(NZ,:,:) = 0.0d0
       divv(0, :,:) = 0.0d0
       divv(NZ,:,:) = 0.0d0


  end subroutine div__b
!
!
! =================================================================
  subroutine integral(a,f)
! =================================================================
  implicit none
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: a
  real(DP)                              :: f
  integer                               :: i, j, k
!------
  f = 0.0d0

  do i = 0, NX-1
  do j = 0, NY-1
  do k = 1, NZ-1
     f = f + (a(k  ,j  ,i)+a(k  ,j  ,i+1) &
             +a(k  ,j+1,i)+a(k  ,j+1,i+1))*0.25D0*dx(i)*dy(j)*dz(k)
  end do
  end do
  end do
  end subroutine integral
!
!
! =================================================================
  subroutine phys_arrange
! =================================================================
  implicit none
  real(DP),dimension(0:nx_r,0:ny_r,0:nz_r) :: cb
  integer :: i,j,k

  cb(:,:,:) =   cx_r(:,:,:)*bx_r(:,:,:)   & 
              + cy_r(:,:,:)*by_r(:,:,:)   &
              + cz_r(:,:,:)*bz_r(:,:,:) 


  do k = 0,1
     ct(:,:,k) = 0.0d0
  end do

! do k = 0,3
! alf(:,:,k) = 0.0d0 
! end do

  do i = NX_avs ,NX_avs + NX_ARR
  do j = NY_avs ,NY_avs + NY_ARR
  do k = 0, NZ_ARR
      ro_arr(i-NX_avs,j-NY_avs,k) =    ro_r(i,j,k)
      bx_arr(i-NX_avs,j-NY_avs,k) =    bx_r(i,j,k)
      by_arr(i-NX_avs,j-NY_avs,k) =    by_r(i,j,k)
      bz_arr(i-NX_avs,j-NY_avs,k) =    bz_r(i,j,k)
      bt_arr(i-NX_avs,j-NY_avs,k) =      bt(i,j,k)
      vx_arr(i-NX_avs,j-NY_avs,k) =    vx_r(i,j,k)
      vy_arr(i-NX_avs,j-NY_avs,k) =    vy_r(i,j,k)
      vz_arr(i-NX_avs,j-NY_avs,k) =    vz_r(i,j,k)
      vt_arr(i-NX_avs,j-NY_avs,k) =      vt(i,j,k)
      cx_arr(i-NX_avs,j-NY_avs,k) =    cx_r(i,j,k)
      cy_arr(i-NX_avs,j-NY_avs,k) =    cy_r(i,j,k)
      cz_arr(i-NX_avs,j-NY_avs,k) =    cz_r(i,j,k)
      ct_arr(i-NX_avs,j-NY_avs,k) =      ct(i,j,k)

      ct_bt_arr(i-NX_avs,j-NY_avs,k) = ct(i,j,k)/bt(i,j,k)
      divv_arr (i-NX_avs,j-NY_avs,k) = divv_r(i,j,k)
 end do
 end do
 end do

! 100  format(1e14.6)
 100  format(1e25.16)
         open(FILE_COORD_R_X,file=dir_a//'coord.xgc',form='formatted')
         do i = NX_avs, NX_avs + NX_ARR
            write(FILE_COORD_R_X,100) xc_w(i*LEAPX) 
         end do
         close(FILE_COORD_R_X)
!
         open(FILE_COORD_R_Y,file=dir_a//'coord.ygc',form='formatted')
         do i = NY_avs, NY_avs + NY_ARR
            write(FILE_COORD_R_Y,100) yc_w(i*LEAPY) 
         end do
         close(FILE_COORD_R_Y)
!
         open(FILE_COORD_R_Z,file=dir_a//'coord.zgc',form='formatted')
         do i = 0, NZ_ARR
            write(FILE_COORD_R_Z,100) zc_w(i*LEAPZ)
         end do
         close(FILE_COORD_R_Z)


  end subroutine phys_arrange

! ===========================================================================
  subroutine make_avs_field(cseries)
! ===========================================================================
  implicit none
  integer, parameter                               :: k_s = 1 ! bottom 
  real(DP),dimension(0:nx_r,0:ny_r,0:nz_r-k_s)     :: bx_avs,  by_avs, bz_avs
  real(DP),dimension(0:nx_r,0:ny_r,0:nz_r-k_s)     :: vx_avs,  vy_avs, vz_avs
  real(DP),dimension(0:nx_r,0:ny_r,0:nz_r-k_s)     :: cx_avs,  cy_avs, cz_avs
  real(DP),dimension(0:NX_ARR, 0:NY_ARR, 0:NZ_ARR)         :: bx_vp, by_vp, bz_vp,   & !changed from real to real(DP)
                                                      vx_vp, vy_vp, vz_vp,   &
                                                      cz_vp, ct_vp, ro_vp,   &
                                                      ct_bt_vp, divv_vp,vt_vp
!
  integer                           :: i,j,k
  character*3, intent(in)           :: cseries
  character*100                     :: filename
  filename = trim(cfile_3d_field)//'+'//cseries//'_R.field'

! ----------------------------------------------------------------------------
! Physical value on AVS
! ----------------------------------------------------------------------------

  cbt(:,:,0)    = 0
  cbt(0,:,:)    = 0
  cbt(nx_r,:,:) = 0
  cbt(:,0,   :) = 0
  cbt(:,ny_r,:) = 0

  ct(:   ,:,0)  = 0
  ct(0   ,:,:)  = 0
  ct(nx_r,:,:)  = 0
  ct(:,0,   :)  = 0
  ct(:,ny_r,:)  = 0

  ctf(:   ,:,0) = 0
  ctf(0   ,:,:) = 0
  ctf(nx_r,:,:) = 0
  ctf(:,0,   :) = 0
  ctf(:,ny_r,:) = 0

! ---------------------------------------------------------------------------
! Rearrange for AVS
! ---------------------------------------------------------------------------
  do i = 0,  nx_r
  do j = 0,  ny_r
  do k = k_s,nz_r
     bx_avs(i,j,k-k_s)  =  bx_r(i,j,k)
     by_avs(i,j,k-k_s)  =  by_r(i,j,k)
     bz_avs(i,j,k-k_s)  =  bz_r(i,j,k)
!
     cx_avs(i,j,k-k_s)  =  cx_r(i,j,k)
     cy_avs(i,j,k-k_s)  =  cy_r(i,j,k)
     cz_avs(i,j,k-k_s)  =  cz_r(i,j,k)
!
     vx_avs(i,j,k-k_s)  =  vx_r(i,j,k)
     vy_avs(i,j,k-k_s)  =  vy_r(i,j,k)
     vz_avs(i,j,k-k_s)  =  vz_r(i,j,k)
  end do
  end do
  end do

! ------------------------------------------------------------------
!                         VAPOR
! ------------------------------------------------------------------
  bt(:,:,:) = sqrt(bx_r(:,:,:)**2 + by_r(:,:,:)**2 + bz_r(:,:,:)**2) 

  ro_vp(:,:,:) = LOG( MAX(real(ro_arr(:,:,:)), 1.0E-30) )  ! plotting log of ro, taking care of negatives with max
  bx_vp(:,:,:) = real(bx_arr(:,:,:))
  by_vp(:,:,:) = real(by_arr(:,:,:))
  bz_vp(:,:,:) = real(bz_arr(:,:,:))
  vx_vp(:,:,:) = real(vx_arr(:,:,:))
  vy_vp(:,:,:) = real(vy_arr(:,:,:))
  vz_vp(:,:,:) = real(vz_arr(:,:,:))
  cz_vp(:,:,:) = real(cz_arr(:,:,:))
  ct_vp(:,:,:) = real(ct_arr  (:,:,:))
  vt_vp(:,:,:) = sqrt(vx_vp(:,:,:)**2 + vy_vp(:,:,:)**2 + vz_vp(:,:,:)**2)
  
  ct_bt_vp(:,:,:) = real(ct_arr(:,:,:)/bt_arr(:,:,:))
!  divv_vp (:,:,:) = real(divv_r(:,:,:))
! ------------------------------------------------------------------
!                   AVS File Output
! ------------------------------------------------------------------ 
  open(FILE_AVS_FIELD,file=dir_a//trim(filename)//'_'//trim(cloop), &
       status='UNKNOWN',form='formatted')

       write(FILE_AVS_FIELD,88) 
 88    format(1H#,' AVS field file')
       write(FILE_AVS_FIELD,*) 'ndim = 3'
       write(FILE_AVS_FIELD,*) 'dim1 = ',NX_ARR+1
       write(FILE_AVS_FIELD,*) 'dim2 = ',NY_ARR+1
       write(FILE_AVS_FIELD,*) 'dim3 = ',NZ_ARR+1
       write(FILE_AVS_FIELD,*) 'nspace = 3'
       write(FILE_AVS_FIELD,*) 'veclen = 6'
       write(FILE_AVS_FIELD,*) 'data = double'
       write(FILE_AVS_FIELD,*) 'field = rectilinear'
       write(FILE_AVS_FIELD,*) 'coord 1 file=coord.xgc filetype=ascii'
       write(FILE_AVS_FIELD,*) 'coord 2 file=coord.ygc filetype=ascii'
       write(FILE_AVS_FIELD,*) 'coord 3 file=coord.zgc filetype=ascii'

       write(FILE_AVS_FIELD,*) 'label=Bx'
       write(FILE_AVS_FIELD,*) 'label=By'
       write(FILE_AVS_FIELD,*) 'label=Bz'
       write(FILE_AVS_FIELD,*) 'label=Ro'
       write(FILE_AVS_FIELD,*) 'label=Vx'
       write(FILE_AVS_FIELD,*) 'label=Vy'
       write(FILE_AVS_FIELD,*) 'label=Vz'
       write(FILE_AVS_FIELD,*) 'label=Vt'
       write(FILE_AVS_FIELD,*) 'label=Jx'
       write(FILE_AVS_FIELD,*) 'label=Jy'
       write(FILE_AVS_FIELD,*) 'label=Jz'
       write(FILE_AVS_FIELD,*) 'label=CB2'
       write(FILE_AVS_FIELD,*) 'label=CT_BT'
       write(FILE_AVS_FIELD,*) 'label=DIVV'

100 format(1e25.16) !to correctly format the data

  filename = trim(cfile_3d_field)//'.'//cseries//'.BX.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 1 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) bx_vp !100 stands for label of format specifier
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.BY.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 2 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) by_vp
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.BZ.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 3 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) bz_vp
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.RO'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 4 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) ro_vp
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.VX.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
      'variable 5 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) vx_vp
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.VY.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 6 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) vy_vp
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.VZ.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 5 file='//trim(filename)//' filetype=unformatted'
   open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) vz_vp
  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.VT.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 5 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) vt_vp
  close(FILE_3D_SPARSE)
!
!  filename = trim(cfile_3d_field)//'.'//cseries//'.CX.R'//'.'//trim(cloop)
!  write(FILE_AVS_FIELD,*) &
!       'variable 9 file='//trim(filename)//' filetype=unformatted'
!  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
!       write(FILE_3D_SPARSE) cx_r
!  close(FILE_3D_SPARSE)
!
!  filename = trim(cfile_3d_field)//'.'//cseries//'.CY.R'//'.'//trim(cloop)
!  write(FILE_AVS_FIELD,*) &
!       'variable 10 file='//trim(filename)//' filetype=unformatted'
!  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
!       write(FILE_3D_SPARSE) cy_r
!  close(FILE_3D_SPARSE)
!
!  filename = trim(cfile_3d_field)//'.'//cseries//'.CZ.R'//'.'//trim(cloop)
!  write(FILE_AVS_FIELD,*) &
!       'variable 11 file='//trim(filename)//' filetype=unformatted'
!  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
!       write(FILE_3D_SPARSE) cz_r
!  close(FILE_3D_SPARSE)
!
  filename = trim(cfile_3d_field)//'.'//cseries//'.CB2.R'//'.'//trim(cloop) 
  write(FILE_AVS_FIELD,*) &
       'variable 6 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) ct_vp
  close(FILE_3D_SPARSE)

  filename = trim(cfile_3d_field)//'.'//cseries//'.CT_BT.R'//'.'//trim(cloop)
  write(FILE_AVS_FIELD,*) &
       'variable 6 file='//trim(filename)//' filetype=unformatted'
  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
       write(FILE_3D_SPARSE) ct_bt_vp
  close(FILE_3D_SPARSE)

!  filename = trim(cfile_3d_field)//'.'//cseries//'.DIVV.R'//'.'//trim(cloop)
!  write(FILE_AVS_FIELD,*) &
!       'variable 6 file='//trim(filename)//' filetype=unformatted'
!  open(FILE_3D_SPARSE,file=dir_a//trim(filename),form='unformatted')
!       write(FILE_3D_SPARSE) divv_vp
!  close(FILE_3D_SPARSE)

!
  close(FILE_AVS_FIELD)
  end subroutine make_avs_field

  end module subroutines
  







