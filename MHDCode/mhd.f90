  module mhd
! ----------------------------------------------------------------------------
!                    MODULE mhd
! ----------------------------------------------------------------------------
  use common
  use mpiut
  use pset
  implicit none
  integer :: i,j,k,l,kp,km,jp,jm,ip,im
  integer :: ipp,imm,jpp,jmm,kpp,kmm
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dxx_bx, dxy_by, dxz_bz
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dyx_bx, dyy_by, dyz_bz
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: dzx_bx, dzy_by, dzz_bz
  public
  contains

! =============================================================================
  subroutine mhd__main
! =============================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: rovx,      rovy,       rovz 
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: div_rov,   ro_lap,     ro_delta
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: curl_e_x,  curl_e_y,   curl_e_z
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: grad_ph_x, grad_ph_y,  grad_ph_z
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: grad_pr_x, grad_pr_y,  grad_pr_z
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: grad_eta_x,grad_eta_y, grad_eta_z
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: gradx_vx,  gradx_vy,   gradx_vz
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: grady_vx,  grady_vy,   grady_vz
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: gradz_vx,  gradz_vy,   gradz_vz
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: bx_lap,    by_lap,     bz_lap
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: vx_lap,    vy_lap,     vz_lap
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: div_v
!
  real(DP) :: rov_dif
  real(DP) :: roi, cbx, cby, cbz, vx_dxvx, vy_dyvx, vz_dzvx,                 &
              vx_dxvy, vy_dyvy, vz_dzvy, vx_dxvz, vy_dyvz, vz_dzvz,          &
              vis_vx, vis_vy, vis_vz, grad_prx, grad_pry, grad_prz
  real(DP) :: grad_eta_jx, grad_eta_jy, grad_eta_jz, eta_bx, eta_by, eta_bz, &
              rot_ex, rot_ey, rot_ez, g_phix, g_phiy, g_phiz
  real(DP) :: div_prv, gam_pr_divv
  real(DP) :: grav_z
  real(DP) :: v_grad_pr, gam_pr_div_v
  real(DP) :: c_divb, c_phi

! ----------------------------------------------------------------------------
! Preparing to calculate the Dedner equation
! ----------------------------------------------------------------------------
  call mhd__grad(ph, grad_ph_x, grad_ph_y, grad_ph_z)

! ---------------------------------------------------------------------------- 
! Modify the Density; sustained by a minimum value                             
! ----------------------------------------------------------------------------
!  update for rho  
   ro(:,:,:) = sqrt(bx(:,:,:)**2 + by(:,:,:)**2 + bz(:,:,:)**2)

!  Set lower limit for rho   
!  do i = -1,nx
!  do j = -1,ny
! do k =  0,nz
!    2D if(ro(k,j,i).ge.1.e-05) then
!        ro(k,j,i) = ro(k,j,i)
!     else
!        ro(k,j,i) = 1.e-05
!     end if
!  end do
!  end do
!  end do

! ----------------------------------------------------------------------------
! Preparing for calculating a continuous equation
! ----------------------------------------------------------------------------
  rovx(:,:,:) = ro(:,:,:)*vx(:,:,:)
  rovy(:,:,:) = ro(:,:,:)*vy(:,:,:)
  rovz(:,:,:) = ro(:,:,:)*vz(:,:,:)
  call mhd__make_div(rovx,rovy,rovz,div_rov)

  ro_delta(:,:,:) = ro(:,:,:) - ro_org(:,:,:)
  call mhd__make_laplace(ro_delta,ro_delta,ro_delta,ro_lap,ro_lap,ro_lap)

! ----------------------------------------------------------------------------
! Preparing for calculate an equation of motion
! ----------------------------------------------------------------------------
  call mhd__grad(vx, gradx_vx,  grady_vx,  gradz_vx)  
  call mhd__grad(vy, gradx_vy,  grady_vy,  gradz_vy)  
  call mhd__grad(vz, gradx_vz,  grady_vz,  gradz_vz)
  call mhd__grad(pr, grad_pr_x, grad_pr_y, grad_pr_z)
  call mhd__make_laplace(vx,vy,vz,vx_lap,vy_lap,vz_lap)

! ----------------------------------------------------------------------------
! Preparing for Induction equation
! ----------------------------------------------------------------------------
  call mhd__make_curl(ex,ey,ez,curl_e_x,curl_e_y,curl_e_z)
  call mhd__make_laplace(bx,by,bz,bx_lap,by_lap,bz_lap)
  call mhd__grad(eta, grad_eta_x, grad_eta_y, grad_eta_z)

! ---------------------------------------------------------------------------
! Preparing for Pressure equation
! ---------------------------------------------------------------------------
  call mhd__make_div(vx,vy,vz,div_v) 

  do i = -1, NX 
  do j = -1, NY
  do k =  0, NZ

! ----------------------------------------------------------------------------
! Sum of each term for continous equation
! ----------------------------------------------------------------------------
  rov_dif = - div_rov(k,j,i)             &
            + vis_dif*ro_lap(k,j,i)

! ----------------------------------------------------------------------------
! Each term for Equation of Motion
! ----------------------------------------------------------------------------
! Density for Equation of Motion 
  roi = 1.0/ro(k,j,i)

! Inertia term                                                                
  vx_dxvx = - vx(k,j,i)*gradx_vx(k,j,i)
  vy_dyvx = - vy(k,j,i)*grady_vx(k,j,i)
  vz_dzvx = - vz(k,j,i)*gradz_vx(k,j,i)

  vx_dxvy = - vx(k,j,i)*gradx_vy(k,j,i)
  vy_dyvy = - vy(k,j,i)*grady_vy(k,j,i)
  vz_dzvy = - vz(k,j,i)*gradz_vy(k,j,i)

  vx_dxvz = - vx(k,j,i)*gradx_vz(k,j,i)
  vy_dyvz = - vy(k,j,i)*grady_vz(k,j,i)
  vz_dzvz = - vz(k,j,i)*gradz_vz(k,j,i)

! Laplace terms                                       
  vis_vx  = vis*vx_lap(k,j,i)
  vis_vy  = vis*vy_lap(k,j,i)
  vis_vz  = vis*vz_lap(k,j,i)

! Pressure
!  grad_prx = - grad_pr_x(k,j,i)
!  grad_pry = - grad_pr_y(k,j,i)
!  grad_prz = - grad_pr_z(k,j,i)

! Lorents Force                     
  cbx = cy(k,j,i)*bz(k,j,i)-cz(k,j,i)*by(k,j,i)
  cby = cz(k,j,i)*bx(k,j,i)-cx(k,j,i)*bz(k,j,i)
  cbz = cx(k,j,i)*by(k,j,i)-cy(k,j,i)*bx(k,j,i)

! Gravity
!  grav_z = - grav*ro(k,j,i)

! -----------------------------------------------------------------------------
! Each term for an indauction equation
! -----------------------------------------------------------------------------
! Convective Term                               
  rot_ex = - curl_e_x(k,j,i)
  rot_ey = - curl_e_y(k,j,i)
  rot_ez = - curl_e_z(k,j,i)

! Diffusion term        
  grad_eta_jx = - ( grad_eta_y(k,j,i)*cz(k,j,i)  &
                  - grad_eta_z(k,j,i)*cy(k,j,i))

  grad_eta_jy = - ( grad_eta_z(k,j,i)*cx(k,j,i)  &
                  - grad_eta_x(k,j,i)*cz(k,j,i))

  grad_eta_jz = - ( grad_eta_x(k,j,i)*cy(k,j,i)  &
                  - grad_eta_y(k,j,i)*cx(k,j,i))

! Diffusion term                                 
  eta_bx =  eta(k,j,i)*bx_lap(k,j,i)
  eta_by =  eta(k,j,i)*by_lap(k,j,i)
  eta_bz =  eta(k,j,i)*bz_lap(k,j,i)

! Divergence Cleaning                             
  g_phix = - grad_ph_x(k,j,i)
  g_phiy = - grad_ph_y(k,j,i)
  g_phiz = - grad_ph_z(k,j,i)

! ------------------------------------------------------------------------
! Each Term for Pressure Equation
! ------------------------------------------------------------------------
!  v_grad_pr     = - ( vx(k,j,i)*grad_pr_x(k,j,i)  &
!                    + vy(k,j,i)*grad_pr_y(k,j,i)  &
!                    + vz(k,j,i)*grad_pr_z(k,j,i))
!
!  gam_pr_div_v  = - gamma*pr(k,j,i)*div_v(k,j,i)

! ------------------------------------------------------------------------ 
! Each Term for Equation on a Divergence Cleaning        
! ------------------------------------------------------------------------ 
  c_divb = - (c_h**2)*divb  (k,j,i)
  c_phi  = -((c_h**2)/c_p)*ph(k,j,i)

! ========================================================================
! Continous Equation
! ========================================================================
  dro(k,j,i) = 0.0d0
!  dro(k,j,i) = dtstep * (rov_dif)

! ========================================================================
! Equation of Motion
! ========================================================================
  dvx(k,j,i) = dtstep * (                &
                            vx_dxvx      &
                          + vy_dyvx      &
                          + vz_dzvx      &
                          + roi*cbx      &
!                          + roi*grad_prx &
                          + vis_vx       &
                        )
!                                                          
  dvy(k,j,i) = dtstep * (                &
                            vx_dxvy      &
                          + vy_dyvy      &
                          + vz_dzvy      &
                          + roi*cby      &
!                          + roi*grad_pry &
                          + vis_vy       &
                        )
!                                              
  dvz(k,j,i) = dtstep * (                &
                          + vx_dxvz      &
                          + vy_dyvz      &
                          + vz_dzvz      &
                          + roi*cbz      &
!                          + roi*grad_prz &
!                          + grav_z       &
                          + vis_vz       &
                          )

! ============================================================================
! Induction Equation                                                 
! ============================================================================ 
  dbx(k,j,i) = dtstep * (             &
                          rot_ex      &
                        + eta_bx      &
                        + grad_eta_jx &
                        + g_phix      &
                        )

  dby(k,j,i) = dtstep * (             &
                          rot_ey      &
                        + eta_by      &
                        + grad_eta_jy &
                        + g_phiy      &
                        )

  dbz(k,j,i) = dtstep * (             &
                          rot_ez      &
                        + eta_bz      &
                        + grad_eta_jz &
                        + g_phiz      &
                        )

! =======================================================================
! Pressure Equation
! =======================================================================
!  dpr(k,j,i) = dtstep * (  v_grad_pr      &
!                         + gam_pr_div_v   &
!                        )
   dpr(k,j,i) = 0.0d0

! =======================================================================     
! The modefied of Divergence Free                                            
! =======================================================================  
  dph(k,j,i) = dtstep * (  c_divb     &
                         + c_phi      &
                        )


  end do
  end do
  end do

! ----------------------------------------------------------------------------
! Boundary Condition at BOTTOM
! ----------------------------------------------------------------------------
! Magnetic Field on the bottom boundary
!  dbx(0,:,:) = dbx(0,:,:)             
!  dby(0,:,:) = dby(0,:,:)
  dbz(0,:,:) = 0.0d0
  dro(0,:,:) = 0.0d0
  dvx(0,:,:) = 0.0d0
  dvy(0,:,:) = 0.0d0
  dvz(0,:,:) = 0.0d0
  dpr(0,:,:) = 0.0d0
  dph(0,:,:) = dph(1,:,:)

! ----------------------------------------------------------------------------
! Boundary Condition at Top
! ----------------------------------------------------------------------------
  dro(NZ,:,:) = 0.0d0
!  dbx(NZ,:,:) = 0.0d0
!  dby(NZ,:,:) = 0.0d0
  dbz(NZ,:,:) = 0.0d0
  dvx(NZ,:,:) = 0.0d0
  dvy(NZ,:,:) = 0.0d0
  dvz(NZ,:,:) = 0.0d0
  dpr(NZ,:,:) = 0.0d0
  dph(NZ,:,:) = dph(NZ-1,:,:)

! -----------------------------------------------------------------------------
! Memory for side boundaries
! -----------------------------------------------------------------------------
  if(index_x == 0) then
     do j =  -1,ny
     do k =   0,nz
        dbx_b(k,j) = dbx(k,j,0)
        dby_b(k,j) = dby(k,j,0)
        dbz_b(k,j) = dbz(k,j,0)

        dvx_b(k,j) = dvx(k,j,0)
        dvy_b(k,j) = dvy(k,j,0)
        dvz_b(k,j) = dvz(k,j,0)
     end do
     end do
!
  else if(index_x == nproc_x-1) then
     do j =  -1,ny
     do k =   0,nz
        dbx_f(k,j) = dbx(k,j,nx)
        dby_f(k,j) = dby(k,j,nx)
        dbz_f(k,j) = dbz(k,j,nx)

        dvx_f(k,j) = dvx(k,j,nx)
        dvy_f(k,j) = dvy(k,j,nx)
        dvz_f(k,j) = dvz(k,j,nx)
     end do
     end do
  end if

  if(index_y==0) then
     do i =  -1,nx
     do k =   0,nz
        dbx_l(k,i) = dbx(k,0,i)
        dby_l(k,i) = dby(k,0,i)
        dbz_l(k,i) = dbz(k,0,i)

        dvx_l(k,i) = dvx(k,0,i)
        dvy_l(k,i) = dvy(k,0,i)
        dvz_l(k,i) = dvz(k,0,i)
     end do
     end do
!
   else if(index_y == nproc_y-1) then
     do i =  -1,nx
     do k =   0,nz
        dbx_r(k,i) = dbx(k,ny,i)
        dby_r(k,i) = dby(k,ny,i)
        dbz_r(k,i) = dbz(k,ny,i)

        dvx_r(k,i) = dvx(k,ny,i)
        dvy_r(k,i) = dvy(k,ny,i)
        dvz_r(k,i) = dvz(k,ny,i)
 
     end do
     end do
   end if

! =============================================================================
! data exchage for Y direction 
! =============================================================================
  if(nproc_y == 1) then
     call mhd__period_y(dbx)
     call mhd__period_y(dby)
     call mhd__period_y(dbz)
     call mhd__period_y(dvx)
     call mhd__period_y(dvy)
     call mhd__period_y(dvz)
     call mhd__period_y(dro)
     call mhd__period_y(dpr)
     call mhd__period_y(dph)
  else
     call mpiut__exchange_y(dbx,dby,dbz)
     call mpiut__exchange_y(dvx,dvy,dvz)
     call mpiut__exchange_y(dro)
     call mpiut__exchange_y(dpr)
     call mpiut__exchange_y(dph)
     call nlff_bound_y            ! to avoid the periodic B.C. 
  end if

! =============================================================================
! data exchage for X direction 
! =============================================================================
  if(nproc_x == 1) then
     call mhd__period_x(dbx)
     call mhd__period_x(dby)
     call mhd__period_x(dbz)
     call mhd__period_x(dvx)
     call mhd__period_x(dvy)
     call mhd__period_x(dvz)
     call mhd__period_x(dro)
     call mhd__period_x(dpr)
     call mhd__period_x(dph)
  else
     call mpiut__exchange_x(dbx,dby,dbz)
     call mpiut__exchange_x(dvx,dvy,dvz)
     call mpiut__exchange_x(dro)
     call mpiut__exchange_x(dpr)
     call mpiut__exchange_x(dph)
     call nlff_bound_x            ! to avoid periodic boundary B.C. 
  end if

! ============================================================================
! set courner
! ============================================================================
  call nlff_bound_xy

  end subroutine mhd__main


! ============================================================================
  subroutine nlff_bound_x
! ============================================================================
  if(xc(0)==0.0) then
     do j = -1,ny
     do k =  0,nz
        dro(k,j,0) = 0.0d0

        dbx(k,j,0) = 0.0d0
        dby(k,j,0) = dby_b(k,j)
        dbz(k,j,0) = dbz_b(k,j)
!        dby(k,j,0) = 0.0d0
!        dbz(k,j,0) = 0.0d0

        dvx(k,j,0) = 0.0d0
        dvy(k,j,0) = 0.0d0
        dvz(k,j,0) = 0.0d0

        dpr(k,j,0) = 0.0d0

        dph(k,j,0) = dph(k,j,1)
      end do
      end do
!
  
  else if(index_x == nproc_x-1) then
     do j = -1,ny
     do k =  0,nz
        dro(k,j,nx) = 0.0d0

        dbx(k,j,nx) = 0.0d0
        dby(k,j,nx) = dby_f(k,j)
        dbz(k,j,nx) = dbz_f(k,j)
!        dby(k,j,nx) = 0.0d0
!        dbz(k,j,nx) = 0.0d0

        dvx(k,j,nx) = 0.0d0
        dvy(k,j,nx) = 0.0d0
        dvz(k,j,nx) = 0.0d0

        dpr(k,j,nx) = 0.0d0

        dph(k,j,nx) = dph(k,j,NX-1)
      end do
      end do
  end if

  end subroutine nlff_bound_x

! =============================================================================
  subroutine nlff_bound_y
! =============================================================================
  if(yc(0)==0.0) then
     do i = -1,nx
     do k =  0,nz
        dro(k,0,i) = 0.0d0
     
        dbx(k,0,i) = dbx_l(k,i)
        dby(k,0,i) = 0.0d0
        dbz(k,0,i) = dbz_l(k,i)
!        dbx(k,0,i) = 0.0d0
!        dbz(k,0,i) = 0.0d0

        dvx(k,0,i) = 0.0d0
        dvy(k,0,i) = 0.0d0
        dvz(k,0,i) = 0.0d0

        dpr(k,0,i) = 0.0d0

        dph(k,0,i) = dph(k,1,i)
     end do
     end do  
!
  else if(index_y == nproc_y-1) then
     do i = -1,nx
     do k =  0,nz
        dro(k,ny,i) = 0.0d0

        dbx(k,ny,i) = dbx_r(k,i)
        dby(k,ny,i) = 0.0d0
        dbz(k,ny,i) = dbz_r(k,i)
!        dbx(k,ny,i) = 0.0d0 
!        dbz(k,ny,i) = 0.0d0

        dvx(k,ny,i) = 0.0d0
        dvy(k,ny,i) = 0.0d0
        dvz(k,ny,i) = 0.0d0

        dpr(k,ny,i) = 0.0d0

        dph(k,ny,i) = dph(k,ny-1,i)
     end do
     end do
  end if  
!
  end subroutine nlff_bound_y
!
! =============================================================================
  subroutine nlff_bound_xy
! =============================================================================
  do i = -1,nx
  do j = -1,ny
  if(xc(i).eq.0.and.yc(j).eq.0) then
     do k = 0,nz
        dro(k,0,0) = 0.0d0
        dbx(k,0,0) = 0.0d0
        dby(k,0,0) = 0.0d0
        dbz(k,0,0) = 0.0d0
        dvx(k,0,0) = 0.0d0
        dvy(k,0,0) = 0.0d0
        dvz(k,0,0) = 0.0d0
        dpr(k,0,0) = 0.0d0
        dph(k,0,0) = 0.0d0
     end do
  end if
  end do
  end do

  do i = -1,nx
  do j = -1,ny
  if(xc(i).eq.xl.and.yc(j).eq.0) then
     do k = 0,nz
        dro(k,0,nx) = 0.0d0
        dbx(k,0,nx) = 0.0d0
        dby(k,0,nx) = 0.0d0
        dbz(k,0,nx) = 0.0d0
        dvx(k,0,nx) = 0.0d0
        dvy(k,0,nx) = 0.0d0
        dvz(k,0,nx) = 0.0d0
        dpr(k,0,nx) = 0.0d0
        dph(k,0,nx) = 0.0d0
      end do
  end if
  end do
  end do

  do i = -1,nx
  do j = -1,ny
  if(xc(i).eq.0.and.yc(j).eq.yl) then
     do k = 0,nz
        dro(k,ny,0) = 0.0d0
        dbx(k,ny,0) = 0.0d0
        dby(k,ny,0) = 0.0d0
        dbz(k,ny,0) = 0.0d0
        dvx(k,ny,0) = 0.0d0
        dvy(k,ny,0) = 0.0d0
        dvz(k,ny,0) = 0.0d0
        dpr(k,ny,0) = 0.0d0
        dph(k,ny,0) = 0.0d0
      end do
  end if
  end do
  end do

  do i = -1,nx
  do j = -1,ny
  if(xc(i).eq.xl.and.yc(j).eq.yl) then
     do k = 0,nz
        dro(k,ny,nx) = 0.0d0
        dbx(k,ny,nx) = 0.0d0
        dby(k,ny,nx) = 0.0d0
        dbz(k,ny,nx) = 0.0d0
        dvx(k,ny,nx) = 0.0d0
        dvy(k,ny,nx) = 0.0d0
        dvz(k,ny,nx) = 0.0d0
        dpr(k,ny,nx) = 0.0d0
        dph(k,ny,nx) = 0.0d0
      end do
  end if
  end do
  end do


  end subroutine nlff_bound_xy  
!
! =============================================================================
  subroutine mhd__sub
! =============================================================================
  real(DP), parameter                :: timea = 1.0
  real(DP), parameter                :: timew = 0.5
  real(DP), parameter                :: t_cri = 2.5
  real(DP), parameter                :: c_cri = 2.0
  real(DP)                           :: grow
  real(DP)                           :: ct
 

! -----------------------------------------------------------------------------
! div B
! -----------------------------------------------------------------------------
  call mhd__make_div(bx,by,bz,divb)
  divb2(:,:,:) = divb(:,:,:)**2
  call pset__integrate(divb2,ave_divb)

! -----------------------------------------------------------------------------
! curl_B
! -----------------------------------------------------------------------------
  call mhd__make_curl(bx,by,bz,cx,cy,cz)

! -----------------------------------------------------------------------------
! set parameter for resistivity
! -----------------------------------------------------------------------------
!  grow = -0.5*tanh(2.0*((atime - t_cri) - timea)/timew) + 0.5
  grow = 1.0

  do i = -1, NX
  do j = -1, NY
  do k =  0, NZ

! -----------------------------------------------------------------------------
! Resistivity
! -----------------------------------------------------------------------------
!  ct   = sqrt(cx(k,j,i)**2 + cy(k,j,i)**2 + cz(k,j,i)**2)
!  if(ct.gt.c_cri.and.k.gt.3) then
!     eta(k,j,i) = 1.0e-05 + 5.0e-04*grow*((ct - c_cri)/c_cri)**2
!  else 
     eta(k,j,i) = 5.e-05     !1.e-05
!  end if

! Bottom & Top
  eta(0, j,i) = 0.0d0
  eta(nz,j,i) = 0.0d0

! Side for x
  if(index_x == 0.0)  then
     eta(k,j,0)=  0.0d0
  else if(index_x == nproc_x-1) then
    eta(k,j,nx) = 0.0d0
  end if

! Side for y
  if(index_y == 0.0) then
     eta(k,0,i)= 0.0d0
  else if(index_y == nproc_y-1) then
    eta(k,ny,i) = 0.0d0
  end if

! -----------------------------------------------------------------------------
! Electric Field
! -----------------------------------------------------------------------------
     ex(k,j,i) = - (vy(k,j,i)*bz(k,j,i)-vz(k,j,i)*by(k,j,i))
     ey(k,j,i) = - (vz(k,j,i)*bx(k,j,i)-vx(k,j,i)*bz(k,j,i))
     ez(k,j,i) = - (vx(k,j,i)*by(k,j,i)-vy(k,j,i)*bx(k,j,i))

  end do
  end do
  end do

  end subroutine mhd__sub

! =============================================================================   
  subroutine mhd__fluxemergence
! =============================================================================   
  real(DP) :: x_r,y_r, z_r, r_a, r_b, r_c, dt, r_top, ts, rxy, xx_r, yy_r, zz_r
  real(DP), parameter :: x_center = 0.500d0              ! Location of X  
  real(DP), parameter :: y_center = 0.500d0              ! Location of Y   
  real(DP), parameter :: r_cri    = 2.0e-02              ! Radius of EF    
  real(DP), parameter :: B_0      = 0.3d0                ! Strength of EF  
  real(DP), parameter :: VE       = 1.0e-02              ! Speed of EF         
  real(DP), parameter :: theta    = pi*(180.0d0)/180.0d0  ! Angle of EF to PIL 
  logical             :: bool_emerging  

  if (R_CRI - VE*(atime-atime0) > 0.d0) then
       bool_emerging = .true.
  else
       bool_emerging = .false.
  endif

  if (bool_emerging) then 
  do j =  0, ny
  do i =  0, nx
     x_r = xc(i) - X_CENTER
     y_r = yc(j) - Y_CENTER
     z_r = max(R_CRI - VE*(atime-atime0), 0.0)

 ! Rotation of coordinate to theat                                                                   
     xx_r = x_r*cos(theta) - y_r*sin(theta)
     yy_r = x_r*sin(theta) + y_r*cos(theta)
     zz_r = z_r

     rxy = sqrt(R_CRI**2 - zz_r**2)
         r_b = sqrt(xx_r**2 + zz_r**2)
         r_c = sqrt(xx_r**2 + yy_r**2)
     if (r_c <= rxy) then
         bx(0,j,i)   =  - zz_r/r_b*B_0*cos(theta)
         by(0,j,i)   =    zz_r/r_b*B_0*sin(theta)
         bz(0,j,i)   =  + xx_r/r_b*B_0
         vz(0,j,i)   =  VE !! <== ??   
     endif
   enddo
   enddo
   
   else
         vz(0,:,:)  = 0.0d0
   end if

  end subroutine mhd__fluxemergence

     
! =============================================================================
  Subroutine mhd__bc_velocity
! =============================================================================
  real(DP), dimension(-1:NY, -1:NX)  :: psi
  real(DP)                           :: grow_1
  real(DP), parameter                :: timea      = 1.0
  real(DP), parameter                :: timew      = 0.5

  grow_1 = -0.5*tanh(2.0*((atime - time_cri_t) - timea)/timew) + 0.5
  vx(0,:,:) = grow_1*vx(0,:,:)
  vy(0,:,:) = grow_1*vy(0,:,:)

  if(myrank == root) then
    open(31,file = dir//'GROW',form = 'formatted')       
     write(31,*) atime, time_cri_t, grow_1   
  endif

  end subroutine mhd__bc_velocity

! =============================================================================
  subroutine mhd__make_div(ax,ay,az,diva)
! =============================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(in)  :: ax, ay, az
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(out) :: diva
  integer :: i,j,k,ip,im,jp,jm,kp,km
! 
  do i = 0, NX-1
     ip=i+1
     im=i-1
  do j = 0, NY-1
     jp=j+1
     jm=j-1
  do k = 1, NZ-1
     kp=k+1
     km=k-1
!
     diva(k,j,i) =                          &
                  ax(k,j,ip)*d1x(i,+1)      &
                + ax(k,j,i )*d1x(i, 0)      &
                + ax(k,j,im)*d1x(i,-1)      &
                + ay(k,jp,i)*d1y(j,+1)      &
                + ay(k,j ,i)*d1y(j, 0)      &
                + ay(k,jm,i)*d1y(j,-1)      &
                + az(kp,j,i)*d1z(k,+1)      &
                + az(k ,j,i)*d1z(k, 0)      &
                + az(km,j,i)*d1z(k,-1)
!
  end do
  end do
  end do

! ============================================================================
! the top, botom and sides  boundaries
! ============================================================================
  diva(0 ,:,:) = 0.0
  diva(NZ,:,:) = 0.0

! ============================================================================
! data exchage for Y direction 
! ============================================================================
  if(nproc_y == 1) then
     call mhd__period_y(diva)
  else
     call mpiut__exchange_y(diva)
     if(index_y == 0) then
        diva(:,0,:) = 0.0d0
     end if
     if(index_y == nproc_y-1) then 
        diva(:,ny,:) = 0.0d0
     end if
  end if

! ============================================================================
! data exchage for X direction 
! ============================================================================
  if(nproc_x == 1) then
     call mhd__period_x(diva)
  else
     call mpiut__exchange_x(diva)
     if(index_x == 0) then
        diva(:,:,0) = 0.0d0
     end if
     if(index_x == nproc_x-1) then 
        diva(:,:,nx) = 0.0d0
     end if
  end if

  end subroutine mhd__make_div     
!
!
! ============================================================================
! Grad Physical Value 
! ============================================================================
  subroutine mhd__grad(ax,grad_x,grad_y,grad_z)
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(in)  :: ax
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(out) :: grad_x, grad_y, grad_z

  do i = 0, NX-1
     ip=i+1
     im=i-1
  do j = 0, NY-1
     jp=j+1
     jm=j-1
  do k = 1, NZ-1
     kp=k+1
     km=k-1
!
     grad_x(k,j,i) =  (ax(k, j, ip )*d1x(i,+1)            &
                      +ax(k, j, i  )*d1x(i, 0)            &
                      +ax(k, j, im )*d1x(i,-1))             

     grad_y(k,j,i) =  (ax(k, jp, i )*d1y(j,+1)            &
                      +ax(k, j , i )*d1y(j, 0)            &
                      +ax(k, jm, i )*d1y(j,-1))             
 
     grad_z(k,j,i) =  (ax(kp ,j, i )*d1z(k,+1)            &
                      +ax(k  ,j, i )*d1z(k, 0)            &
                      +ax(km ,j, i )*d1z(k,-1))             
!
  end do
  end do
  end do

! -----------------------------------------------------------------------------
! top and botom boundaries
! -----------------------------------------------------------------------------
  do  i  = 0, nxm1
      ip = i+1
      im = i-1
  do  j  = 0, nym1
      jp = j+1
      jm = j-1
!                 
      k  = 0
      kp = k+1
      kpp= k+2
!
      grad_x(k,j,i) = (ax(k  ,j  ,ip )*d1x(i,+1)           &
                     + ax(k  ,j  ,i  )*d1x(i, 0)           &
                     + ax(k  ,j  ,im )*d1x(i,-1))              
!
      grad_y(k,j,i) = (ax(k  ,jp ,i  )*d1y(j,+1)           &
                     + ax(k  ,j  ,i  )*d1y(j, 0)           &
                     + ax(k  ,jm ,i  )*d1y(j,-1))              
!
      grad_z(k,j,i) = (ax(kpp,j  ,i  )*d1z(k,+1)           &
                     + ax(kp ,j  ,i  )*d1z(k, 0)           &
                     + ax(k  ,j  ,i  )*d1z(k,-1))              

        k = nz
       km = k-1
      kmm = k-2
!
      grad_x(k,j,i) = (ax(k  ,j  ,ip )*d1x(i,+1)          &
                     + ax(k  ,j  ,i  )*d1x(i, 0)          &
                     + ax(k  ,j  ,im )*d1x(i,-1))               
!
      grad_y(k,j,i) = (ax(k  ,jp ,i  )*d1y(j,+1)          &
                     + ax(k  ,j  ,i  )*d1y(j, 0)          &
                     + ax(k  ,jm ,i  )*d1y(j,-1))         
!
      grad_z(k,j,i) = (ax(k  ,j  ,i  )*d1z(k,+1)          &
                     + ax(km ,j  ,i  )*d1z(k, 0)          &
                     + ax(kmm,j  ,i  )*d1z(k,-1))         
!
  end do
  end do

! ============================================================================
! data exchage for Y direction 
! ============================================================================
!  if(nproc_y == 1) then
!     call mhd__period_y(grad_x)
!     call mhd__period_y(grad_y)
!     call mhd__period_y(grad_z)
!  else
     call mpiut__exchange_y(grad_x,grad_y,grad_z)
        
! ----------------------------------------------------------------------------
  if(index_y == 0) then
! ----------------------------------------------------------------------------
  do i  = 0, nxm1
     ip = i+1
     im = i-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     j  = 0
     jp = j+1
     jpp= j+2

     grad_x(k,j,i) = ( ax(k  ,j  ,ip )*d1x(i,+1)           &
                     + ax(k  ,j  ,i  )*d1x(i, 0)           &
                     + ax(k  ,j  ,im )*d1x(i,-1))              
!
     grad_y(k,j,i) = ( ax(k  ,jpp ,i )*d1y(j,+1)           &
                     + ax(k  ,jp  ,i )*d1y(j, 0)           &
                     + ax(k  ,j   ,i )*d1y(j,-1))              
!
     grad_z(k,j,i) = ( ax(kp ,j  ,i  )*d1z(k,+1)           &
                     + ax(k  ,j  ,i  )*d1z(k, 0)           &
                     + ax(km ,j  ,i  )*d1z(k,-1))        


  end do
  end do

! ----------------------------------------------------------------------------
  else if(index_y == nproc_y-1) then
! ----------------------------------------------------------------------------
  do i  = 0, nxm1
     ip = i+1
     im = i-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     j  = ny
     jm = j-1
     jmm= j-2

     grad_x(k,j,i) = ( ax(k  ,j  ,ip )*d1x(i,+1)           &
                     + ax(k  ,j  ,i  )*d1x(i, 0)           &
                     + ax(k  ,j  ,im )*d1x(i,-1))              
!
     grad_y(k,j,i) = ( ax(k  ,j   ,i )*d1y(j,+1)           &
                     + ax(k  ,jm  ,i )*d1y(j, 0)           &
                     + ax(k  ,jmm ,i )*d1y(j,-1))              
!
     grad_z(k,j,i) = ( ax(kp ,j  ,i  )*d1z(k,+1)           &
                     + ax(k  ,j  ,i  )*d1z(k, 0)           &
                     + ax(km ,j  ,i  )*d1z(k,-1))      

  end do
  end do
      
  end if

!  end if

! ============================================================================
! data exchage for X direction 
! ============================================================================
!  if(nproc_x == 1) then
!     call mhd__period_x(grad_x)
!     call mhd__period_x(grad_y)
!     call mhd__period_x(grad_z)
!  else
     call mpiut__exchange_x(grad_x,grad_y,grad_z)

! ----------------------------------------------------------------------------
  if(index_x == 0) then
! ----------------------------------------------------------------------------
  do j  = 0, nym1
     jp = j+1
     jm = j-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     i  = 0
     ip = i+1
     ipp= i+2
!
     grad_x(k,j,i) = ( ax(k  ,j  ,ipp)*d1x(i,+1)           &
                     + ax(k  ,j  ,ip )*d1x(i, 0)           &
                     + ax(k  ,j  ,i  )*d1x(i,-1))              
!
     grad_y(k,j,i) = (ax(k  ,jp ,i  )*d1y(j,+1)           &
                    + ax(k  ,j  ,i  )*d1y(j, 0)           &
                    + ax(k  ,jm ,i  )*d1y(j,-1))              
!
     grad_z(k,j,i) = (ax(kp ,j  ,i  )*d1z(k,+1)           &
                    + ax(k  ,j  ,i  )*d1z(k, 0)           &
                    + ax(km ,j  ,i  )*d1z(k,-1))        
!
!
  end do
  end do
!
! ----------------------------------------------------------------------------
  else if(index_x == nproc_x-1) then
! ----------------------------------------------------------------------------
  do j  = 0, nym1
     jp = j+1
     jm = j-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     i  = nx
     im = i-1
     imm= i-2
!
     grad_x(k,j,i) = (ax(k  ,j  ,i  )*d1x(i,+1)           &
                     +ax(k  ,j  ,im )*d1x(i, 0)           &
                     +ax(k  ,j  ,imm)*d1x(i,-1))              
!
     grad_y(k,j,i) = (ax(k  ,jp  ,i )*d1y(j,+1)           &
                     +ax(k  ,j   ,i )*d1y(j, 0)           &
                     +ax(k  ,jm  ,i )*d1y(j,-1))              
!
     grad_z(k,j,i) = (ax(kp ,j  ,i  )*d1z(k,+1)           &
                     +ax(k  ,j  ,i  )*d1z(k, 0)           &
                     +ax(km ,j  ,i  )*d1z(k,-1))      
!
  end do
  end do
!       
  end if
!
!  end if

  end subroutine mhd__grad

! ============================================================================
! Curl Physical Value
! ============================================================================
  subroutine mhd__make_curl(ax,ay,az,curlx,curly,curlz)
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(in)  :: ax, ay, az
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(out) :: curlx, curly, curlz
!
  do i = 0, nxm1
     ip = i+1
     im = i-1
  do j = 0, nym1
     jp = j+1
     jm = j-1
  do k = 1, nzm1
     kp = k+1
     km = k-1
!
     curlx(k,j,i) = (az(k ,jp,i )*d1y(j,+1)  &
                    +az(k, j ,i )*d1y(j, 0)  &
                    +az(k ,jm,i )*d1y(j,-1)) &
                   -(ay(kp,j ,i )*d1z(k,+1)  &
                    +ay(k ,j ,i )*d1z(k, 0)  &
                    +ay(km,j ,i )*d1z(k,-1)) 


     curly(k,j,i) = (ax(kp,j ,i )*d1z(k,+1)  &
                    +ax(k ,j ,i )*d1z(k, 0)  &
                    +ax(km,j ,i )*d1z(k,-1)) &
                   -(az(k ,j ,ip)*d1x(i,+1)  &
                    +az(k ,j ,i )*d1x(i, 0)  &
                    +az(k ,j ,im)*d1x(i,-1)) 


     curlz(k,j,i) = (ay(k ,j ,ip)*d1x(i,+1)  &
                    +ay(k ,j ,i )*d1x(i, 0)  &
                    +ay(k ,j ,im)*d1x(i,-1)) &
                   -(ax(k ,jp,i )*d1y(j,+1)  &
                    +ax(k ,j ,i )*d1y(j, 0)  &
                    +ax(k ,jm,i )*d1y(j,-1)) 

!
  end do
  end do
  end do

! ----------------------------------------------------------------------------
! top & bottom boundary
! ----------------------------------------------------------------------------
  do i = 0, nxm1
     ip = i+1
     im = i-1
  do j = 0, nym1
     jp = j+1
     jm = j-1
!
     k  = 0
     kp = k+1
     kpp= k+2
!
     curlx(k,j,i) = (az(k  ,jp ,i  )*d1y(j,+1)     &
                    +az(k  , j ,i  )*d1y(j, 0)     &
                    +az(k  ,jm ,i  )*d1y(j,-1))    &
                   -(ay(kpp,j  ,i  )*d1z(k,+1)     &
                    +ay(kp, j  ,i  )*d1z(k, 0)     &
                    +ay(k , j  ,i  )*d1z(k,-1))    
  

     curly(k,j,i) = (ax(kpp,j  ,i  )*d1z(k,+1)  &
                    +ax(kp ,j  ,i  )*d1z(k, 0)  &
                    +ax(k  ,j  ,i  )*d1z(k,-1)) &
                   -(az(k  ,j  ,ip )*d1x(i,+1)  &
                    +az(k  ,j  ,i  )*d1x(i, 0)  &
                    +az(k  ,j  ,im )*d1x(i,-1)) 

     curlz(k,j,i) = (ay(k  ,j  ,ip )*d1x(i,+1)  &
                    +ay(k  ,j  ,i  )*d1x(i, 0)  &
                    +ay(k  ,j  ,im )*d1x(i,-1)) &
                   -(ax(k  ,jp ,i  )*d1y(j,+1)  &
                    +ax(k  ,j  ,i  )*d1y(j, 0)  &
                    +ax(k  ,jm ,i  )*d1y(j,-1)) 

!
     k = nz
     km = k-1
     kmm= k-2
!
     curlx(k,j,i) = (az(k  ,jp ,i  )*d1y(j,+1)  &
                    +az(k,  j,  i  )*d1y(j, 0)  &
                    +az(k  ,jm ,i  )*d1y(j,-1)) &
                   -(ay(k  ,j  ,i  )*d1z(k,+1)  &
                    +ay(km ,j  ,i  )*d1z(k, 0)  &
                    +ay(kmm,j  ,i  )*d1z(k,-1)) 

     curly(k,j,i) = (ax(k  ,j  ,i  )*d1z(k,+1)  &
                    +ax(km ,j  ,i  )*d1z(k, 0)  &
                    +ax(kmm,j  ,i  )*d1z(k,-1)) &
                   -(az(k  ,j  ,ip )*d1x(i,+1)  &
                    +az(k  ,j  ,i  )*d1x(i, 0)  &
                    +az(k  ,j  ,im )*d1x(i,-1)) 

     curlz(k,j,i) = (ay(k  ,j  ,ip )*d1x(i,+1)  &
                    +ay(k  ,j  ,i  )*d1x(i, 0)  &
                    +ay(k  ,j  ,im )*d1x(i,-1)) &
                   -(ax(k  ,jp ,i  )*d1y(j,+1)  &
                    +ax(k,  j  ,i  )*d1y(j, 0)  &
                    +ax(k  ,jm ,i  )*d1y(j,-1))

!
  end do
  end do

! ============================================================================
! data exchage for Y direction and side boundary
! ============================================================================
!  if(nproc_y == 1) then
!     call mhd__period_y(curlx)
!     call mhd__period_y(curly) 
!     call mhd__period_y(curlz)
!  else
     call mpiut__exchange_y(curlx,curly,curlz)

! ----------------------------------------------------------------------------
  if(index_y == 0) then
! ----------------------------------------------------------------------------
  do i  = 0, nxm1
     ip = i+1
     im = i-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     j  = 0
     jp = j+1
     jpp= j+2
!
     curlx(k,j,i) = (az(k  ,jpp,i )*d1y(j,+1)   &
                    +az(k  ,jp ,i )*d1y(j, 0)   &
                    +az(k  ,j  ,i )*d1y(j,-1))  &
                   -(ay(kp ,j  ,i )*d1z(k,+1)   &
                    +ay(k  ,j  ,i )*d1z(k, 0)   &
                    +ay(km ,j  ,i )*d1z(k,-1))  

     curly(k,j,i) = (ax(kp ,j ,i  )*d1z(k,+1)   &
                    +ax(k  ,j ,i  )*d1z(k, 0)   &
                    +ax(km ,j ,i  )*d1z(k,-1))  &
                   -(az(k  ,j ,ip )*d1x(i,+1)   &
                    +az(k  ,j ,i  )*d1x(i, 0)   &
                    +az(k  ,j ,im )*d1x(i,-1))   

     curlz(k,j,i) = (ay(k ,j  ,ip )*d1x(i,+1)   &
                    +ay(k ,j  ,i  )*d1x(i, 0)   &
                    +ay(k ,j  ,im )*d1x(i,-1))  &
                   -(ax(k ,jpp,i  )*d1y(j,+1)   &
                    +ax(k ,jp ,i  )*d1y(j, 0)   &
                    +ax(k ,j  ,i  )*d1y(j,-1))   
!
 end do
 end do

! ----------------------------------------------------------------------------
  else if(index_y == nproc_y-1) then
! ---------------------------------------------------------------------------
  do i  = 0, nxm1
     ip = i+1
     im = i-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     j  = ny
     jm = j-1
     jmm= j-2
!
     curlx(k,j,i) = (az(k ,j  ,i )*d1y(j,+1)     &
                    +az(k ,jm ,i )*d1y(j, 0)     &
                    +az(k ,jmm,i )*d1y(j,-1))    &
                   -(ay(kp,j  ,i )*d1z(k,+1)     &
                    +ay(k ,j  ,i )*d1z(k, 0)     &
                    +ay(km,j  ,i )*d1z(k,-1))  

     curly(k,j,i) = (ax(kp,j ,i  )*d1z(k,+1)     &
                    +ax(k ,j ,i  )*d1z(k, 0)     &
                    +ax(km,j ,i  )*d1z(k,-1))    & 
                   -(az(k ,j ,ip )*d1x(i,+1)     &
                    +az(k ,j ,i  )*d1x(i, 0)     &
                    +az(k ,j ,im )*d1x(i,-1))  
!
     curlz(k,j,i) = (ay(k ,j  ,ip )*d1x(i,+1)    &
                    +ay(k ,j  ,i  )*d1x(i, 0)    &
                    +ay(k ,j  ,im )*d1x(i,-1))   &
                   -(ax(k ,j  ,i  )*d1y(j,+1)    &
                    +ax(k, jm ,i  )*d1y(j, 0)    &
                    +ax(k ,jmm,i  )*d1y(j,-1))
!
  end do
  end do         

  end if

!  end if

! ============================================================================
! data exchage for X direction 
! ============================================================================
!  if(nproc_x == 1) then
!     call mhd__period_x(curlx)
!     call mhd__period_x(curly)
!     call mhd__period_x(curlz)
!  else
     call mpiut__exchange_x(curlx,curly,curlz)

! ----------------------------------------------------------------------------
  if(index_x == 0) then
! ----------------------------------------------------------------------------
  do j  = 0, nym1
     jp = j+1
     jm = j-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     i  = 0
     ip = i+1
     ipp= i+2
!
     curlx(k,j,i) = ( az(k ,jp ,i )*d1y(j,+1 )  &
                    + az(k , j ,i )*d1y(j, 0 )  &
                    + az(k ,jm ,i )*d1y(j,-1 )) &
                   -( ay(kp ,j, i )*d1z(k,+1 )  &
                    + ay(k , j ,i )*d1z(k, 0 )  &
                    + ay(km ,j ,i )*d1z(k,-1 ))  
!
     curly(k,j,i) = ( ax(kp,j ,i  )*d1z(k,+1 )  &
                    + ax(k ,j ,i  )*d1z(k, 0 )  &
                    + ax(km,j ,i  )*d1z(k,-1 )) &
                   -( az(k ,j ,ipp)*d1x(i,+1 )  &
                    + az(k ,j ,ip )*d1x(i, 0 )  &
                    + az(k ,j ,i  )*d1x(i,-1 ))   
!
     curlz(k,j,i) = ( ay(k ,j ,ipp)*d1x(i,+1 )  &
                    + ay(k ,j ,ip )*d1x(i, 0 )  &
                    + ay(k ,j ,i  )*d1x(i,-1 )) &
                   -( ax(k ,jp,i  )*d1y(j,+1 )  &
                    + ax(k ,j ,i  )*d1y(j, 0 )  &
                    + ax(k ,jm,i  )*d1y(j,-1 ))   
!
  end do
  end do
!
! ----------------------------------------------------------------------------
  else if(index_x == nproc_x-1) then
! ----------------------------------------------------------------------------
  do j  = 0, nym1
     jp = j+1
     jm = j-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     i  = nx
     im = i-1
     imm= i-2
!
     curlx(k,j,i) = ( az(k  ,jp ,i  )*d1y(j,+1 )     &
                    + az(k  ,j  ,i  )*d1y(j, 0 )     &
                    + az(k  ,jm ,i  )*d1y(j,-1 ))    &
                   -( ay(kp ,j  ,i  )*d1z(k,+1 )     &
                    + ay(k  ,j  ,i  )*d1z(k, 0 )     &
                    + ay(km ,j  ,i  )*d1z(k,-1 ))  

     curly(k,j,i) = ( ax(kp ,j  ,i  )*d1z(k,+1 )     &
                    + ax(k  ,j  ,i  )*d1z(k, 0 )     &
                    + ax(km ,j  ,i  )*d1z(k,-1 ))    & 
                   -( az(k  ,j  ,i  )*d1x(i,+1 )     &
                    + az(k  ,j  ,im )*d1x(i, 0 )     &
                    + az(k  ,j  ,imm)*d1x(i,-1 ))  

     curlz(k,j,i) = ( ay(k  ,j  ,i  )*d1x(i,+1 )     &
                    + ay(k  ,j  ,im )*d1x(i, 0 )     &
                    + ay(k  ,j  ,imm)*d1x(i,-1 ))    &
                   -( ax(k  ,jp ,i  )*d1y(j,+1 )     &
                    + ax(k  ,j  ,i  )*d1y(j, 0 )     &
                    + ax(k  ,jm ,i  )*d1y(j,-1 ))
!
  end do
  end do         
!
  end if

!  end if

  end subroutine mhd__make_curl

! ============================================================================
  subroutine mhd__make_laplace(ax,ay,az,grad2_x,grad2_y,grad2_z)
! ============================================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(in)  :: ax, ay, az
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(out) :: grad2_x, & 
                                                        grad2_y, & 
                                                        grad2_z 
  do i = 0, nxm1
     ip = i+1
     im = i-1
  do j = 0, nym1
     jp = j+1
     jm = j-1
  do k = 1, nzm1
     kp = k+1
     km = k-1

     grad2_x(k,j,i) = (ax(k ,j ,ip)*d2x(i,+1)  &
                     + ax(k ,j ,i )*d2x(i, 0)  &
                     + ax(k ,j ,im)*d2x(i,-1)) &
                     +(ax(k ,jp,i )*d2y(j,+1)  &
                     + ax(k ,j ,i )*d2y(j, 0)  &
                     + ax(k ,jm,i )*d2y(j,-1)) & 
                     +(ax(kp,j ,i )*d2z(k,+1)  &
                     + ax(k ,j ,i )*d2z(k, 0)  &
                     + ax(km,j ,i )*d2z(k,-1)) 

     grad2_y(k,j,i) = (ay(k ,j ,ip)*d2x(i,+1)  &
                     + ay(k ,j ,i )*d2x(i, 0)  &
                     + ay(k ,j ,im)*d2x(i,-1)) &
                     +(ay(k ,jp,i )*d2y(j,+1)  &
                     + ay(k ,j ,i )*d2y(j, 0)  &
                     + ay(k ,jm,i )*d2y(j,-1)) &
                     +(ay(kp,j ,i )*d2z(k,+1)  &
                     + ay(k ,j ,i )*d2z(k, 0)  &
                     + ay(km,j ,i )*d2z(k,-1)) 

     grad2_z(k,j,i) = (az(k ,j ,ip)*d2x(i,+1)  &
                     + az(k ,j ,i )*d2x(i, 0)  &
                     + az(k ,j ,im)*d2x(i,-1)) &
                     +(az(k ,jp,i )*d2y(j,+1)  &
                     + az(k ,j ,i )*d2y(j, 0)  &
                     + az(k ,jm,i )*d2y(j,-1)) &
                     +(az(kp,j ,i )*d2z(k,+1)  &
                     + az(k ,j ,i )*d2z(k, 0)  &
                     + az(km,j ,i )*d2z(k,-1)) 

!
  end do
  end do
  end do
!
! -----------------------------------------------------------------------------
!                     the top & bottom boundary 
! -----------------------------------------------------------------------------
  do i = 0, nxm1
     ip = i+1
     im = i-1
  do j = 0, nym1
     jp = j+1
     jm = j-1
!
     k  = 0
     kp = k+1
     kpp= k+2
!
     grad2_x(k,j,i) = (ax(k  ,j  ,ip )*d2x(i,+1)   &
                      +ax(k  ,j  ,i  )*d2x(i, 0)   &
                      +ax(k  ,j  ,im )*d2x(i,-1))  &
                     +(ax(k  ,jp ,i  )*d2y(j,+1)   &
                      +ax(k  ,j  ,i  )*d2y(j, 0)   &
                      +ax(k  ,jm ,i  )*d2y(j,-1))  &
                     +(ax(kpp,j  ,i  )*d2z(k,+1)   &
                      +ax(kp ,j  ,i  )*d2z(k, 0)   &
                      +ax(k  ,j  ,i  )*d2z(k,-1))  

     grad2_y(k,j,i) = (ay(k  ,j  ,ip )*d2x(i,+1)  &
                      +ay(k  ,j  ,i  )*d2x(i, 0)  &
                      +ay(k  ,j  ,im )*d2x(i,-1)) &
                     +(ay(k  ,jp ,i  )*d2y(j,+1)  &
                      +ay(k  ,j  ,i  )*d2y(j, 0)  &
                      +ay(k  ,jm ,i  )*d2y(j,-1)) &
                     +(ay(kpp,j  ,i  )*d2z(k,+1)  &
                      +ay(kp ,j  ,i  )*d2z(k, 0)  &
                      +ay(k  ,j  ,i  )*d2z(k,-1)) 

     grad2_z(k,j,i) = (az(k  ,j  ,ip )*d2x(i,+1)  &
                      +az(k  ,j  ,i  )*d2x(i, 0)  &
                      +az(k  ,j  ,im )*d2x(i,-1)) &
                     +(az(k  ,jp ,i  )*d2y(j,+1)  &
                      +az(k  ,j  ,i  )*d2y(j, 0)  &
                      +az(k  ,jm ,i  )*d2y(j,-1)) &
                     +(az(kpp,j  ,i  )*d2z(k,+1)  &
                      +az(kp ,j  ,i  )*d2z(k, 0)  &
                      +az(k  ,j  ,i  )*d2z(k,-1)) 
!
     k = nz
     km = k-1
     kmm= k-2
!
     grad2_x(k,j,i) = (ax(k  ,j  ,ip )*d2x(i,+1)  &
                      +ax(k  ,j  ,i  )*d2x(i, 0)  &
                      +ax(k  ,j  ,im )*d2x(i,-1)) &
                     +(ax(k  ,jp ,i  )*d2y(j,+1)  &
                      +ax(k  ,j  ,i  )*d2y(j, 0)  &
                      +ax(k  ,jm ,i  )*d2y(j,-1)) &
                     +(ax(k  ,j  ,i  )*d2z(k,+1)  &
                      +ax(km ,j  ,i  )*d2z(k, 0)  &
                      +ax(kmm,j  ,i  )*d2z(k,-1)) 

     grad2_y(k,j,i) = (ay(k  ,j  ,ip )*d2x(i,+1)  &
                      +ay(k  ,j  ,i  )*d2x(i, 0)  &
                      +ay(k  ,j  ,im )*d2x(i,-1)) &
                     +(ay(k  ,jp ,i  )*d2y(j,+1)  &
                      +ay(k  ,j  ,i  )*d2y(j, 0)  &
                      +ay(k  ,jm ,i  )*d2y(j,-1)) &
                     +(ay(k  ,j  ,i  )*d2z(k,+1)  &
                      +ay(km ,j  ,i  )*d2z(k, 0)  &
                      +ay(kmm,j  ,i  )*d2z(k,-1)) 

     grad2_z(k,j,i) = (az(k  ,j  ,ip )*d2x(i,+1)  &
                      +az(k  ,j  ,i  )*d2x(i, 0)  &
                      +az(k  ,j  ,im )*d2x(i,-1)) &
                     +(az(k  ,jp ,i  )*d2y(j,+1)  &
                      +az(k,  j  ,i  )*d2y(j, 0)  &
                      +az(k  ,jm ,i  )*d2y(j,-1)) &
                     +(az(k  ,j  ,i  )*d2z(k,+1)  &
                      +az(km ,j  ,i  )*d2z(k, 0)  &
                      +az(kmm,j  ,i  )*d2z(k,-1)) 

   end do
   end do
!
! ============================================================================
! data exchage for Y direction and side boundary
! ============================================================================
!  if(nproc_y == 1) then
!     call mhd__period_y(grad2_x)
!     call mhd__period_y(grad2_y)
!     call mhd__period_y(grad2_z)
!  else
     call mpiut__exchange_y(grad2_x,grad2_y,grad2_z)

! ----------------------------------------------------------------------------
  if(index_y == 0) then
! ----------------------------------------------------------------------------
  do i  = 0, nxm1
     ip = i+1
     im = i-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     j  = 0
     jp = j+1
     jpp= j+2

!
     grad2_x(k,j,i) = (ax(k  ,j  ,ip)*d2x(i,+1)   &
                      +ax(k  ,j  ,i )*d2x(i, 0)   &
                      +ax(k  ,j  ,im)*d2x(i,-1))  &
                     +(ax(k  ,jpp,i )*d2y(j,+1)   &
                      +ax(k  ,jp ,i )*d2y(j, 0)   &
                      +ax(k  ,j  ,i )*d2y(j,-1))  & 
                     +(ax(kp ,j  ,i )*d2z(k,+1)   &
                      +ax(k  ,j  ,i )*d2z(k, 0)   &
                      +ax(km ,j  ,i )*d2z(k,-1)) 
!
     grad2_y(k,j,i) = (ay(k  ,j  ,ip)*d2x(i,+1)   &
                      +ay(k  ,j  ,i )*d2x(i, 0)   &
                      +ay(k  ,j  ,im)*d2x(i,-1))  &
                     +(ay(k  ,jpp,i )*d2y(j,+1)   &
                      +ay(k  ,jp ,i )*d2y(j, 0)   &
                      +ay(k  ,j  ,i )*d2y(j,-1))  &
                     +(ay(kp ,j  ,i )*d2z(k,+1)   &
                      +ay(k  ,j  ,i )*d2z(k, 0)   &
                      +ay(km ,j  ,i )*d2z(k,-1)) 

!
     grad2_z(k,j,i) = (az(k  ,j  ,ip)*d2x(i,+1)   &
                      +az(k  ,j  ,i )*d2x(i, 0)   &
                      +az(k  ,j  ,im)*d2x(i,-1))  &
                     +(az(k  ,jpp,i )*d2y(j,+1)   & 
                      +az(k  ,jp ,i )*d2y(j, 0)   &
                      +az(k  ,j  ,i )*d2y(j,-1))  &  
                     +(az(kp ,j  ,i )*d2z(k,+1)   &
                      +az(k  ,j  ,i )*d2z(k, 0)   &
                      +az(km ,j  ,i )*d2z(k,-1)) 
!
  end do
  end do

! ----------------------------------------------------------------------------
  else if(index_y == nproc_y-1) then
! ----------------------------------------------------------------------------
  do i  = 0, nxm1
     ip = i+1
     im = i-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     j  = ny
     jm = j-1
     jmm= j-2

     grad2_x(k,j,i) = (ax(k ,j  ,ip)*d2x(i,+1)   &
                     + ax(k ,j  ,i )*d2x(i, 0)   &
                     + ax(k ,j  ,im)*d2x(i,-1))  &
                     +(ax(k ,j  ,i )*d2y(j,+1)   &
                     + ax(k ,jm ,i )*d2y(j, 0)   &
                     + ax(k ,jmm,i )*d2y(j,-1))  & 
                     +(ax(kp,j  ,i )*d2z(k,+1)   &
                     + ax(k ,j  ,i )*d2z(k, 0)   &
                     + ax(km,j  ,i )*d2z(k,-1))
!
     grad2_y(k,j,i) = (ay(k ,j  ,ip)*d2x(i,+1)   &
                     + ay(k ,j  ,i )*d2x(i, 0)   &
                     + ay(k ,j  ,im)*d2x(i,-1))  & 
                     +(ay(k ,j  ,i )*d2y(j,+1)   &
                     + ay(k ,jm ,i )*d2y(j, 0)   &
                     + ay(k ,jmm,i )*d2y(j,-1))  &
                     +(ay(kp,j  ,i )*d2z(k,+1)   &
                     + ay(k ,j  ,i )*d2z(k, 0)   &
                     + ay(km,j  ,i )*d2z(k,-1))  

     grad2_z(k,j,i) = (az(k ,j  ,ip)*d2x(i,+1)   &
                     + az(k ,j  ,i )*d2x(i, 0)   &
                     + az(k ,j  ,im)*d2x(i,-1))  &
                     +(az(k ,j  ,i )*d2y(j,+1)   &
                     + az(k ,jm ,i )*d2y(j, 0)   &
                     + az(k ,jmm,i )*d2y(j,-1))  &
                     +(az(kp,j  ,i )*d2z(k,+1)   &
                     + az(k ,j  ,i )*d2z(k, 0)   &
                     + az(km,j  ,i )*d2z(k,-1)) 
!
  end do
  end do         
  end if
!
!  end if
! =============================================================================
! set boundary & data exchage for X direction 
! =============================================================================
!  if(nproc_x == 1) then
!     call mhd__period_x(grad2_x)
!     call mhd__period_x(grad2_y)
!     call mhd__period_x(grad2_z)
!  else
     call mpiut__exchange_x(grad2_x,grad2_y,grad2_z)

! -----------------------------------------------------------------------------
  if(index_x == 0) then
! -----------------------------------------------------------------------------
  do j  = 0, nym1
     jp = j+1
     jm = j-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     i  = 0
     ip = i+1
     ipp= i+2

     grad2_x(k,j,i) = (ax(k ,j  ,ipp)*d2x(i,+1)  &
                     + ax(k ,j  ,ip )*d2x(i, 0)  &
                     + ax(k ,j  ,i  )*d2x(i,-1)) &
                     +(ax(k ,jp ,i  )*d2y(j,+1)  &
                     + ax(k ,j  ,i  )*d2y(j, 0)  &
                     + ax(k ,jm ,i  )*d2y(j,-1)) &
                     +(ax(kp,j  ,i  )*d2z(k,+1)  &
                     + ax(k ,j  ,i  )*d2z(k, 0)  &
                     + ax(km,j  ,i  )*d2z(k,-1))
!
     grad2_y(k,j,i) = (ay(k ,j  ,ipp)*d2x(i,+1)  &
                     + ay(k ,j  ,ip )*d2x(i, 0)  &
                     + ay(k ,j  ,i  )*d2x(i,-1)) &
                     +(ay(k ,jp ,i  )*d2y(j,+1)  &
                     + ay(k ,j  ,i  )*d2y(j, 0)  &
                     + ay(k ,jm ,i  )*d2y(j,-1)) &   
                     +(ay(kp,j  ,i  )*d2z(k,+1)  &
                     + ay(k ,j  ,i  )*d2z(k, 0)  &
                     + ay(km,j  ,i  )*d2z(k,-1))   

!
     grad2_z(k,j,i) = (az(k ,j  ,ipp)*d2x(i,+1)  &
                     + az(k ,j  ,ip )*d2x(i, 0)  &
                     + az(k ,j  ,i  )*d2x(i,-1)) &
                     +(az(k ,jp ,i  )*d2y(j,+1)  &
                     + az(k ,j  ,i  )*d2y(j, 0)  &
                     + az(k ,jm ,i  )*d2y(j,-1)) &
                     +(az(kp,j  ,i  )*d2z(k,+1)  &
                     + az(k ,j  ,i  )*d2z(k, 0)  &
                     + az(km,j  ,i  )*d2z(k,-1))
   
!
  end do
  end do

! -----------------------------------------------------------------------------
  else if(index_x == nproc_x-1) then
! -----------------------------------------------------------------------------
  do j  = 0, nym1
     jp = j+1
     jm = j-1
  do k =  1, nzm1
     kp = k+1
     km = k-1
!
     i  = nx
     im = i-1
     imm= i-2
!
     grad2_x(k,j,i) = (ax(k  ,j  ,i  )*d2x(i,+1)   &
                     + ax(k  ,j  ,im )*d2x(i, 0)   &
                     + ax(k  ,j  ,imm)*d2x(i,-1))  &
                     +(ax(k  ,jp ,i  )*d2y(j,+1)   &
                     + ax(k  ,j  ,i  )*d2y(j, 0)   &
                     + ax(k  ,jm ,i  )*d2y(j,-1))  &
                     +(ax(kp ,j  ,i  )*d2z(k,+1)   &
                     + ax(k  ,j  ,i  )*d2z(k, 0)   &
                     + ax(km ,j  ,i  )*d2z(k,-1))  
!
     grad2_y(k,j,i) = (ay(k  ,j  ,i  )*d2x(i,+1)   &
                     + ay(k  ,j  ,im )*d2x(i, 0)   &
                     + ay(k  ,j  ,imm)*d2x(i,-1))  & 
                     +(ay(k  ,jp ,i  )*d2y(j,+1)   &
                     + ay(k  ,j  ,i  )*d2y(j, 0)   &
                     + ay(k  ,jm ,i  )*d2y(j,-1))  &
                     +(ay(kp ,j  ,i  )*d2z(k,+1)   &
                     + ay(k  ,j  ,i  )*d2z(k, 0)   &
                     + ay(km ,j  ,i  )*d2z(k,-1))  
!
     grad2_z(k,j,i) = (az(k  ,j  ,i  )*d2x(i,+1)   &
                     + az(k  ,j  ,im )*d2x(i, 0)   &
                     + az(k  ,j  ,imm)*d2x(i,-1))  &
                     +(az(k  ,jp ,i  )*d2y(j,+1)   &
                     + az(k  ,j  ,i  )*d2y(j, 0)   &
                     + az(k  ,jm ,i  )*d2y(j,-1))  &
                     +(az(kp ,j  ,i  )*d2z(k,+1)   &
                     + az(k  ,j  ,i  )*d2z(k, 0)   &
                     + az(km ,j  ,i  )*d2z(k,-1))

  end do
  end do         
!
  end if
!
!  end if
   
  end subroutine mhd__make_laplace

! =============================================================================
  subroutine mhd__period_x(a)
! =============================================================================
  use constants
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(INOUT) :: a
     a(:,:,NX) = a(:,:,0)
     a(:,:,-1) = a(:,:,NX-1)
  return
  end subroutine mhd__period_x

! =============================================================================
  subroutine mhd__period_y(a)
! =============================================================================
  use constants
  real(DP), dimension(0:NZ,-1:NY,-1:NX), intent(INOUT) :: a
      a(:,NY,:) = a(:,0,:)
      a(:,-1,:) = a(:,NY-1,:)
  return
  end subroutine mhd__period_y

! ============================================================================
  subroutine mhd__set_integral ! Physical Values are integrated
! ============================================================================
  real(DP),dimension(0:nz,-1:ny,-1:nx) :: ene_mag, ene_kin

! ----------------------------------------------------------------------------
! Magnetic & Kinetic Energy
! ----------------------------------------------------------------------------
  do i = -1,nx
  do j = -1,ny
  do k =  0,nz
     ene_mag(k,j,i) = 0.5*   (bx(k,j,i)**2 +  by(k,j,i)**2 +   bz(k,j,i)**2)
     ene_kin(k,j,i) = 0.5*ro(k,j,i)*(vx(k,j,i)**2 +  & 
                                     vy(k,j,i)**2 +  &
                                     vz(k,j,i)**2)
  end do
  end do
  end do

  call pset__integrate(ene_mag, ave_ene_mag)
  call pset__integrate(ene_kin, ave_ene_kin)

  end subroutine mhd__set_integral

! -----------------------------------------------------------------------------
  end module mhd















