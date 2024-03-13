  program main
! =====================PROGRAM merge========================================
! This is used to marge the data in sub-domained data.
! ==========================================================================
  use subroutines
  use ana
  character*3   :: cseries
  character*4   :: cproc
  character*100 :: filename
  integer       :: series=1, proc, p_x, p_y

! ===========================================================================
!                       READ NAMELIST
! ===========================================================================
  call pset

! ===========================================================================
!                     check proc number 
! ===========================================================================
  nproc = nproc_x*nproc_y
  if(N_PROC .ne. nproc) then
  write(*,*) '##ERR(main): N_PROC \= nproc_x*nproc_y',      & 
              N_PROC, nproc_x,nproc_y
  stop
  end if

  if(N_PROC_X .ne. nproc_x) then
  write(*,*) '##ERR(main): N_PROC_X \= nproc_x',            &
              N_PROC_X, nproc_x
  stop
  end if

  if(N_PROC_Y .ne. nproc_y) then
  write(*,*) '##ERR(main): N_PROC_Y \= nproc_y',            &
              N_PROC_Y, nproc_y
  stop
  end if

! ============================================================================
!                 READ DATA
! ============================================================================
  write(cseries,'(i3.3)') series
  write(*,*) '## read file #',series
!
  do loop = 1,11
  write(cloop,'(i3.3)') loop
  write(6,*)cloop
!
  do proc = 0, N_PROC-1
  call domain_coord(proc, p_x, p_y)
  write(cproc,'(i4.4)') proc

! ---------------------------------------------------------------------------
! Magnetic Field
! ---------------------------------------------------------------------------
  open(FILE_RESTART,file=dir_r//trim(cfile_restart)//'_'//                 &
                                trim(cloop)//'.'//cproc,form='unformatted')
  read(FILE_RESTART) iwrite, nloop, atime, dtstep, bx, by, bz, ph
  close(FILE_RESTART)

! ---------------------------------------------------------------------------
! Velocity Field
! ---------------------------------------------------------------------------
  open(11,file=dir_r//'B3D_RESTART2'//'_'//                 &
                                trim(cloop)//'.'//cproc,form='unformatted')
  read(11) vx, vy, vz, ro
  close(11)

! ===========================================================================
! calculate J & DivB
! ===========================================================================
  call   calJ(p_x,p_y)
  call div__b(p_x,p_y)

! ===========================================================================
! copy to the global array
! ===========================================================================
  call copy_to_sparseglobal(bx,bx_r,p_x,p_y)
  call copy_to_sparseglobal(by,by_r,p_x,p_y)
  call copy_to_sparseglobal(bz,bz_r,p_x,p_y)
  call copy_to_sparseglobal(vx,vx_r,p_x,p_y)
  call copy_to_sparseglobal(vy,vy_r,p_x,p_y)
  call copy_to_sparseglobal(vz,vz_r,p_x,p_y)
  call copy_to_sparseglobal(cx,cx_r,p_x,p_y)
  call copy_to_sparseglobal(cy,cy_r,p_x,p_y)
  call copy_to_sparseglobal(cz,cz_r,p_x,p_y)
  call copy_to_sparseglobal(ro,ro_r,p_x,p_y)
  call copy_to_sparseglobal(ph,ph_r,p_x,p_y)
  call copy_to_sparseglobal(divb,divb_r,p_x,p_y)
  call copy_to_sparseglobal(divv,divv_r,p_x,p_y)
  call copy_to_sparseglobal(cb2 ,cb2_r ,p_x,p_y)

  end do


! ========================================================================
! find maximum value
! ========================================================================
  write(FILE_SYSOUT,*) ' bx      (max, min) = ',maxval(bx_r),minval(bx_r)
  write(FILE_SYSOUT,*) ' by      (max, min) = ',maxval(by_r),minval(by_r)
  write(FILE_SYSOUT,*) ' bz      (max, min) = ',maxval(bz_r),minval(bz_r)

  write(FILE_SYSOUT,*) ' vx      (max, min) = ',maxval(vx_r),minval(vx_r)
  write(FILE_SYSOUT,*) ' vy      (max, min) = ',maxval(vy_r),minval(vy_r)
  write(FILE_SYSOUT,*) ' vz      (max, min) = ',maxval(vz_r),minval(vz_r)

  write(FILE_SYSOUT,*) ' ro      (max, min) = ',maxval(ro_r),minval(ro_r)
  write(FILE_SYSOUT,*) ' ph      (max, min) = ',maxval(ph_r),minval(ph_r)
  write(FILE_SYSOUT,*) ' time               = ',atime


! ========================================================================
! Analysis & Visualization
! ========================================================================
  call ana_1
  call phys_arrange
  call make_avs_field(cseries)

! ============================================================
!                   3D OUTPUT
! ===========================================================
  do i = 0,nx_r
  do j = 0,ny_r
     bx_p1(i,j,:) = bx_r(i,j,:)
     by_p1(i,j,:) = by_r(i,j,:)
     bz_p1(i,j,:) = bz_r(i,j,:)
     vx_p1(i,j,:) = vx_r(i,j,:)
     vy_p1(i,j,:) = vy_r(i,j,:)
     vz_p1(i,j,:) = vz_r(i,j,:)
     ro_p1(i,j,:) = ro_r(i,j,:)
     ph_p1(i,j,:) = ph_r(i,j,:)
  end do
  end do
!
  do i = 0,nx_r
     bx_p1(i,-1,:) = bx_r(i,0,:)
     by_p1(i,-1,:) = by_r(i,0,:)
     bz_p1(i,-1,:) = bz_r(i,0,:)
     vx_p1(i,-1,:) = vx_r(i,0,:)
     vy_p1(i,-1,:) = vy_r(i,0,:)
     vz_p1(i,-1,:) = vz_r(i,0,:)
     ro_p1(i,-1,:) = ro_r(i,0,:)
     ph_p1(i,-1,:) = ph_r(i,0,:)
  enddo

  do j = 0,ny_r
     bx_p1(-1,j,:) = bx_r(0,j,:)
     by_p1(-1,j,:) = by_r(0,j,:)
     bz_p1(-1,j,:) = bz_r(0,j,:)
     vx_p1(-1,j,:) = vx_r(0,j,:)
     vy_p1(-1,j,:) = vy_r(0,j,:)
     vz_p1(-1,j,:) = vz_r(0,j,:)
     ro_p1(-1,j,:) = ro_r(0,j,:)
     ph_p1(-1,j,:) = ph_r(0,j,:)
  enddo
    
! ===========================================================
!                    2D OUTPUT
! ===========================================================
  do i = 0,nx_r
  do j = 0,ny_r
     bx2d   (i,j) = bx_r(i,j,0)
     by2d   (i,j) = by_r(i,j,0)
     bz2d   (i,j) = bz_r(i,j,0)

     bx2d_P1(i,j) = bx_r(i,j,0)
     by2d_P1(i,j) = by_r(i,j,0)
     bz2d_P1(i,j) = bz_r(i,j,0)

     vx2d_P1(i,j) = vx_r(i,j,0)
     vy2d_P1(i,j) = vy_r(i,j,0)
     vz2d_P1(i,j) = vz_r(i,j,0) 
  end do
  end do

  do i = 0,nx_r
     bx2d_P1(i,-1) = bx2d(i,0)
     by2d_P1(i,-1) = by2d(i,0)
     bz2d_P1(i,-1) = bz2d(i,0)
!
     vx2d_P1(i,-1) = vx2d_P1(i,0)
     vy2d_P1(i,-1) = vy2d_P1(i,0)
     vz2d_P1(i,-1) = vz2d_P1(i,0)    
  enddo

  do j = 0,ny_r
     bx2d_P1(-1,j) = bx2d(0,j)
     by2d_P1(-1,j) = by2d(0,j)
     bz2d_P1(-1,j) = bz2d(0,j)
!
     vx2d_P1(-1,j) = vx2d_P1(0,j)
     vy2d_P1(-1,j) = vy2d_P1(0,j)
     vz2d_P1(-1,j) = vz2d_P1(0,j)     
  enddo

! ===========================================================
!               Modefied Magnetic Field
! ===========================================================
! do k = 0,nz_r-10
!    bx_rr(:,:,k) = bx_r(:,:,k)
!    by_rr(:,:,k) = by_r(:,:,k)
!    bz_rr(:,:,k) = bz_r(:,:,k)
! end do
! ===========================================================
!                   OutPut data
! ===========================================================
!   bx_sp(:,:,:) = real(bx_r(:,:,:))
!   by_sp(:,:,:) = real(by_r(:,:,:))
!   bz_sp(:,:,:) = real(bz_r(:,:,:))

!   vx_sp(:,:,:) = real(vx_r(:,:,:))
!   vy_sp(:,:,:) = real(vy_r(:,:,:))
!   vz_sp(:,:,:) = real(vz_r(:,:,:))

!! 3D data
!   open(11,file = dir//'B3D_MHD.'//trim(cloop),form='unformatted')
!   write(11)bx_r, by_r, bz_r
!   close(11) 
!
!   open(11,file = dir//'B3D_MHD_RS.'//trim(cloop),form='unformatted')
!   write(11)bx_rs, by_rs, bz_rs
!   close(11) 
 
! ===========================================================
! Debug
! ===========================================================
  open(11,file = dir//'BOTTOM_BZ.'//trim(cloop),form='formatted')
  do i = 0, nx_r
  do j=  0, ny_r
     write(11,*) xc_r(i),yc_r(j),bz_r(i,j,0)
  end do
     write(11,*)' '
  end do
  close(11)
 
  open(11,file = dir//'BOTTOM_BX.'//trim(cloop),form='formatted')
  do i = 0, nx_r
  do j=  0, ny_r
     write(11,*) xc_r(i),yc_r(j),bx_r(i,j,0)
  end do
     write(11,*)' '
  end do
  close(11)      

  open(11,file = dir//'BOTTOM_BY.'//trim(cloop),form='formatted')
  do i = 0, nx_r
  do j=  0, ny_r
     write(11,*) xc_r(i),yc_r(j),by_r(i,j,0)
  end do
     write(11,*)' '
  end do
  close(11)    
!
!  open(11,file = dir//'BOTTOM_VX.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),vx_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)   
!
!  open(11,file = dir//'BOTTOM_VY.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!    write(11,*) xc_r(i),yc_r(j),vy_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)
!
!  open(11,file = dir//'BOTTOM_VZ.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
! do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),vz_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)
!
  
!  open(11,file = dir//'BOTTOM_CX.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),cx_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)
 !
!  open(11,file = dir//'BOTTOM_CY.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),cy_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)
!
!  open(11,file = dir//'BOTTOM_CZ.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),cz_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
 ! close(11)
!
!  open(11,file = dir//'BOTTOM_FZ.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),cx_r(i,j,1)*by_r(i,j,1)-cy_r(i,j,1)*bx_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)  
!
! open(11,file = dir//'BOTTOM_FZ1.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),cx_r(i,j,1)*by_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)
!                                                        
!  open(11,file = dir//'BOTTOM_FZ2.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do j=  0, ny_r
!     write(11,*) xc_r(i),yc_r(j),-cy_r(i,j,1)*bx_r(i,j,1)
!  end do
!     write(11,*)' '
!  end do
!  close(11)  
!
!  open(11,file = dir//'VZRO'//'.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,vz_r(i,ny_r/2,k),ro_r(i,ny_r/2,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)
!
!  open(11,file = dir//'PH'//'.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,ph_r(i,ny_r/2,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'BOTTOM_CZ',form='formatted')
!  do i = 0, nx_r
!  do j = 0, ny_r
!     write(11,*) xc_r(i),yc_r(j),cz_r(i,j,0)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'VZCT_YZ'//'.'//trim(cloop),form='formatted')
!  do j = 0, ny_r
!  do k=  0, nz_r
!     write(11,*) j,k,vz_r(132,j,k),ct(132,j,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'VZ_XZ'//'.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,vz_r(i,165,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'RO_XZ'//'.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,ro_r(i,165,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)


!  open(11,file = dir//'BY_YZ'//'.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,by_r(i,0,k),by_r(i,ny_r,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(11,file = dir//'BZ_YZ'//'.'//trim(cloop),form='formatted')
!  do i = 0, nx_r
!  do k=  0, nz_r
!     write(11,*) i,k,bz_r(i,0,k),bz_r(i,ny_r,k)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

!  open(13,file = dir//'DTSTEP',form='formatted')
!  write(13,*) dtstep

  
!  open(11,file = dir//'BX_1D.'//trim(cloop),form='formatted')
!  do k = 0, nz_r
!     write(11,*) zc_r(k),bx_r(nx_r/2,ny_r/2,k)
!  end do
!  close(11)
!
!  open(11,file = dir//'BY_1D.'//trim(cloop),form='formatted')
!  do k = 0, nz_r
!     write(11,*) zc_r(k),by_r(nx_r/2,ny_r/2,k)
!  end do
!  close(11)  
!
!  open(11,file = dir//'BZ_1D.'//trim(cloop),form='formatted')
!  do k = 0, nz_r
!     write(11,*) zc_r(k),bz_r(nx_r/2,ny_r/2,k)
!  end do
!  close(11)
! 
!  open(11,file = dir//'VZ_1D.'//trim(cloop),form='formatted')
!  do k = 0, nz_r
!     write(11,*) zc_r(k),vz_r(nx_r/2,ny_r/2,k)
!  end do
!  close(11)  

  open(11,file = dir//'TIME.'//trim(cloop),form='formatted') 
     write(11,*) atime                                                                                   
  close(11) 
  
  end do

  stop
9999 write(*,*) '## file ',trim(filename),' not found ##'
  stop
  end program main










