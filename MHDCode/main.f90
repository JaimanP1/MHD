! ===========================================================
!              Magnetohydrodynamic code Ver1 
! ===========================================================
  program main
! -----------------------------------------------------------
!                      PROGRAM MAIN
! -----------------------------------------------------------
  use common
  use job
  use mpi
  use mpiut
  use pset
  use iset
  use rkg
  use mhd
  use out
  implicit none
  integer :: istep, iset_err

! ==========================================================
! read NAMELIST 
! ==========================================================
  call pset__namelist

! ==========================================================
! MPI init  
! ==========================================================
  call mpiut__init
  cmyrank = '.'//chari4(myrank)

! ==========================================================
! set Initial conditions
! ==========================================================
  call pset__init 
  call iset__initial(iset_err) 

  if(iset_err /= 0) go to 900

  call mhd__sub
  call rkg__set_coefficient 

  if(start_from_initial) then
     call out__restart
  end if

! =============================================================================
! MAIN LOOP 
! =============================================================================
  do while(job__is_fine())
     nloop = nloop + 1

! -----------------------------------------------------------------------------
! set boundary condition
! -----------------------------------------------------------------------------
  call mhd__bc_velocity

! -----------------------------------------------------------------------------
! Mag E. and Kin. E    
! -----------------------------------------------------------------------------
  call  mhd__set_integral

! -----------------------------------------------------------------------------
! Emergency Stop
! -----------------------------------------------------------------------------
  if(ave_ene_mag.gt.1.0e+00) then
  go to 90
  end if     

! -----------------------------------------------------------------------------
! RKG
! -----------------------------------------------------------------------------
  do istep = 1, 4

!  call mhd__fluxemergence
  call mhd__sub

   if(istep.eq.1) then
      call pset__dt
   end if

   call mhd__main
   if(istep.eq.1) call rkg__prog1
   if(istep.eq.2) call rkg__prog2
   if(istep.eq.3) call rkg__prog3
   if(istep.eq.4) call rkg__prog4

  end do
  
! =============================================================================
! OUTPUT
! =============================================================================
90 if(mod(nloop,nloop_output) == 0) call out__restart

  end do

! =================== close files ============================
900 call job__finalize

  stop
  end program main









