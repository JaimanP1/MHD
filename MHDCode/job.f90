  module job
! -----------------------------------------------------------
!                      MODULE job
! -----------------------------------------------------------
  use common
  use mpiut
  implicit none
  contains

! -----------------------------------------------------------
!
 function job__is_fine()
     logical :: job__is_fine
!
     if(nloop < nloop_end) then
        job__is_fine = .true.
     else
        job__is_fine = .false.
     end if
!
  end function job__is_fine

! -----------------------------------------------------------
  subroutine job__check_equal(index1, index2, message)
! -----------------------------------------------------------
  integer, INTENT(IN) :: index1, index2
  character(len=5), INTENT(IN) :: message
!
  if(index1 /= index2) then
     write(*,*) '##ERR:job__check_equal ',index1,index2,message
     call job__finalize
  end if
  return
  end subroutine job__check_equal

! -----------------------------------------------------------
  subroutine job__close_files
! -----------------------------------------------------------
  close(FILE_SYSOUT)
  close(FILE_COORDINATE_X)
  close(FILE_COORDINATE_Y)
  close(FILE_COORDINATE_Z)
  close(FILE_OUTPUT_LIST)
  if(myrank == root) close(FILE_TIME_LIST)

  return

  end subroutine job__close_files

! -----------------------------------------------------------
  subroutine job__finalize
! -----------------------------------------------------------
!
   call job__close_files
   call mpiut__finalize
!
   stop
 end subroutine job__finalize

end module job




