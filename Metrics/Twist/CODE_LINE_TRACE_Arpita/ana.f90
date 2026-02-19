! ==========================================================================
  subroutine ana(bz0,cz0,ct0,x0,y0,z0,xi,yi,zi,sbx,sby,sbz,bt,ct,ff_alpha0, &
                 ff_alpha, run_number)
! ==========================================================================
  implicit none
  include "common.f90"
  real*8,intent(inout) :: bz0,cz0,ct0,x0,y0,z0,xi,yi,zi,sbx,sby,sbz,bt,ct,  &
                          ff_alpha0,ff_alpha
  real*8               :: d_alpha, ave_alpha
  integer              :: run_number
  integer              :: i,j,k
  character*2          :: I_NUMBER

! --------------------------------------------------------------------------
! Output of data for the final value
! --------------------------------------------------------------------------
! magnetic twist Tw = (1/4.pi)*alpha.dl = (1/4.pi)*alpha_total                             
   open(75,file = dir//'TWIST_MAP',form='formatted')
   write(75,*) x0,y0,alpha_total/pi/4.
   
!   Open the TWIST_CHECK file once 
!   open(76, file = dir//'TWIST_CHECK',form='formatted',status='unknown',position ='append')

   ! If the Twist value is greater than 10, write to TWIST_CHECK file
!   if (alpha_total/pi/4. > 10.0d0) then
!      write(76,*) x0, y0, alpha_total/pi/4.
!   endif

!   close(76)            ! Close the TWIST_CHECK after writing 
   end subroutine ana
