      subroutine r2f(ar,cf)
      use common
!
      real*8    ar(NX_W,NY_W)
      real*8    wkx(NY_W,2*NX_W),wky(NX_W,2*NY_W)
      complex*16 cf(NX_W,NY_W)

      do i = 1, NX_W
      do j = 1, NY_W
         wkx(j,2*i-1) = ar(i,j)/(NX_W*NY_W)
         wkx(j,2*i  ) = 0.0
      end do
      end do

      call dcft2(wkx,NY_W,NX_W,1)
 
      do i = 1, NX_W
      do j = 1, NY_W
         wky(i,2*j-1) = wkx(j,2*i-1)
         wky(i,2*j  ) = wkx(j,2*i  )
      end do
      end do

      call dcft2(wky,NX_W,NY_W,1)

      do i = 1, NX_W
      do j = 1, NY_W
         cf(i,j) = wky(i,2*j-1)+IUNIT*wky(i,2*j)
      end do
      end do

      return
! ========================
      entry f2r(cf,ar)
! ========================

      do i = 1, NX_W
      do j = 1, NY_W
         wkx(j,2*i-1) = real(cf(i,j))
         wkx(j,2*i  ) = dimag(cf(i,j))
      end do
      end do

      call dcft2(wkx,NY_W,NX_W,-1)

      do i = 1, NX_W
      do j = 1, NY_W
         wky(i,2*j-1) = wkx(j,2*i-1)
         wky(i,2*j  ) = wkx(j,2*i  )
      end do
      end do

      call dcft2(wky,NX_W,NY_W,-1)

      do i = 1, NX_W
      do j = 1, NY_W
         ar(i,j) = wky(i,2*j-1)
      end do
      end do

      return
      end
! ==========================
      subroutine dealia(cf)
! ==========================

      use common
      complex*16 cf(NX_W,NY_W)
!
      do i = 1, NX_W
         m = md(i)
      do j = 1, NY_W
         n = nd(j)

         if(abs(m).gt.NX_W/3.or.abs(n).gt.NY_W/3) then
            cf(i,j) = 0.0
         end if

      end do
      end do
          

      return
      end

         
      
