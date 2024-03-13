c***********************************************************************
      subroutine dcft2(data,ll,nn,isign)
c***********************************************************************
c
c this is a (1-d) fft routine which transforms, simultaneously, ll-set
c  of nn-dimensional double precision complex data 'data'.
c
c nn must be a power of 2. this is not checked for.
c
c vectorization is taken in the ll-loop: so this routine is effective
c  when you have quite many data with same dimension to be transformed,
c  i.e., ll>>nn.
c
c        data( 1,real(1))
c        data( 2,real(1))
c              !
c        data(ll,real(1))
c        data( 1,imag(1))          do 1  n = 1 , nn
c              !                     do 2  l = 1 , ll
c        data(ll,imag(1))              data(l,2*n-1) = re(n-th data)
c        data( 1,real(2))              data(l,2*n)   = im(n-th data)
c              !                 2   continue
c        data(ll,real(2))        1 continue
c        data( 1,imag(2))
c              !
c        data(ll,imag(2))
c        data( 1,real(3))
c              !
c              !
c              !
c        data( 1,real(nn))
c              !
c        data(ll,real(nn))
c        data( 1,imag(nn))
c              !
c        data(ll,imag(nn))
c
c data is replaced by nn times of its fourier transform (isign=1)
c  or is replaced by its inverse fourier transform     (isign=-1).
c so you must devide data by nn, if necessary.
c
c 1-d fft used in this routine is based on danielson-lanczos algorithm.
c referece:  "numerical recipe in c", press et al.
c  cambridge university press, 1988, p412
c
c  by kageyama, a    1992.04.04
c
      implicit real*8 (a-h,o-z)
c
      parameter(twopi=6.2831853071795864d0)
c
      real*8 data(ll,2*nn)
c
      n2 = 2*nn
      j = 1
c
      do 1010  i = 1 , n2 , 2
        if(j.gt.i) then
          do 1500  l = 1 , ll
                tempr   = data(l,j)
                tempi   = data(l,j+1)
            data(l,j)   = data(l,i)
            data(l,j+1) = data(l,i+1)
            data(l,i)   = tempr
            data(l,i+1) = tempi
 1500     continue
        endif
        m = nn
        do 1000  while( (m.ge.2).and.(j.gt.m) )
          j = j-m
          m = m/2
 1000   continue
        j = j + m
 1010 continue
c
      mmax = 2
c
      do 2000  while( n2.gt.mmax )
        istep = 2*mmax
        theta = -twopi/dfloat(isign*mmax)
        wpr   = -2.d0*dsin(0.5d0*theta)**2
        wpi   = dsin(theta)
        wr    = 1.d0
        wi    = 0.d0
        do 2300 m = 1 , mmax , 2
          do 2200 i = m , n2 , istep
            j = i + mmax
            do 2100 l = 1 , ll
              tempr = wr*data(l,j)-wi*data(l,j+1)
              tempi = wr*data(l,j+1)+wi*data(l,j)
              data(l,j)   = data(l,i) - tempr
              data(l,j+1) = data(l,i+1) - tempi
              data(l,i)   = data(l,i) + tempr
              data(l,i+1) = data(l,i+1) + tempi
 2100       continue
 2200     continue
          wtemp = wr
          wr = wr*wpr -    wi*wpi + wr
          wi = wi*wpr + wtemp*wpi + wi
 2300   continue
        mmax = istep
 2000 continue
c
      return
      end







