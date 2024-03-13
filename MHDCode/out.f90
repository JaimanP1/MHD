  module out
! -----------------------------------------------------------
!                    MODULE out
! -----------------------------------------------------------
  use common
  use mpiut
  use mpi
  implicit none
  real(DP), dimension(2) :: eng, sendbuf, recvbuf

  integer,private :: i,j,k
  public
  private :: make_avs_field
  contains

! ============================================================
  subroutine out__disk
! ============================================================
  real(DP), dimension(0:NX,0:NY,0:NZ) :: array3d
!
! ----------------------
  iwrite = iwrite + 1   ! sequential data number
  cwrite = '_'//chari3(iwrite)
! ----------------------
! write 3d data (magnetic field)

      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'BX'//cwrite//cmyrank, &
       form='unformatted')
      call reform(bx,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)

      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'BY'//cwrite//cmyrank, &
       form='unformatted')
      call reform(by,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)

      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'BZ'//cwrite//cmyrank, &
       form='unformatted')
      call reform(bz,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)

! write 3d data (velocity) 

      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'VX'//cwrite//cmyrank, &
       form='unformatted')
      call reform(vx,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)

      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'VY'//cwrite//cmyrank, &
       form='unformatted')
      call reform(vy,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)

      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'VZ'//cwrite//cmyrank, &
       form='unformatted')
      call reform(vz,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)

! +++ DATA for diagnosis +++
      open(FILE_3D_FIELD, &
       file = trim(cfile_3d_field)//'DIVB'//cwrite//cmyrank, &
       form='unformatted')
      call reform(divb,array3d)
      write(FILE_3D_FIELD)  array3d
      close(FILE_3D_FIELD)


! ++++ AVS_Field_File

      call make_avs_field

  return
  end subroutine out__disk

! =========================================================== 
  subroutine out__restart
! ===========================================================
  iwrite = iwrite + 1   ! sequential data number
  cwrite = '_'//chari3(iwrite)

  open(FILE_RESTART,file = dir//trim(cfile_restart)//cwrite//cmyrank, &
       form='unformatted')
  write(FILE_RESTART) iwrite, nloop, atime, dtstep, bx, by, bz, ph
  close(FILE_RESTART)

  open(11,file = dir//'B3D_RESTART2'//cwrite//cmyrank,form='unformatted')
  write(11) vx,vy,vz,ro,ro_org
  close(11)

! -------------------
! Debug
! -------------------
!  do i = 0,nx
!  if(xc(i).eq.0) then
!     open(11,file = 'TEST_B'//cmyrank,form = 'formatted')
!     do j = 0,NY
!     do k = 0,NZ
!        write(11,*) yc(j),zc(k),by(k,j,i),bz(k,j,i)
!     end do
!        write(11,*)' '
!     end do
!     close(11)
!  end if
!  end do
!
!  do i = 0,nx
!  if(xc(i).eq.0) then
!     open(11,file = 'TEST_C'//cmyrank,form = 'formatted')
!     do j = 0,NY
!     do k = 0,NZ
!        write(11,*) yc(j),zc(k),cy(k,j,i),cz(k,j,i)
!     end do
!        write(11,*)' '
!     end do
!     close(11)
!  end if
!  end do
!

!  open(11,file = 'TEST_V'//cmyrank,form = 'formatted')
!  do i = 0,NX
!  do j = 0,NY
!     write(11,*) xc(i),yc(j),vx(0,j,i),vy(0,j,i)
!  end do
!     write(11,*)' '
!  end do
!  close(11)

   end subroutine out__restart


! ===========================================================
  subroutine make_avs_field
! ===========================================================
! === write 3D field file ====

      open(FILE_AVS_FIELD, &
      file = 'avs'//cwrite//cmyrank//'.fld',form='formatted')
      write(FILE_AVS_FIELD,'(a)') '# AVS field file for an active region'
      write(FILE_AVS_FIELD,'(a,i5)') 'ndim = ',3
      write(FILE_AVS_FIELD,'(a,i5)') 'dim1 = ',NX+1
      write(FILE_AVS_FIELD,'(a,i5)') 'dim2 = ',NY+1
      write(FILE_AVS_FIELD,'(a,i5)') 'dim3 = ',NZ+1
      write(FILE_AVS_FIELD,'(a,i5)') 'nspace = ',3
      write(FILE_AVS_FIELD,'(a,i5)') 'veclen = ',9
      write(FILE_AVS_FIELD,'(a)')    'data = double'
      write(FILE_AVS_FIELD,'(a)')    'field = rectilinear'
      write(FILE_AVS_FIELD,'(a)')    'label = b_x'
      write(FILE_AVS_FIELD,'(a)')    'label = b_y'
      write(FILE_AVS_FIELD,'(a)')    'label = b_z'
      write(FILE_AVS_FIELD,'(a)')    'label = cb2'
      write(FILE_AVS_FIELD,'(a)')    'label = divb'
      write(FILE_AVS_FIELD,'(a)')    'label = v_x'
      write(FILE_AVS_FIELD,'(a)')    'label = v_y'
      write(FILE_AVS_FIELD,'(a)')    'label = v_z'
      write(FILE_AVS_FIELD,'(a)')    'label = jxb_jb'
      write(FILE_AVS_FIELD,'(a)')    '#'
      write(FILE_AVS_FIELD,'(a)') 'coord 1 file='//trim(cfile_coordinate_x)//cmyrank//' filetype=ascii'
      write(FILE_AVS_FIELD,'(a)') 'coord 2 file='//trim(cfile_coordinate_y)//cmyrank//' filetype=ascii'
      write(FILE_AVS_FIELD,'(a)') 'coord 3 file='//trim(cfile_coordinate_z)//cmyrank//' filetype=ascii'
      write(FILE_AVS_FIELD,'(a)')    '#'
      write(FILE_AVS_FIELD,'(a)') 'variable 1 file='// &
                  trim(cfile_3d_field)//'BX'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 2 file='// &
                  trim(cfile_3d_field)//'BY'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 3 file='// &
                  trim(cfile_3d_field)//'BZ'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 4 file='// &
                  trim(cfile_3d_field)//'CB2'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 5 file='// &
                  trim(cfile_3d_field)//'DIVB'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 6 file='// &
                  trim(cfile_3d_field)//'VX'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 7 file='// &
                  trim(cfile_3d_field)//'VY'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 8 file='// &
                  trim(cfile_3d_field)//'VZ'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)') 'variable 9 file='// &
                  trim(cfile_3d_field)//'JXB_JB'//cwrite//cmyrank// &
                  ' filetype=unformatted'
      write(FILE_AVS_FIELD,'(a)')    '#'
      close(FILE_AVS_FIELD)

  end subroutine make_avs_field

! ============================================================
  subroutine reform(a,b)
! ============================================================
  real(DP), dimension(0:NZ,-1:NY,-1:NX) :: a
  real(DP), dimension(0:NX,0:NY,0:NZ) :: b
  integer :: i,j,k

  do i = 0, NX
  do j = 0, NY
  do k = 0, NZ      
     b(i,j,k) = a(k,j,i)
  end do
  end do
  end do

  end subroutine reform
  end module out






