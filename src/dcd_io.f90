!
! module dcd_io
!
! Contains four routines:
!
!   subroutine dcd_nframes: reads the number of frames of the dcd file.
!   subroutine dcd_natoms: reads the number of atoms of the dcd file.
!   subroutine read_dcd: reads a dcd file.
!   subroutine wrtie_dcd: writes a dcd file.
!
!   Details on each routine preceed the actual routine, below.
!
! Version 16.323.2
!

module dcd_io
  
  implicit none

  contains 

  !
  ! dcd_nframes: This subroutine reads the header of a dcd file and reads the
  !              number of frames
  !
  !   Input: name of the dcd file.
  !   Output: number of frames in this dcd file.
  !
  
  subroutine dcd_nframes(dcdfile,nframes)
  
    implicit none
    character(len=*), intent(in) :: dcdfile
    integer, intent(out) :: nframes
    integer :: ioerror
    character(len=4) :: dummyc
    
    open(10,file=dcdfile,action='read',status='old',form='unformatted',iostat=ioerror)
    if ( ioerror /= 0 ) then
      nframes = 0
      return
    end if
    read(10) dummyc, nframes
    close(10)
    
    return

  end subroutine dcd_nframes 
  
  !
  !  dcd_natoms: This subroutine reads the number of atoms of a dcd file from the
  !              header
  !
  !    Input: name of the dcd file.
  !    Output: number of atoms of the dcd file.
  !
  
  subroutine dcd_natoms(dcdfile,natoms)
  
    implicit none
    character(len=*), intent(in) :: dcdfile
    integer, intent(out) :: natoms
    integer :: ioerror, dummyi
    character(len=4) :: dummyc
    
    open(10,file=dcdfile,action='read',status='old',form='unformatted',iostat=ioerror)
    if ( ioerror /= 0 ) then
      natoms = 0
      return
    end if
    read(10) dummyc
    read(10) dummyi
    read(10) natoms
    close(10)
    
    return

  end subroutine dcd_natoms
  
  !
  ! read_dcd: This routine reads a dcd file and returns a vector
  !           containing the coordinates of the atoms for each frame.
  !           
  ! L. Martinez, Jun 29, 2011.
  !
  ! How to use it:
  !
  ! Calling the subroutine:
  !
  ! call read_dcd('/home/user/dcdfile.dcd',stride,firstframe,lastframe,firstatom,lastatom,
  !               unitcell,sides,x,y,z,error)
  !
  ! On input:
  !           dcdfile.dcd: Name of the dcd file to be read (with path).
  !           stride [integer]: return only one of each 'stride' frames.
  !                   for example, if stride=5, one of each 5 frames will be
  !                   returned.
  !           firstframe [integer]: Read from firstframe on (1 is first, not 0)
  !           lastframe [integer]: Read until the lastframe frame is reached.
  !           firstatom [integer]: index of the first atom os interest (1 is first, not 0)
  !           lastatom [integer]: index of the last atom of interest. 
  !           unitcell: logical variable that determines if the dcd file contains the
  !                     unitcell information or not
  ! 
  !     These vectors are empty on input, the dimensions are as specified:
  !           sides: Single precision vector of size sides(3*nframes), where
  !                  nframes=lastframe-firstframe+1. It will contain the sides
  !                  of the unit cell on each read frame, if "unitcell=.true.". If
  !                  "unitcell=.false." this vector is not used at all and, thus,
  !                  a dummy single precision variable can be passed.
  !           x,y,z: Single precision floating point vectors of dimension
  !                  x(nframes*natoms),y(nframes*natoms),z(nframes*natoms),
  !                  that will store the coordinates of each atom, on each frame.
  !                  where natoms=lastatom-firstatom+1 (the number of atoms for which
  !                  the coordinates will be read.
  !
  ! On return:
  !           sides: If "unitcell=.true.", this vector will contain the sides of the
  !                  unitcell organized as follows
  !                  sides(1) = x-size of the unit cell of frame 1
  !                  sides(2) = y-size of the unit cell of frame 1
  !                  sides(3) = z-size of the unit cell of frame 1
  !                  sides(4) = x-size of the unit cell of frame 2
  !                  etc. 
  !           x,y,z: The matrices above, filled with the coordinates of the atoms
  !                  at each point of the trajectory. Thus, for example, if
  !                  10 atoms were read:
  !                  x(1) is the 'x' coordinate of atom 1, of frame 1.
  !                  x(2) is the 'x' coordinate of atom 2, of frame 1.
  !                  etc.
  !                  x(11) is the 'x' coordiante of atom 2, of frame 2.
  !                  etc.
  !           error: Integer number with the flag of the error (or not), with the code:
  !                  error=0: No error.
  !                  error=1: Stride is less than 1.  
  !                  error=2: First frame is greater than last frame.
  !                  error=3: Firsta atom is greater than last atom.
  !                  error=4: Last atom required is greater than the number of atoms of the dcd file.
  !                  error=5: Last frame required is greater than the number of frames.
  !                  error=6: Error openning the the dcd file.
  !                  error=7: First frame or first atom are less than 1
  !
  
  subroutine read_dcd(filename,stride,firstframe,lastframe,firstatom,lastatom,unitcell,&
                      sides,x,y,z,error)
  
    implicit none
    
    ! Subroutine parameters and intents 
    
    logical, intent(in) :: unitcell
    character(len=*), intent(in) :: filename
    integer, intent(in) :: stride, firstframe, lastframe, firstatom, lastatom
    real, intent(out) :: x((lastframe-firstframe+1)*(lastatom-firstatom+1)),&
                         y((lastframe-firstframe+1)*(lastatom-firstatom+1)),&
                         z((lastframe-firstframe+1)*(lastatom-firstatom+1))
    real, intent(out) :: sides(3*(lastframe-firstframe+1))
    integer, intent(out) :: error

    ! Internal variables
    
    character(len=4) :: dummyc
    integer :: dummyi, index, i, j, ntotframes, ntotat, ioerror, iatom, natoms
    real :: dummyr
    double precision :: readsidesx, readsidesy, readsidesz, dummyd

    ! Checking for simple input errors

    error = 0 
    if(stride < 1) then
      error = 1
      return
    end if
    if(lastframe < firstframe) then
      error = 2
      return
    end if
    if(lastatom < firstatom) then
      error = 3
      return
    end if
    if(firstframe < 1 .or. firstatom < 1) then
      error = 7
      return
    end if
    
    ! Reading the header of the DCD file
    
    open(10,file=filename,action='read',form='unformatted',iostat=ioerror)
    if ( ioerror /= 0 ) then
      error = 6
      return
    end if
    read(10) dummyc, ntotframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(10) dummyi, dummyr
    read(10) ntotat
    
    if ( lastatom > ntotat ) then
      error = 4
      return
    end if
    if ( lastframe > ntotframes ) then
      error = 5
      return
    end if

    ! Number of atoms that will be read per frame

    natoms = lastatom - firstatom + 1

    ! Reading the coordinates
    
    index = 0
    do i = 1, lastframe
    
    ! If the first frame was not reached or stride says so, just cycle
    
      if ( i < firstframe .or. mod(i-firstframe,stride) /= 0 ) then 
        if(unitcell) read(10) dummyd
        read(10) dummyr
        read(10) dummyr
        read(10) dummyr
        cycle
      end if
    
    ! Read the sizes of the box
    
      index = index + 1
      j = 3*(index-1)
      if(unitcell) then
        read(10) readsidesx, dummyd, readsidesy, dummyd, dummyd, readsidesz
        sides(j+1) = sngl(readsidesx)
        sides(j+2) = sngl(readsidesy)
        sides(j+3) = sngl(readsidesz)
      end if
    
    ! Reading coordinates
      
      iatom = natoms*(index - 1)                                           
      read(10) (dummyr, j = 1, firstatom - 1 ), (x(j), j = iatom + 1, iatom + natoms)
      read(10) (dummyr, j = 1, firstatom - 1 ), (y(j), j = iatom + 1, iatom + natoms)            
      read(10) (dummyr, j = 1, firstatom - 1 ), (z(j), j = iatom + 1, iatom + natoms)           
    
    end do
    close(10)

    return

  end subroutine read_dcd
  
  
  !
  ! write_dcd: This routine writes a dcd given the coordinate vectors
  !            and number of atoms
  !           
  ! L. Martinez, Jul 4, 2011.
  !
  ! How to use it:
  !
  ! Calling the subroutine:
  !
  ! call write_dcd('dcdfile.dcd',natoms,nframes,x,y,z,unitcell,sides,error)
  !
  ! On input:
  !           dcdfile.dcd: Name of the dcd file to be written (with path).
  !           natoms [integer]: number of atoms
  !           nframes [integer]: number of frames
  !           x,y,z: Single precision floating point vectors of dimension
  !                  x(nframes*natoms),y(nframes*natoms),z(nframes*natoms),
  !                  that store the coordinates of each atom, on each frame.
  !           unitcell: Logical variable which must be set to "true" if the sides
  !                     of the unit cell per frame will be provided in the sides()
  !                     vector, or "false" if not.
  !           sides: Single precision vector of dimension sides(3*nframes), containing
  !                  the sides of the unit cell in each frame, as:
  !                  sides(1) = x-size of the unit cell of frame 1
  !                  sides(2) = y-size of the unit cell of frame 1
  !                  sides(3) = z-size of the unit cell of frame 1
  !                  sides(4) = x-size of the unit cell of frame 2
  !                  etc.
  !           error: Integer indicating the status of the output of the routine.
  !                  If no error occurred, error=0. Else:
  !                  error=1: Number of atoms was set to less than 1. 
  !                  error=2: Number of frames was set to less than 1. 
  !                  error=3: Error opening dcd file for writting.
  !
  ! On return: nothing. The routine will create (and overwrite if
  !            existent!) the DCD file indicated.
  !
  
  subroutine write_dcd(filename,natoms,nframes,x,y,z,unitcell,sides,error)
  
    implicit none
    
    ! Subroutine parameters and intents 
    
    logical, intent(in) :: unitcell
    character(len=*), intent(in) :: filename
    integer, intent(in) :: natoms, nframes
    real, intent(in) :: x(natoms*nframes),&
                        y(natoms*nframes),&
                        z(natoms*nframes)
    double precision :: readsidesx, readsidesy, readsidesz
    real :: sides(3*nframes)
    integer, intent(out) :: error
    
    ! Internal variables
    
    integer :: i, j, ioerror, iatom
    character(len=80), parameter :: title='DCDfile'
    
    ! Checking some data
    
    error = 0
    if ( natoms < 1 ) then
      error = 1
      return
    end if
    if ( nframes < 1 ) then
      error = 2
      return
    end if
    
    ! Writting the header of the DCD file
    
    open(10,file=filename,action='write',status='replace',form='unformatted',iostat=ioerror)
    if ( ioerror /= 0 ) then
      error = 3
      return
    end if
    
    write(10) 'CORD',nframes,1,1,nframes,0,0,0,3*natoms,0,0.001,1,0,0,0,0,0,0,0,0,24
    write(10) 1, title
    write(10) natoms
    
    ! Writting the coordinates
    
    do i = 1, nframes
    
    ! Read the sizes of the box
    
      j = 3*(i-1)
      if ( unitcell ) then 
        readsidesx = dble(sides(j+1))
        readsidesy = dble(sides(j+2))
        readsidesz = dble(sides(j+3))
        write(10) readsidesx, 0.d0, readsidesy, 0.d0, 0.d0, readsidesz
      end if
    
    ! Writting the coordinates
      
      iatom = natoms*(i-1)                                          
      write(10) (x(j), j = iatom + 1, iatom + natoms)
      write(10) (y(j), j = iatom + 1, iatom + natoms)            
      write(10) (z(j), j = iatom + 1, iatom + natoms)           
    
    end do
    
    close(10)
    
    return

  end subroutine write_dcd

end module dcd_io
