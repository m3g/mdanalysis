!
! globalstructure.f90: Given a trajectory in DCD format and the corresponding
!                      PDB reference structure, computes the following
!                      structural factors:
!                      1) RMSD 
!                      2) Fraction of native contacts of the atoms,
!                         deffined as the number of atoms which are closer than
!                         some distance in each frame, relative to the
!                         number of contacts in the reference frame
!                         this contact distance is deffined by the 'contact' 
!                         keyword
!                      3) Gyration radius
!
!  Usage: ./globalstructure < globalstructure.inp
!
!  Input file example (it must have all keywords):
!
!                      pdbfile     ./ATRI_ready.pdb
!                      dcdfile     ./ATRI_desnat.dcd
!                      reference   DCD
!                      firstframe  1
!                      lastframe   113
!                      timestep    1.0
!                      scaletime   0.001
!                      firstatom   1
!                      lastatom    4278
!                      atomtype    CA
!                      segment     ALL
!                      contact     7.0
!                      rmsdout     ATRI_rmsd.dat
!                      natiout     ATRI_nativecontacts.dat
!                      gyraout     ATRI_gyration.dat
!
!  The reference option may be 'DCD', 'PDB' or 'AVG'. If 'PDB' is chosen, the
!  coordinates of the PDB file are chosen as the reference for alignment
!  and rmsd calculation. If 'DCD' is chosen, the coordinates of the
!  first frame of the DCD file are chosen. Finally, if 'AVG' is chosen,
!  the reference will be the average coordinates of the atoms in
!  the dcd file.
!
!  contact is the distance between atoms that are considered "in
!  contact". For CAs in proteins this may be about 7.0.
!
!  lastframe and lastatom may be set to LAST if all frames or atoms
!  are to be considered.
!
!  atomtype and segment may be set to ALL if all atoms or all segments
!  are to be considered.
! 
!  timestep is the timestep between to frames of the dcd file. It 
!  may have any unit, it is only used for the outputing the time
!  instead of the frame number in the output files
!  
!  scaletime will scale the x-axis output, for example for if the timestep
!  is 2 fs and you want the output in picoseconds
!
!  Changing maximum dimensions:
!     maxfrm: Maximum number of frames of the DCD file (1 change)
!     maxatm: Maximum number of atoms to be considered (2 changes)
!     maxcontact: Maximum number of contacts per atom (1 change)
! 
!  Leandro Martinez, IQ-UNICAMP, June 15, 2007
!
!  Version 16.323
!

program globalstructure

  implicit none
  integer, parameter :: maxfrm = 12000, maxatm = 10000, maxcontact = 1000
  
  character(len=4) :: dummyc, atomtype, reference, segment
  character(len=200) :: pdbfile, dcdfile, record, &
                keyword, lines(maxatm), rmsdout, natiout, &
                gyraout, value
  integer :: i, j, k, length, nset, latom, iatom, lframe, natom, &
          ntotat, index(maxatm), iframe, nframes, dummyi, &
          nc(maxatm), contact(maxatm,maxcontact), &
          nctot, npres, ifr, ioerr
  real :: x(maxfrm,maxatm), y(maxfrm,maxatm), z(maxfrm,maxatm), &
       aref(maxatm,3), t, dummyr, xvar(maxatm,3), &
       timestep, d, rmsd, xm, ym, zm, gyr, dcon, scaletime
  logical :: dcdunitcell, readfromdcd
  external :: version

! Input data from input file 

  scaletime = 1.
  do while(.true.) 
    read(*,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    if(keyword(record).eq.'pdbfile') then
      pdbfile = value(record)
    else if(keyword(record).eq.'dcdfile') then
      dcdfile = value(record)
    else if(keyword(record).eq.'reference') then
      reference = value(record)
    else if(keyword(record).eq.'scaletime') then
      record = value(record)
      read(record,*) scaletime
    else if(keyword(record).eq.'firstframe') then
      record = value(record)
      read(record,*) iframe
    else if(keyword(record).eq.'lastframe') then
      record = value(record)
      if(record.eq.'LAST') then
        lframe = 0
      else
        read(record,*) lframe
      end if
    else if(keyword(record).eq.'timestep') then
      record = value(record)
      read(record,*) timestep
    else if(keyword(record).eq.'firstatom') then
      record = value(record)
      read(record,*) iatom
    else if(keyword(record).eq.'lastatom') then
      record = value(record)
      if(record.eq.'LAST') then
        latom = 0
      else
        read(record,*) latom
      end if
    else if(keyword(record).eq.'contact') then
      record = value(record)
      read(record,*) dcon
    else if(keyword(record).eq.'atomtype') then
      record = value(record)
      read(record,*) atomtype
    else if(keyword(record).eq.'segment') then
      record = value(record)
      read(record,*) segment
    else if(keyword(record).eq.'rmsdout') then
      rmsdout = value(record)
    else if(keyword(record).eq.'natiout') then
      natiout = value(record)
    else if(keyword(record).eq.'gyraout') then
      gyraout = value(record)
    else if(keyword(record) .eq. "#") then
      cycle
    else if(keyword(record).gt.' ') then
      record = keyword(record)
      write(*,*) ' keyword: ', record(1:length(record)),' was not recognized. '
      stop
    end if
  end do

! Output input parameters for checking

  write(*,*)
  write(*,"( '  ',52('#') )") 
  write(*,*) ' GLOBALSTRUCTURE '
  write(*,*) ' Compute global structural factors from MD '
  write(*,"( '  ',52('#') )") 
  write(*,*)

  call version

  write(*,*) ' PDBFILE     = ', pdbfile(1:length(pdbfile))
  write(*,*) ' DCDFILE     = ', dcdfile(1:length(dcdfile))
  write(*,*) ' FIRST FRAME = ', iframe
  if(lframe.ne.0) then
    nframes = lframe - iframe + 1
    write(*,*) ' LAST FRAME  = ', lframe
    write(*,*) ' NUMBER OF FRAMES = ', nframes
  else 
    write(*,*) ' LAST FRAME  = LAST '
  end if
  write(*,*) ' TIME STEP   = ', timestep
  write(*,*) ' SCALETIME   = ', scaletime
  write(*,*) ' FIRST ATOM  = ', iatom
  if(latom.ne.0) then
    natom = latom - iatom + 1
    write(*,*) ' LAST ATOM   = ', latom
    write(*,*) ' NUMBER OF ATOMS IN RANGE  = ', natom
  else 
    write(*,*) ' LAST ATOM   = LAST '
  end if
  write(*,*) ' ATOM TYPE   = ', atomtype
  write(*,*) ' SEGMENT  = ', segment 
  write(*,*) ' CONTACT DISTANCE = ', dcon
  write(*,*) ' REFERENCE: ', reference
  write(*,*) ' RMSD FILE  = ', rmsdout(1:length(rmsdout))
  write(*,*) ' NATIVE CONTACTS FILE  = ', natiout(1:length(natiout))
  write(*,*) ' GYRATION RADIUS FILE  = ', gyraout(1:length(gyraout))
  write(*,*)

! Opening the pdb file and reading coordinates

  write(*,*) ' READING THE PDB FILE... '
  record = '###'
  open(10,file=pdbfile,status='old')
  do while(record(1:4).ne.'ATOM'.and.record(1:6).ne.'HETATM')
    read(10,"( a200 )") record
  end do
  if(iatom.gt.1) then
    do i = 1, iatom - 1 
      read(10,"( a200 )") record
    end do
  end if
  natom = 0
  i = iatom - 1
  do while(i.lt.latom.or.latom.eq.0)
    if(record(1:4).eq.'ATOM'.or.record(1:6).eq.'HETATM') then
      i = i + 1
      read(record(13:16),*) dummyc
      if(dummyc.eq.atomtype.or.atomtype.eq.'ALL') then
        if( segment .ne. 'ALL') then
          read(record(73:76),*,iostat=ioerr) dummyc
          if ( ioerr /= 0 ) then
            write(*,*) ' Error reading segment on PDB file. '
            write(*,*) ' line=',record(1:82)
            stop
          end if
        else 
          dummyc = '000'
        end if
        if(dummyc.eq.segment.or.segment.eq.'ALL') then
          natom = natom + 1
          index(natom) = i
          if(reference.eq.'PDB') read(record(31:54),*) (aref(natom,j), j = 1, 3)
          lines(natom) = record
        end if
      end if
    end if
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
  end do
  close(10)
  if ( natom .eq. 0 ) then
    write(*,*) ' ERROR: No atoms satisfy the selection criteria. ' 
    stop
  end if
  write(*,*) ' ACTUAL NUMBER OF ATOMS CONSIDERING TYPE = ', natom

  ! Opening the dcd file and reading the header

  readfromdcd = .true.
  call chkperiod(dcdfile,dcdunitcell,readfromdcd)

  open(10,file=dcdfile,status='old',form='unformatted')
  read(10) dummyc, nset, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
  read(10) dummyi, dummyr
  read(10) ntotat

  ! Checking if the DCD file contains unit cell information

  write(*,*) ' TOTAL NUMBER OF ATOMS IN DCD FILE  = ', ntotat
  write(*,*) ' TOTAL NUMBER OF FRAMES IN DCD FILE = ', nset

  if(lframe.gt.maxfrm) then
    write(*,*) ' ERROR: Last frame index greater than maximum number of frames possible. Increase maxfrm.'
    stop
  end if
  if(natom.gt.maxatm) then
    write(*,*) ' ERROR: Number of atoms greater than maximum number of frames possible. Increase maxatm.'
    stop
  end if
  if(lframe.gt.nset) then
    write(*,*) ' ERROR: Last frame grater than the number of frames in the dcd file. '
    stop      
  else if(lframe.eq.0) then
    lframe = nset
    nframes = lframe - iframe + 1
    write(*,*) ' DCD FRAMES TO BE READ = ', nframes
  end if

  !
  ! Reading the coordinates
  !

  ! Checking if the DCD file contains unit cell information

  write(*,*) ' READING THE DCD... '
  if(lframe.ge.50) then 
    write(*,"( '  0%',44(' '),'100%',/,'  ',$ )")
  else 
    write(*,"( '  ',$ )")
  end if
  do i = 1, lframe
    if(dcdunitcell) read(10)
    read(10,iostat=ioerr) (t, j=1,index(1)-1),x(i,1),((t,j=index(k-1)+1,index(k)-1),x(i,k),k=2,natom)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read coordinates from DCD file. '
      close(10)
      stop    
    end if
    read(10) (t, j=1,index(1)-1),y(i,1),((t,j=index(k-1)+1,index(k)-1),y(i,k),k=2,natom)
    read(10) (t, j=1,index(1)-1),z(i,1),((t,j=index(k-1)+1,index(k)-1),z(i,k),k=2,natom)
    if(mod(i,max0(1,lframe/50)).eq.0) write(*,"( a1,$ )") '*'
  end do
  write(*,*)
  close(10)

  ! Obtain the reference if it was chosen to be the first frame
      
  if(reference.eq.'DCD') then
    do i = 1, natom
      aref(i,1) = x(iframe,i)
      aref(i,2) = y(iframe,i)
      aref(i,3) = z(iframe,i) 
    end do
  end if

  ! Computing the average structure if it was set to be the reference

  if(reference.eq.'AVG') then
    write(*,*) ' Computing average structure... '
    do i = 1, natom
      do ifr = iframe, lframe
        aref(i,1) = aref(i,1) + x(ifr,i)
        aref(i,2) = aref(i,2) + y(ifr,i)
        aref(i,3) = aref(i,3) + z(ifr,i)
      end do
      aref(i,1) = aref(i,1) / dfloat(lframe-iframe)
      aref(i,2) = aref(i,2) / dfloat(lframe-iframe)
      aref(i,3) = aref(i,3) / dfloat(lframe-iframe)
    end do
  end if

  ! Aligning all frames to reference frame

  write(*,*) ' Aligning all frames to reference ... '
  do ifr = iframe, lframe
    do i = 1, natom
      xvar(i,1) = x(ifr,i)
      xvar(i,2) = y(ifr,i)
      xvar(i,3) = z(ifr,i)
    end do
    call align(natom,aref,xvar)
    do i = 1, natom
      x(ifr,i) = xvar(i,1)
      y(ifr,i) = xvar(i,2)
      z(ifr,i) = xvar(i,3)
    end do
  end do

  !
  ! Computing the rmsd relative to the reference position 
  !              

  write(*,*) ' Computing rmsd... '
  open(10,file=rmsdout)
  write(10,"( a )") '# Time(frame*timestep)    RMSD'
  do ifr = iframe, lframe
    rmsd = 0.d0
    do i = 1, natom
      rmsd = rmsd + ( x(ifr,i) - aref(i,1) )**2 &
                  + ( y(ifr,i) - aref(i,2) )**2 &  
                  + ( z(ifr,i) - aref(i,3) )**2   
    end do
    rmsd = sqrt(rmsd/float(natom))
    write(10,*) ifr*timestep*scaletime, rmsd
  end do
  close(10)

  !
  ! Computing the giration radius of the structure
  !

  write(*,*) ' Computing the gyration radius at each frame... '
  open(10,file=gyraout)
  write(10,"( a )") '# Time(frame*timestep)    Gyration radius'
  do ifr = iframe, lframe

  ! Computing the mean position of the particles

  xm = 0.d0
  ym = 0.d0
  zm = 0.d0
  do i = 1, natom 
    xm = xm + x(ifr,i)
    ym = ym + y(ifr,i)
    zm = zm + z(ifr,i)
  end do
  xm = xm / dfloat(natom)
  ym = ym / dfloat(natom)
  zm = zm / dfloat(natom)

  ! Computing the gyration radius
      
    gyr = 0.d0
    do i = 1, natom
      gyr = gyr + ( x(ifr,i) - xm )**2 + ( y(ifr,i) - ym )**2 + ( z(ifr,i) - zm )**2
    end do
    gyr = sqrt(gyr / float(natom))
    write(10,*) ifr*timestep*scaletime, gyr
  end do
  close(10)

  !
  ! Computing the fraction of native contacts
  !

  write(*,*) ' Computing fraction of native contacts... '
  do i = 1, natom
    nc(i) = 0
  end do
  do i = 1, natom
    do j = 1, natom
      if (abs(i-j).gt.3) then
        d = sqrt((aref(i,1)-aref(j,1))**2 + (aref(i,2)-aref(j,2))**2 + (aref(i,3)-aref(j,3))**2)
        if(d.lt.dcon) then
          nc(i) = nc(i) + 1
          if ( nc(i) .gt. maxcontact ) then 
            write(*,*) ' ERROR: Number of contacts > maxcontact '
            stop
          end if
          contact(i,nc(i)) = j
        end if
      end if
    end do
  end do

  ! Computing the total number of native contacts

  nctot = 0
  do i = 1, natom
    nctot = nctot + nc(i)
  end do

  ! Computing the percentage of native contacts preserved

  open(10,file=natiout)
  write(10,"( a )") '# Native contacts (NC) '
  write(10,"( a )") '# Time(frame*timestep)   NC(fraction) NC(number)'
  do ifr = iframe, lframe
    npres = 0
    do i = 1, natom
      do j = 1, nc(i)
        k = contact(i,nc(i))
        d = sqrt((x(ifr,i)-x(ifr,k))**2 + (y(ifr,i)-y(ifr,k))**2 + (z(ifr,i)-z(ifr,k))**2)
        if(d.le.dcon) npres = npres + 1 
      end do
    end do
    write(10,*) ifr*timestep*scaletime, float(npres)/float(nctot), npres
  end do
  close(10)

  write(*,*) ' FINISHED. '
  write(*,"( '  ',52('#') )") 
  write(*,*) 

end program globalstructure

!
! Subroutine align: given two sets of vectors, finds the best
!                   rotation that to align the two sets. Returns
!                   the second vector aligned to the first vector
!
!                 Method: S. K. Kearsley,
!                         "On the orthogonal transformation used for
!                          structural comparisons"
!                         Acta Cryst. (1989) A45, 208-210
!
!  On input: nat: the number of atoms of the structures
!            aref: the reference coordinates
!            x: the coordinates to be aligned to aref
!
!  On return:
!            x: the coordinates aligned
!
!  Requires: subroutine jacobi
!
!  Leandro Martinez, IQ-UNICAMP, 26/10/2005
!

subroutine align(nat,aref,x)

  implicit none
  integer, parameter :: maxatm = 10000
  integer :: i, j, k, nat, iq

  real :: x(maxatm,3), y(maxatm,3), cmx(3), cmy(3), q(4,4), xp(maxatm), &
          yp(maxatm), zp(maxatm), xm(maxatm), ym(maxatm), zm(maxatm), &
          a(4), u(3,3), xnew(maxatm,3), aref(maxatm,3)

  ! For ssyev
  real :: work(12)
  integer :: info

  ! Copying the reference atoms

  do i = 1, nat
    do j = 1, 3
      y(i,j) = aref(i,j)
    end do
  end do

  ! Computing the centroid of the structures

  do i = 1, 3
    cmx(i) = 0.
    cmy(i) = 0.
  end do

  do i = 1, nat
    do j = 1, 3
      cmx(j) = cmx(j) + x(i,j)
      cmy(j) = cmy(j) + y(i,j)
    end do
  end do

  do i = 1, 3
    cmx(i) = cmx(i) / float(nat)
    cmy(i) = cmy(i) / float(nat)
  end do

  ! Translating the atoms to the origin

  do i = 1, nat
    do j = 1, 3
      x(i,j) = x(i,j) - cmx(j)
      y(i,j) = y(i,j) - cmy(j)
    end do
  end do

  ! Computing the quaternion matrix

  do i = 1, nat
    xm(i) = y(i,1) - x(i,1)
    ym(i) = y(i,2) - x(i,2)
    zm(i) = y(i,3) - x(i,3)
    xp(i) = y(i,1) + x(i,1)
    yp(i) = y(i,2) + x(i,2)
    zp(i) = y(i,3) + x(i,3)
  end do

  do i = 1, 4
    do j = 1, 4
      q(i,j) = 0.
    end do
  end do

  do i = 1, nat
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)

  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix
  ! 'q' contains the eigenvectors associates with eigenvalues in ascending order

  call ssyev('V','U',4,q,4,a,work,12,info)

  ! Computing the rotation matrix

  iq = 1
  u(1,1) = q(1,iq)**2 + q(2,iq)**2 - q(3,iq)**2 - q(4,iq)**2
  u(1,2) = 2. * ( q(2,iq)*q(3,iq) + q(1,iq)*q(4,iq) )
  u(1,3) = 2. * ( q(2,iq)*q(4,iq) - q(1,iq)*q(3,iq) )
  u(2,1) = 2. * ( q(2,iq)*q(3,iq) - q(1,iq)*q(4,iq) )
  u(2,2) = q(1,iq)**2 + q(3,iq)**2 - q(2,iq)**2 - q(4,iq)**2
  u(2,3) = 2. * ( q(3,iq)*q(4,iq) + q(1,iq)*q(2,iq) )
  u(3,1) = 2. * ( q(2,iq)*q(4,iq) + q(1,iq)*q(3,iq) )
  u(3,2) = 2. * ( q(3,iq)*q(4,iq) - q(1,iq)*q(2,iq) )
  u(3,3) = q(1,iq)**2 + q(4,iq)**2 - q(2,iq)**2 - q(3,iq)**2

  ! Applying the rotation to the vectors in x

  do i = 1, nat
    do j = 1, 3
      xnew(i,j) = 0.
      do k = 1, 3
        xnew(i,j) = xnew(i,j) + u(j,k) * x(i,k)
      end do
    end do
  end do

  ! Translating to the centroid of the reference points

  do i = 1, nat
    do j = 1, 3
      xnew(i,j) = xnew(i,j) + cmy(j)
    end do
  end do

  ! Updating the 'x' vector

  do i = 1, nat
    do j = 1, 3
      x(i,j) = xnew(i,j)
    end do
  end do

  return
end subroutine align

