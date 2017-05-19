!
! Version 16.323.2
!

!
! Controls the size of the character variables
!
module charsize
  integer, parameter :: charsize1 = 6
end module charsize

! Function that prints the version

subroutine version()
  
  write(*,*) 
  write(*,*) ' Version 17.138 '
  write(*,*) 

end subroutine version

!
! This file contains common functions and subroutines for
! all programs. Must be compiled with each of them.
!
! L. Martinez, I. Pasteur, Jun 18, 2008.
!

!
! Gets keyword from input file
!

function keyword(string)

  implicit none
  integer :: if, il
  character(len=200) :: keyword, string
  
  if = 1
  do while(string(if:if) <= ' '.and. if < 200) 
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do
  il = il - 1
  keyword = string(if:il)

return
end function keyword    

!
! Gets keyword value from input file
!

function value(string)

  implicit none
  integer :: if, il, length
  character(len=200) :: value, string
  
  if = 1
  do while(string(if:if) <= ' '.and.if < 200) 
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do
  il = il - 1
  if = il + 1 
  do while(string(if:if) <= ' '.and.if < 200) 
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do
  value = string(if:il)
  if(length(value) == 0) then
    write(*,*) ' ERROR: Some keyword without value: '
    write(*,*) string(1:length(string))
    stop
  end if

return
end function value

!
! subroutine getdim: Simply return the number of atoms of the system
!                    as specified in the psf and the number of
!                    classes in parameter files to allocate arrays.
!
! On input: psffile: name of the psf file
!           inputfile: name of the input file
!
! On return: ndim: minimum required dimension
!

subroutine getdim(psffile,inputfile,ndim)   

  implicit none
  integer :: ndim, nclass, status
  character(len=200) :: psffile, inputfile, record, keyword,&
                        value, file
  character :: firstchar
  
  open(10,file=psffile,action='read',iostat=status)
  if ( status /= 0 ) then
    write(*,*) ' ERROR: Could not open PSF file: ', trim(adjustl(psffile))
    stop
  end if
  read(10,"( a200 )") record
  do while(record /= '!NATOM')
    read(10,*,iostat=status) ndim, record
  end do
  close(10)
  
  nclass = 0
  open(99,file=inputfile,action='read')
  do while(.true.)
    read(99,"( a200 )",iostat=status) record
    if(status /= 0) exit
    if(keyword(record) == 'par') then
      file = value(record)
      record(1:9) = '#########'
      open(10,file=file,action='read',status='old',iostat=status)
      if(status /= 0) then
        write(*,*) ' ERROR: Parameter file not found: '
        write(*,*) trim(file)
        stop
      end if
      do while (record(1:9) /= 'NONBONDED')
        read(10,"( a200 )",iostat=status) record   
        if(status /= 0) then
          write(*,*) ' ERROR: Error reading parameter file (after NONBONDED mark): '
          write(*,*) trim(file)
          stop
        end if
      end do
      do while(.true.)
        read(10,"( a200 )",iostat=status) record
        if(status /= 0) exit
        if(firstchar(record) /= '!'.and.firstchar(record) > ' ') then
          if(status /= 0) cycle
          nclass = nclass + 1
        end if
      end do
      close(10)
    end if
  end do
  close(99)  
  
  write(*,*) ' Number of atoms of the system: ', ndim
  write(*,*) ' Number of classes in parameter files: ', nclass
  if(nclass > ndim) ndim = nclass

return
end subroutine getdim 
                           
!
! subroutine readpsf: Reads a PSF file to get charges and assign
!                     Lennard-Jones parameters that are provided 
!                     as input (which were read previously using
!                     subroutine readpar). In this program the
!                     only used parameter are the masses, actually.
!
! On input: psffile: name of the psf file
!           nclass: Number of atom classes read before by readpar
!           class: array containing the classes of each atom read
!                  by readpar      
!           eps: array containing the epsilon parameter of each type
!                of atom in each residue
!           sig: array containing the sigma parameter of each type
!                of atom in each residue
!           assignpars: .true. if the parameter file was read before
!                      and the parameters are needed, .false. 
!                      otherwise.
!
! On return: natom: total number of atoms of the PSF file
!            segat: name of the segment of each atom
!            resat: name of the residue of each atom 
!            resid: index of the residue of the atom
!            classat: class of each atom
!            typeat: type of each atom
!            q: array containing the charge of each atom
!            e: array containing the epsilon parameter for each atom
!            s: array containing the sigma parameter for each atom
!            mass: array containing the mass of each atom
!

subroutine readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
                   resid,classat,typeat,q,e,s,mass,assignpars)   

  use charsize
  implicit none
  integer :: i, j, natom, length, nclass, resid(*), status, iresid,&
             irlast
  real :: eps(*), sig(*), q(*), e(*), s(*), mass(*)
  character(len=charsize1) :: segat(*), resat(*), classat(*),& 
                              typeat(*), class(*)
  character(len=200) :: psffile, record
  logical :: found, assignpars
  
  ! Reading PSF file
  
  open(10,file=psffile,action='read')
  read(10,*,iostat=status) natom, record
  do while( record /= '!NATOM')
    read(10,*,iostat=status) natom, record
  end do
  
  resid(1) = 1
  do i = 1, natom
    read(10,"( a200 )") record
    read(record,*,iostat=status) j, segat(i), iresid, resat(i),&
                                 typeat(i), classat(i), q(i), mass(i)
    if ( status /= 0 ) then
      write(*,*) ' ERROR: Reading atom line in psf file: '
      write(*,*) record(1:length(record))
      write(*,*) ' Expected: number, segment, residue number, type,',&
                 ' class, charge, mass '
      stop
    end if
    if(i.gt.1) then
      if(segat(i) == segat(i-1) .and. &
         resat(i) == resat(i-1) .and. &
         iresid == irlast ) then
        resid(i) = resid(i-1)
      else
        resid(i) = resid(i-1) + 1
      end if   
    end if
    irlast = iresid
  end do
  close(10)
  
  ! Assigning paramaters for each atom
  
  if(.not.assignpars) return
  
  do i = 1, natom
    j = 0                             
    found = .false.
    do while(.not.found.and.j < nclass)
      j = j + 1
      if(class(j) == classat(i)) then
        e(i) = eps(j)
        s(i) = sig(j)
        found = .true.
      end if
    end do   
    if(.not.found) then
      write(*,*) ' ERROR: Could not find Lennard-Jones parameters for atom:',&
                   segat(i),' ',resat(i),' ',typeat(i),' ',classat(i)
      stop
    end if
  end do
  
return
end subroutine readpsf

!
! Subroutine that checks if dcd file contains periodic cell information
!

subroutine chkperiod(dcdfile,dcdaxis,readfromdcd) 

  implicit none
  real :: dummyr, x
  double precision :: side(6)
  integer :: dummyi, ntotat, i, status
  character(len=200) :: dcdfile
  character(len=4) :: dummyc
  logical :: dcdaxis, readfromdcd
  
  dcdaxis = .false.
  open(10,file=dcdfile,action='read',form='unformatted',status='old',iostat=status)
  if( status /= 0 ) then
    write(*,*) ' ERROR: Error opening dcd file: ' 
    write(*,*) trim(dcdfile)
    stop
  end if
  read(10,iostat=status) dummyc, dummyi, (dummyi,i=1,8), dummyr,&
                         (dummyi,i=1,9)
  if(status /= 0) then
    write(*,*) ' ERROR: Error reading dcd file: '
    write(*,*) trim(dcdfile)
    stop
  end if
  read(10,iostat=status) dummyi, dummyr
  if(status /= 0) then
    write(*,*) ' ERROR: Error reading dcd file: '
    write(*,*) trim(dcdfile)
    stop
  end if
  read(10,iostat=status) ntotat
  if(status /= 0) then
    write(*,*) ' ERROR: Error reading dcd file: '
    write(*,*) trim(dcdfile)
    stop
  end if
  read(10,iostat=status) (x,i=1,ntotat)
  close(10)
  
  ! If periodic cell information was found:
  
  if(status /= 0) then
    open(10,file=dcdfile,action='read',form='unformatted')
    read(10) dummyc, dummyi, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
    read(10) dummyi, dummyr
    read(10) ntotat
    read(10,iostat=status) (side(i),i=1,6)
    close(10)
    if( status /= 0 ) then
      write(*,"(' ERROR: Could not read either coordinates nor ',/,&
               &'        periodic cell sizes in first line of ',/,&
               &'        first frame of DCD file.' )")
      stop
    end if
    dcdaxis = .true.
    write(*,*) ' DCD file appears to contain periodic cell information. '
    write(*,"( '  Sides in first frame: ',3(f8.2) )") side(1), side(3), side(6)
    if(.not.readfromdcd) then
      write(*,*) ' Warning: Calculation will not use this information! '
      write(*,*) ' To use it set: periodic readfromdcd'
    end if
    return 
  
  ! If periodic cell information was not found:
  
  else
    dcdaxis = .false.
    write(*,*) ' DCD file does not contain periodic cell information. '
    if(readfromdcd) then
      write(*,*) ' ERROR: Periodic cell information not found in dcd&
                   & file and periodic=readfromdcd '
      stop
    end if
  end if

return
end subroutine chkperiod

!
! Subroutine image: computes the minimum image given the coordinates
!                   and the size of the box. Modifies the coordinates
!                   given
!

subroutine image(x,y,z,sidex,sidey,sidez)

  implicit none
  real :: x, y, z, sidex, sidey, sidez
  
  if(sidex < 1.e-5 .or. sidey < 1.e-5 .or. sidez < 1.e-5) then
    write(*,*) ' ERROR: Periodic box side too short. '
    stop
  end if
  x = mod(x,sidex)
  y = mod(y,sidey)
  z = mod(z,sidez)
  if(x > sidex/2) x = x - sidex
  if(y > sidey/2) y = y - sidey
  if(z > sidez/2) z = z - sidez
  if(x < -sidex/2) x = x + sidex
  if(y < -sidey/2) y = y + sidey
  if(z < -sidez/2) z = z + sidez
  
return
end subroutine image


!
! Subroutine readpar: Reads the parameter files in charmm format
!                     and assign for each class of atom its Lennard-Jones
!                     parameters
!
! On input: parfile: name of the parameter file to be read
!           nclass: number of classes of atom read including previous
!                   parameter files
!           class: array containing the class of each atom of each
!                  residue up to nclass
!
! On return: eps: array containing the epsilon parameter class of
!                 atom
!            sig: array containing the sigma parameter class of atom
!            nclass: number of classes of atoms including the current
!                    file
!            class: array containing the atom classes including
!                   the current classes 
!

subroutine readpar(parfile,nclass,class,eps,sig)

  use charsize
  implicit none
  integer :: nclass, status
  real :: eps(*), sig(*)
  real :: dummye
  character :: firstchar
  character(len=charsize1) :: class(*)
  character(len=200) :: record, parfile
  
  ! Reading the parameter file
  
  record(1:9) = '#########'
  open(10,file=parfile,action='read')
  do while (record(1:9) /= 'NONBONDED')
    read(10,"( a200 )") record   
  end do
  do while(.true.)
    read(10,"( a200 )",iostat=status) record
    if(status /= 0) exit
    if(firstchar(record) /= '!'.and.firstchar(record) > ' ') then
      read(record,*,iostat=status) class(nclass+1), dummye,&
                                   eps(nclass+1), sig(nclass+1)
      if(status /= 0) cycle
      nclass = nclass + 1
    end if
  end do
  close(10)

return
end subroutine readpar

!
! Get the first non-blanck character from a line
!

function firstchar(record)

  implicit none
  character(len=200) :: record
  character :: firstchar
  integer :: i

  i = 1
  do while(record(i:i) <= ' ')
    i = i + 1
    if ( i == 200 ) exit
  end do
  firstchar = record(i:i)

return
end function firstchar

!
! Function that sets the length of a string
!

function length(string)

  implicit none
  integer :: length
  character(len=200) :: string
  
  length = 200
  do while(string(length:length) <= ' ')
    length = length - 1
    if ( length == 0 ) exit
  end do

return
end function length      
 
!
! Compute the norm of a vector
!

function xnorm(x1,x2,x3)

  implicit none
  real :: xnorm, x1, x2, x3
  
  xnorm = sqrt(x1*x1 + x2*x2 + x3*x3)

return 
end function xnorm

!
! Given two vectors returns the cosine of the
! angle between them 
!

function cosang(x1,y1,z1,x2,y2,z2)

  implicit none
  real :: cosang, x1, y1, z1, x2, y2, z2, xnorm
  
  cosang = x1*x2 + y1*y2 + z1*z2
  cosang = cosang / ( xnorm(x1,y1,z1) * xnorm(x2,y2,z2) )

return
end function cosang

! 
! Subroutine that computes the center of mass of a group
! of atoms from a sequential vector of coordinates.
! 
! On input:
!
! n_atoms: Number of atoms of the group.
! i_atoms: Vector with the atom indices in the structure.
! mass: Vector with atom masses of all atoms in the structure.
! mass_group: Total mass of the group.
! x_at, y_at, z_at: Current dcd coordinates.
! 
! On return:
!
! cmx, cmy, cmz: The center of mass of the group.
!

subroutine compute_cm(n_atoms,i_atoms,mass,mass_group,&
                      x_at,y_at,z_at,cmx,cmy,cmz)

  implicit none
  integer :: i, i_atoms(*), n_atoms
  real :: x_at(*), y_at(*), z_at(*)
  real :: mass(*), mass_group, cmx, cmy, cmz
  
  cmx = 0.d0
  cmy = 0.d0
  cmz = 0.d0
  do i = 1, n_atoms
    cmx = cmx + x_at(i)*mass(i_atoms(i))
    cmy = cmy + y_at(i)*mass(i_atoms(i))
    cmz = cmz + z_at(i)*mass(i_atoms(i))
  end do
  cmx = cmx / mass_group
  cmy = cmy / mass_group
  cmz = cmz / mass_group 
  
return
end subroutine compute_cm 

! 
! Subroutine that computes the center of mass of a group
! of atoms from a frame of a dcd file.
! 
! On input:
!
! n_atoms: Number of atoms of the group.
! i_atoms: Vector with the atom indices in the structure.
! mass: Vector with atom masses of all atoms in the structure.
! mass_group: Total mass of the group.
! xdcd, ydcd, zdcd: Current dcd coordinates.
! iatom: Index of the first atom on this frame - 1
! 
! On return:
!
! cmx, cmy, cmz: The center of mass of the group.
!

subroutine compute_cm_dcd(n_atoms,i_atoms,mass,mass_group,&
                          xdcd,ydcd,zdcd,iatom,&
                          cmx,cmy,cmz)

  implicit none
  integer :: i, ii, n_atoms, i_atoms(*), iatom
  real :: mass(*), mass_group, cmx, cmy, cmz
  real :: xdcd(*), ydcd(*), zdcd(*)
  
  cmx = 0.d0
  cmy = 0.d0
  cmz = 0.d0
  do i = 1, n_atoms
    ii = iatom + i_atoms(i)
    cmx = cmx + xdcd(ii)*mass(i_atoms(i))
    cmy = cmy + ydcd(ii)*mass(i_atoms(i))
    cmz = cmz + zdcd(ii)*mass(i_atoms(i))
  end do
  cmx = cmx / mass_group
  cmy = cmy / mass_group
  cmz = cmz / mass_group 
  
return
end subroutine compute_cm_dcd 

! Subroutine that corrects the coordinates of the atoms in phantom
! cells for computation

subroutine movephantomcoor(xin,yin,zin,xout,yout,zout,&
                           nboxes,ibox,jbox,kbox,axis)

  integer :: nboxes(3), ibox, jbox, kbox
  real :: xin, yin, zin, xout, yout, zout, axis(3)

  xout = xin
  yout = yin
  zout = zin

  if ( ibox == 0 ) then
    xout = xin - axis(1)
  else if ( ibox == nboxes(1)+1 ) then
    xout = xin + axis(1) 
  else 
    xout = xin
  end if

  if ( jbox == 0 ) then
    yout = yin - axis(2)
  else if ( jbox == nboxes(2)+1 ) then
    yout = yin + axis(2)
  else
    yout = yin
  end if 

  if ( kbox == 0 ) then 
    zout = zin - axis(3)
  else if ( kbox == nboxes(3)+1 ) then
    zout = zin + axis(3)
  else
    zout = zin 
  end if

  return
end subroutine movephantomcoor

! 
! Subroutine that fills the phatnom boxes for periodic boundary conditions
! with the linked cell method
!
! On input: nboxes(3): is the vector containing the number of cells in each
!                     direction (except phantom ones)
!           iatomfirst: Is the first atom of each cell already filled up for
!                       real cells. 
!
! On output: iatomfirst will be filled up completely in the borders with replicas
!            of the corresponding iatomfirst parameters
!

subroutine phantomcells(nboxes,iatomfirst,nbdim) 

  implicit none
  integer :: i, j, k
  integer :: nboxes(3), nbdim(3), &
             iatomfirst(0:nbdim(1)+1,0:nbdim(2)+1,0:nbdim(3)+1)

  ! Vertices

  iatomfirst(0,0,0) = iatomfirst(nboxes(1),nboxes(2),nboxes(3))

  iatomfirst(nboxes(1)+1,0,0) = iatomfirst(1,nboxes(2),nboxes(3))
  iatomfirst(0,nboxes(2)+1,0) = iatomfirst(nboxes(1),1,nboxes(3))
  iatomfirst(0,0,nboxes(3)+1) = iatomfirst(nboxes(1),nboxes(2),1)

  iatomfirst(0,nboxes(2)+1,nboxes(3)+1) = iatomfirst(nboxes(1),1,1)
  iatomfirst(nboxes(1)+1,0,nboxes(3)+1) = iatomfirst(1,nboxes(2),1)
  iatomfirst(nboxes(1)+1,nboxes(2)+1,0) = iatomfirst(1,1,nboxes(3))

  iatomfirst(nboxes(1)+1,nboxes(2)+1,nboxes(3)+1) = iatomfirst(1,1,1)

  ! Axes 

  do i = 1, nboxes(1)
    iatomfirst(i,0,0) = iatomfirst(i,nboxes(2),nboxes(3))
    iatomfirst(i,nboxes(2)+1,0) = iatomfirst(i,1,nboxes(3))
    iatomfirst(i,0,nboxes(3)+1) = iatomfirst(i,nboxes(2),1)
    iatomfirst(i,nboxes(2)+1,nboxes(3)+1) = iatomfirst(i,1,1)
  end do
  do j = 1, nboxes(2)
    iatomfirst(0,j,0) = iatomfirst(nboxes(1),j,nboxes(3))
    iatomfirst(nboxes(1)+1,j,0) = iatomfirst(1,j,nboxes(3))
    iatomfirst(0,j,nboxes(3)+1) = iatomfirst(nboxes(1),j,1)
    iatomfirst(nboxes(1)+1,j,nboxes(3)+1) = iatomfirst(1,j,1)
  end do
  do k = 1, nboxes(3)
    iatomfirst(0,0,k) = iatomfirst(nboxes(1),nboxes(2),k)
    iatomfirst(nboxes(1)+1,0,k) = iatomfirst(1,nboxes(2),k)
    iatomfirst(0,nboxes(2)+1,k) = iatomfirst(nboxes(1),1,k)
    iatomfirst(nboxes(1)+1,nboxes(2)+1,k) = iatomfirst(1,1,k)
  end do

  ! Faces

  do j = 1, nboxes(2)
    do k = 1, nboxes(3)
      iatomfirst(0,j,k) = iatomfirst(nboxes(1),j,k)
      iatomfirst(nboxes(1)+1,j,k) = iatomfirst(1,j,k)
    end do
  end do
  do i = 1, nboxes(1)
    do k = 1, nboxes(3)
      iatomfirst(i,0,k) = iatomfirst(i,nboxes(2),k)
      iatomfirst(i,nboxes(2)+1,k) = iatomfirst(i,1,k)
    end do
  end do
  do i = 1, nboxes(1)
    do j = 1, nboxes(2)
      iatomfirst(i,j,0) = iatomfirst(i,j,nboxes(3))
      iatomfirst(i,j,nboxes(3)+1) = iatomfirst(i,j,1)
    end do
  end do

return
end subroutine phantomcells

!
! Subroutine transrot: Compute the coordinates of a molecule
! using rotations and translations relative to its center of mass 
!
!      beta is a counterclockwise rotation around x axis.
!      gamma is a counterclockwise rotation around y axis.
!      theta is a counterclockwise rotation around z axis.
!

subroutine compcart(natoms,xref,yref,zref,x,y,z,cmx,cmy,cmz,beta,gamma,theta)

  implicit none
  integer :: i, natoms
  real :: xref(*), yref(*), zref(*), x(*), y(*), z(*), &
          cmx, cmy, cmz, beta, gamma, theta, &
          v1(3), v2(3), v3(3),c1, s1, c2, s2, c3, s3
  
  ! Compute rotation matrix
  
  c1 = cos(beta) 
  s1 = sin(beta) 
  c2 = cos(gamma) 
  s2 = sin(gamma) 
  c3 = cos(theta) 
  s3 = sin(theta)
  
  v1(1) = c2*c3
  v1(2) = c1*s3 + c3*s1*s2
  v1(3) = s1*s3 - c1*c3*s2
  v2(1) = -c2*s3
  v2(2) = c1*c3 - s1*s2*s3
  v2(3) = c1*s2*s3 + c3*s1
  v3(1) = s2
  v3(2) = -c2*s1
  v3(3) = c1*c2         
  
  ! Apply rotation and translation
  
  do i = 1, natoms
    x(i) = cmx + xref(i)*v1(1) + yref(i)*v2(1) + zref(i)*v3(1)    
    y(i) = cmy + xref(i)*v1(2) + yref(i)*v2(2) + zref(i)*v3(2)    
    z(i) = cmz + xref(i)*v1(3) + yref(i)*v2(3) + zref(i)*v3(3)    
  end do
  
  return
end subroutine compcart

!
! Subroutine that returns the xmin and xmax vectors that contain the
! maximum and minimum coordinates of a group of atoms
!
! If init = .true. it will initialize the xmax and xmin vectors. Otherwise,
! the input values of xmin and xmax will be considered in the evaluation.
!

subroutine getmaxmin(ngroup,group,iatom,xdcd,ydcd,zdcd,xmin,xmax,init)

  integer :: ngroup, iatom, group(ngroup), ii, i, ifirst
  real :: xdcd(*), ydcd(*), zdcd(*), xmin(3), xmax(3)
  logical :: init

  if ( init ) then
    ii = iatom + group(1)
    xmin(1) = xdcd(ii)
    xmin(2) = ydcd(ii)
    xmin(3) = zdcd(ii)
    xmax(1) = xdcd(ii)
    xmax(2) = ydcd(ii)
    xmax(3) = zdcd(ii)
    ifirst = 2
  else
    ifirst = 1
  end if
  do i = ifirst, ngroup
    ii = iatom + group(i)
    xmin(1) = min(xmin(1),xdcd(ii))
    xmin(2) = min(xmin(2),ydcd(ii))
    xmin(3) = min(xmin(3),zdcd(ii))
    xmax(1) = max(xmax(1),xdcd(ii))
    xmax(2) = max(xmax(2),ydcd(ii))
    xmax(3) = max(xmax(3),zdcd(ii))
  end do

  return
end subroutine getmaxmin

!
! Determines the number of frames of a DCD file
!
! DCD comes open with header read before, and return
! at the same point
!
subroutine getnframes(iunit,nframes,dcdaxis,lastframe)

  integer :: iunit, nframes, status
  logical :: dcdaxis
  integer :: i, lastframe
  real :: x
  double precision :: side
  character(len=4) :: char

  write(*,*) ' Reading DCD file to determine number of frames ... '
  nframes = 0
  do
    if(dcdaxis) then 
      read(iunit,iostat=status) side
      if ( status /= 0 ) exit
    end if
    read(iunit,iostat=status) x
    if ( status /= 0 ) exit
    read(iunit,iostat=status) x
    if ( status /= 0 ) exit
    read(iunit,iostat=status) x
    if ( status /= 0 ) exit
    nframes = nframes + 1
    if ( lastframe /= 0 .and. nframes == lastframe ) exit
  end do
  rewind(iunit)
  read(iunit) char
  read(iunit) i
  read(iunit) i

  write(*,*) ' Total number of frames to read in this dcd file: ', nframes
  if(nframes < lastframe) then
    write(*,*) ' WARNING: lastframe greater than the number of frames of this dcd file. '
  end if

end subroutine getnframes

!
! Module with functions to operate on file names and strings
!

module file_operations

  contains

    !
    ! Function that determines the basename of a file,
    ! removing the path and the extension
    !
    
    character(len=200) function basename(filename)
     
      implicit none
      character(len=200) :: filename
    
      basename = remove_path(filename)
      basename = remove_extension(basename)
    
    end function basename
    
    !
    ! Function that removes the extension of a file name
    !
    
    character(len=200) function remove_extension(filename)
     
      implicit none
      integer :: i, idot
      character(len=200) :: filename
    
      remove_extension = filename
      i = length(remove_extension)
      idot = i+1
      do while(i > 0)
        if ( remove_extension(i:i) == "." ) then
          idot = i
        end if
        i = i - 1
      end do
      i = i + 1
      remove_extension = remove_extension(1:idot-1)
      do i = idot, 200
        remove_extension(i:i) = achar(32)
      end do
    
    end function remove_extension

    !
    ! Function that removes the path from a file name
    !
    
    character(len=200) function remove_path(filename)
    
      implicit none
      integer :: i, ilength
      character(len=200) :: filename
    
      remove_path = trim(adjustl(filename))
      ilength = length(remove_path)
      i = ilength
      do while(remove_path(i:i) /= "/")
        i = i - 1
        if ( i == 0 ) exit
      end do
      i = i + 1
      remove_path(1:ilength-i+1) = remove_path(i:ilength)
      do i = ilength-i+1, 200
        remove_path(i:i) = achar(32)
      end do
    
    end function remove_path
    
    !
    ! Function that return only the extension of the file
    !
    
    character(len=200) function file_extension(filename)
     
      implicit none
      integer :: i, idot
      character(len=200) :: filename
    
      i = length(filename)
      idot = 1
      do while(i > 0)
        if ( filename(i:i) == "." ) then
          idot = i
          exit
        end if
        i = i - 1
      end do
      file_extension(1:length(filename)-idot) = filename(idot+1:length(filename))
      do i = length(filename)-idot+1, 200
        file_extension(i:i) = achar(32)
      end do
    
    end function file_extension

    !
    ! Function that determines the length of a string
    !
    
    integer function length(string)
    
      implicit none
      character(len=200) :: string
      length = 200
      do while( empty_char(string(length:length)) ) 
        length = length - 1
      end do
    
    end function length
    
    !
    ! Function that determines if a character is empty
    !
    
    logical function empty_char(char)
    
      implicit none
      character :: char
      empty_char = .false.
      if ( char == achar(9) .or. &
           char == achar(32) .or. &
           char == '' ) then
        empty_char = .true.
      end if 
    
    end function empty_char
    
    !
    ! Function that checks if a line is a comment line
    !
    
    logical function comment(string)
      
      implicit none
      integer :: i
      character(len=200) :: string
      i = 1
      do while( empty_char(string(i:i)) .and. i < 200 ) 
        i = i + 1
      end do
      comment = .false.
      if ( string(i:i) == "#" .or. i == 200 ) comment = .true.
    
    end function comment

    !
    ! Subroutine that tests if file exists, and tries to open it.
    ! Asks the user if he/she wants the file to be overwritten 
    !

    subroutine checkfile(file)
    
      implicit none
      integer :: ioerr
      character(len=200) :: file
      character(len=1) :: char
    
      open(10,file=file,status='new',action='write',iostat=ioerr)
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Trying to create file: ', trim(adjustl(file)),' but file already exists '
        write(*,"(a,$)") '  Overwrite it? (Y/N): '
        read(*,*) char
        if ( char == "Y" ) then
          open(10,file=file,action='write',iostat=ioerr)
          if ( ioerr /= 0 ) then
            write(*,*) ' Could not open file. Quitting. '
            stop
          end if
          close(10)
        else
          write(*,*) ' Quitting. '
          stop
        end if
      end if
      close(10)
    
    end subroutine checkfile

end module file_operations

!
! Function that returns a real random number
!

real function random()
  implicit none
  call random_number(random)
  return
end function random

!
! Subroutine that returns a random number given a the seed
!

subroutine init_random_number(iseed)
  implicit none
  integer :: i, iseed, size
  integer, allocatable :: seed(:)
  call random_seed(size=size)
  allocate(seed(size))
  do i = 1, size
    seed(i) = i*iseed
  end do
  call random_seed(put=seed)
  deallocate(seed)
  return
end subroutine init_random_number

!
! Subroutine that uses the date to create a random seed
! 

subroutine seed_from_time(seed)

  implicit none
  integer :: seed, value(8)
  character(len=10) :: b(3)
  call date_and_time( b(1), b(2), b(3), value )
  seed = value(1)+value(2)+value(3)+value(4)+value(5)+value(6)+value(7)+value(8)
  seed = seed + value(1)+value(2)+value(3)+value(4)+value(5)/100+value(6)*100+value(7)/10+value(8)*10

end subroutine seed_from_time










