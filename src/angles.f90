!
! angles: A program for analysing angles and dihedrals from groups
!         in molecular dynamics simulations in NAMD DCD format.
!
! Auxiliar dimensions:
!          memory: Controls the amount of memory used for reading dcd
!                  files. If the program ends with segmentation fault
!                  without any other printing, decrease the value of
!                  this parameter.  
!
! L. Martinez, Mar 17, 2008.
!
! Dynamically allocatable arrays and their dimensions:
!
!       group1: number of atoms of group1
!       eps, sig, q, e, s, class, segat, resat, typeat, classat: total
!               number of atoms of the system
!

! Static variables
 
use charsize
implicit none
integer, parameter :: memory=15000000
real, parameter :: pi=3.141592654
integer :: natom, ngroup1,&
           narg, length, firstframe, lastframe, stride, nclass,&
           nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
           j, iframe, icycle, nfrcycle, iatom, k, ii, catom,&
           status, keystatus, ndim, iargc, lastatom, group1(4)
double precision :: readsidesx, readsidesy, readsidesz, t
real :: side(memory,3), mass1
real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
        x1, y1, z1, time0, etime, tarray(2), angle2,& 
        xang(4), yang(4), zang(4), angle, dihedral, xnorm, cosang,&
        scaletime, axis(3)
character(len=200) :: groupfile, line, record, value, keyword,&
                      dcdfile, inputfile, output, psffile
character(len=4) :: dummyc
logical :: periodic, readfromdcd, dcdaxis, centeratom      

! Allocatable arrays

integer, allocatable :: resid(:)
real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:)
character(len=charsize1), allocatable :: class(:), segat(:), resat(:),&
                                         typeat(:), classat(:)

! Compute time

time0 = etime(tarray)

! Output title

write(*,"(/,' ####################################################',&
          &/,/,&
          & '   ANGLES: Compute angles or dihedrals from DCD files.',&
          &/,/,&
          & ' ####################################################',/)")    

call version()

! Some default parameters

firstframe = 1
lastframe = 0
stride = 1
periodic = .true.
readfromdcd = .true.
centeratom = .false.
catom = 0   
scaletime = 1.

! Open input file and read parameters

narg = iargc()
if(narg == 0) then
  write(*,*) ' Run with: ./angles input.inp '
  stop
end if   
call getarg(1,record)

inputfile = record(1:length(record))
open(99,file=inputfile,action='read')
do 
  read(99,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(keyword(record) == 'dcd') then
    dcdfile = value(record)
    write(*,*) ' DCD file name: ', dcdfile(1:length(dcdfile)) 
  else if(keyword(record) == 'groups') then
    groupfile = value(record)
    write(*,*) ' Groups file name: ', groupfile(1:length(groupfile))
  else if(keyword(record) == 'psf') then
    psffile = value(record)
    write(*,*) ' PSF file name: ', psffile(1:length(psffile))
  else if(keyword(record) == 'firstframe') then
    line = value(record)
    read(line,*,iostat=keystatus) firstframe
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'scaletime') then
    line = value(record)
    read(line,*,iostat=keystatus) scaletime
    if(keystatus /= 0) exit
  else if(keyword(record) == 'lastframe') then
    line = value(record)
    if(line(1:length(line)) /= 'last') then
      read(line,*,iostat=keystatus) lastframe
    end if
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'stride') then
    line = value(record)
    read(line,*,iostat=keystatus) stride
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'centeratom') then
    line = value(record)
    read(line,*,iostat=keystatus) catom
    if(keystatus /= 0) exit 
    centeratom = .true.
  else if(keyword(record) == 'periodic') then
    line = value(record)
    read(line,*,iostat=keystatus) line
    if(keystatus /= 0) exit 
    if(line == 'no') then
      periodic = .false.
      readfromdcd = .false.
    else if(line == 'readfromdcd') then
      periodic = .true.
      readfromdcd = .true.
    else
      periodic = .true.
      readfromdcd = .false.
      read(record,*,iostat=keystatus) line, axis(1), axis(2), axis(3)
      if(keystatus /= 0) exit 
    end if 
  else if(keyword(record) == 'output') then
    output = value(record)
    write(*,*) ' Output file name: ', output(1:length(output))
  else if(keyword(record) == 'atoms') then
    write(*,"(a,/,a)") ' ERROR: The option atoms must be used ',&
                       '        with the angles.sh script, not directly. '
    stop
  else if(record(1:1) /= '#'.and.record(1:1) > ' ') then
    write(*,*) ' ERROR: Unrecognized keyword found: ',record(1:length(record))
    stop
  end if
end do               
close(99)

! If some error was found in some keyword value, report error and stop

if(keystatus /= 0) then
  line = keyword(record)
  write(*,*) ' ERROR: Could not read value for keyword: ',line(1:length(line))
  stop
end if

! Reading the header of psf file and parameter files to get dimensions

call getdim(psffile,inputfile,ndim)
allocate( eps(ndim), sig(ndim), q(ndim), e(ndim), s(ndim), mass(ndim),&
          segat(ndim), resat(ndim), classat(ndim), typeat(ndim),&
          class(ndim), resid(ndim) )
       
! Check for simple input errors

if(stride < 1) then
  write(*,*) ' ERROR: stride cannot be less than 1. ' 
  stop
end if
if(lastframe < firstframe.and.lastframe /= 0) then
  write(*,*) ' ERROR: lastframe must be greater or equal to firstframe. '
  stop
end if
if(periodic .and. .not. centeratom) then
  write(*,*) ' ERROR: Periodic options must be accompanied with the '
  write(*,*) '        centeratom option to define the center of the '
  write(*,*) '        periodic cell in this module. '
  stop
end if

! Output some information if not set

write(*,*) ' First frame to be considered: ', firstframe
if(lastframe == 0) then
  write(*,*) ' Last frame to be considered: last '
else
  write(*,*) ' Last frame to be considered: ', lastframe
end if
write(*,*) ' Stride (will jump frames): ', stride

! Read PSF file

write(*,*) ' Reading PSF file: ', psffile(1:length(psffile))
call readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
             resid,classat,typeat,q,e,s,mass,.false.)
write(*,*) ' Number of atoms in PSF file: ', natom
write(*,*)

! Read group information from group file
! Counting the number of atoms with non-zero beta to check for errors

open(10,file=groupfile,action='read')
ngroup1 = 0
do 
  read(10,"( a200 )",iostat=status) line
  if(status /= 0) exit
  if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
    if(line(63:66) /= '0.00') ngroup1 = ngroup1 + 1        
    if(line(57:60) /= '0.00') ngroup1 = ngroup1 + 1        
  end if     
end do
close(10)
if(ngroup1 /= 3 .and. ngroup1 /= 4) then
  write(*,*) ' ERROR: Number of atoms in selection must be 3 or 4.'
  stop
end if

! Now reading reading the group atoms

open(10,file=groupfile,action='read')
mass1 = 0.d0
iatom = 0
do i = 1, 4
  group1(i) = 0
end do
do 
  read(10,"( a200 )",iostat=status) line
  if(status /= 0) exit
  if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
    iatom = iatom + 1
    if(line(63:66) /= '0.00'.or. &
       line(57:66) /= '0.00') then
      if(line(57:60) == '1.00') group1(1) = iatom
      if(line(57:60) == '2.00') group1(2) = iatom
      if(line(63:66) == '3.00') group1(3) = iatom
      if(line(63:66) == '4.00') group1(4) = iatom
      mass1 = mass1 + mass(iatom)
    end if
  end if     
end do
close(10)
lastatom = catom
do i = 1, ngroup1 
  lastatom = max0(lastatom,group1(i))
end do
do i = 1, ngroup1
  if(group1(i) == 0) then
    write(*,*) ' ERROR: Some atom was not correctly defined. '
    stop
  end if
end do

! Output some group properties for testing purposes

write(*,*) ' Number of atoms of group 1: ', ngroup1 
write(*,*) ' First atom of group 1: ', group1(1)
write(*,*) ' Last atom of group 1: ', group1(ngroup1)
write(*,*) ' Mass of group 1: ', mass1
if(ngroup1 < 3) then
  write(*,*) ' ERROR: Selection must contain at least three atoms.'
  stop
end if
write(*,*) ' Time (x-axis) will be scaled by: ', scaletime

 
! Checking if dcd file contains periodic cell information

write(*,"( /,tr2,52('-') )")
write(*,*) ' Periodic cell data: Read carefully. '
call chkperiod(dcdfile,dcdaxis,readfromdcd) 
if(.not.readfromdcd.and.periodic) then
  write(*,*) ' User provided periodic cell dimensions: '
  write(*,*) axis(1), axis(2), axis(3)
end if
if(centeratom) then
  write(*,*) ' All interactions will be computed by wrapping' 
  write(*,*) ' coordinates relative to atom ', catom,' ',typeat(catom),&
             ' ',resat(catom),' ',segat(catom)
  write(*,*) ' This means that this is NOT a usual minimum image computation.'
end if
if(.not.periodic.and.centeratom) then
  write(*,*) ' ERROR: centeratom option can only be used '
  write(*,*) ' if some periodic condition is to be applied. '
  stop
end if   

! Reading the dcd file

write(*,"(tr2,52('-'),/ )")
write(*,*) ' Reading the DCD file header: '
open(10,file=dcdfile,action='read',form='unformatted')
read(10) dummyc, nframes, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
read(10) dummyi, dummyr
read(10) ntotat

write(*,*)
write(*,*) ' Number of atoms as specified in the dcd file: ',ntotat     
write(*,*) ' Total number of frames in this dcd file: ', nframes
if(nframes < lastframe) then
  write(*,*) ' ERROR: lastrame greater than the number of '
  write(*,*) '        frames of the dcd file. '
  stop
end if
if(lastframe == 0) lastframe = nframes    
if(ntotat /= natom) then
  write(*,"(a,/,a)") ' ERROR: Number of atoms in the dcd file does not',&
                    &'        match the number of atoms in the psf file'
  stop
end if
 
! Open output file and writes all information of this run

open(51,file=output(1:length(output)))
write(51,"( '#',/,&
           &'# Output of angles.f90:',/,&
           &'# Input file: ',a,/,& 
           &'# DCD file: ',a,/,& 
           &'# Group file: ',a,/,&
           &'# PSF file: ',a,/,& 
           &'# First frame: ',i5,' Last frame: ',i5,' Stride: ',i5,/,&
           &'#',/,&
           &'# Periodic conditions: ',/,&
           &'# Periodic: ',l1,' Read from DCD: ',l1,&
           &   ' Center atom: ',l1,' ',i6,/,&
           &'#',/,&
           &'# Number of atoms and mass of group 1: ',i6,f12.3,/,&
           &'# First and last atoms of group 1: ',i6,tr1,i6,/,&
           &'#',/ )")&
           &inputfile(1:length(inputfile)),&
           &dcdfile(1:length(dcdfile)),&
           &groupfile(1:length(groupfile)),&
           &psffile(1:length(psffile)),&
           &firstframe, lastframe, stride,&
           &periodic, readfromdcd, centeratom, catom,&
           &ngroup1, mass1, group1(1), group1(ngroup1) 

if(ngroup1.eq.3) write(51,"('#  TIME ',t23,'ANGLE')") 
if(ngroup1.eq.4) write(51,"('#  TIME ',t20,'DIHEDRAL',t36,'COLINEARITY')") 

! Now going to read the dcd file

memframes = memory / ntotat
ncycles = lastframe / memframes + 1
memlast = lastframe - memframes * ( ncycles - 1 )
write(*,*) ' Will read and store in memory at most ', memframes,&
           ' frames per reading cycle. '
write(*,*) ' There will be ', ncycles, ' cycles of reading. '
write(*,*) ' Last cycle will read ', memlast,' frames. '
write(*,*)        

! Beautiful dcd reading output

write(*,"( t3,'Cycle',t17,'|0%' )",advance='no') 
j = 0
do i = 1, min0(memframes,lastframe)
  if(mod(i,max0(1,min0(memframes,lastframe)/35)) == 0) j = j + 1
end do
do i = 1, j - 6
  write(*,"( ' ' )",advance='no')
end do
write(*,"( '100%|' )")

! Reading dcd file and computing angles
 
iframe = 0
do icycle = 1, ncycles 

  write(*,"( t3,i4,' Read...',t17,'|' )",advance='no') icycle

! Each cycle fills the memory as specified by the memory parameter 

  if(icycle == ncycles) then
    nfrcycle = memlast
  else
    nfrcycle = memframes
  end if

  iatom = 0
  do j = 1, nfrcycle    
    if(dcdaxis) then
      read(10) readsidesx, t, readsidesy, t, t, readsidesz
      side(j,1) = sngl(readsidesx)
      side(j,2) = sngl(readsidesy)
      side(j,3) = sngl(readsidesz)
    end if
    read(10) (xdcd(k), k = iatom + 1, iatom + lastatom)
    read(10) (ydcd(k), k = iatom + 1, iatom + lastatom)            
    read(10) (zdcd(k), k = iatom + 1, iatom + lastatom)           
    iatom = iatom + ntotat
    if(mod(j,max0(1,nfrcycle/35)) == 0) write(*,"( '*' )",advance='no') 
  end do
  write(*,"( '|' )")   
  write(*,"( t3,'Computing... ',t17,'|' )",advance='no') 

! Computing the angles

  iatom = 0
  do k = 1, nfrcycle
    iframe = iframe + 1

! Reseting the variables that contains the angles in each frame

    if(mod(iframe,stride) == 0 .and. iframe >= firstframe) then

      do i = 1, ngroup1
        ii = iatom + group1(i)

! Move atom ii if centeratom option is set 

        if(centeratom) then
          x1 = xdcd(ii) - xdcd(iatom+catom)
          y1 = ydcd(ii) - ydcd(iatom+catom)
          z1 = zdcd(ii) - zdcd(iatom+catom)
          if(readfromdcd) then
            call image(x1,y1,z1,side(k,1),side(k,2),side(k,3))
          else
            call image(x1,y1,z1,axis(1),axis(2),axis(3))
          end if
        else
          x1 = xdcd(ii)
          y1 = ydcd(ii)
          z1 = zdcd(ii)
        end if
        xang(i) = x1
        yang(i) = y1
        zang(i) = z1

      end do

! Computing the angle

      if(ngroup1.eq.3) then

        angle = cosang(xang(1)-xang(2),&
                       yang(1)-yang(2),&
                       zang(1)-zang(2),&
                       xang(3)-xang(2),&
                       yang(3)-yang(2),&
                       zang(3)-zang(2))

        angle = 180*acos(angle)/pi
        if(angle > 180) angle = angle - 180
        write(51,"( i8,tr2,f17.6 )") scaletime*iframe, angle

      else if(ngroup1.eq.4) then

        angle = dihedral(xang(1),yang(1),zang(1),&
                         xang(2),yang(2),zang(2),&
                         xang(3),yang(3),zang(3),&
                         xang(4),yang(4),zang(4))

        angle2 = ( (xang(4)-xang(3)) * (xang(2)-xang(1)) ) + &
                 ( (yang(4)-yang(3)) * (yang(2)-yang(1)) ) + &
                 ( (zang(4)-zang(3)) * (zang(2)-zang(1)) )
        angle2 = angle2 / (xnorm(xang(2)-xang(1),   &
                                 yang(2)-yang(1),   &
                                 zang(2)-zang(1)) * &
                           xnorm(xang(4)-xang(3),   &
                                 yang(4)-yang(3),   &
                                 zang(4)-zang(3)))  
        write(51,"( e12.6,2(tr2,f17.6) )") scaletime*iframe, angle, angle2
      
      end if


    end if

! Printing status bar

    if(mod(k,max0(1,nfrcycle/35)) == 0) write(*,"( '*' )",advance='no') 

    iatom = iatom + ntotat
  end do
  write(*,"( '|' )")
         
end do
close(10)
close(51)

! Write final messages with names of output files and their content

time0 = etime(tarray) - time0
write(*,*)
write(*,"( tr2,52('-') )")
write(*,*)
write(*,*) ' OUTPUT FILES: ' 
write(*,*)
write(*,*) ' Wrote output file: ', output(1:length(output))
write(*,*)
write(*,*) ' Which contains the desired angle as a function'
write(*,*) ' of the simulation frame. '
write(*,*) 
write(*,*) ' Running time: ', time0
write(*,*) '####################################################'
write(*,*) 
write(*,*) '  END: Normal termination.  '
write(*,*) 
write(*,*) '####################################################'
write(*,*)        

end

!
! Compute the dihedral angle defined by four atoms
!

function dihedral(x0,y0,z0,&
                  x1,y1,z1,&
                  x2,y2,z2,&
                  x3,y3,z3)

implicit none
real, parameter :: pi=3.141592654
real :: dihedral, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3,&
        n1(3), n2(3), sign, s(3)

call vprod(n1,x1-x0,y1-y0,z1-z0,x2-x1,y2-y1,z2-z1)
call vprod(n2,x2-x1,y2-y1,z2-z1,x3-x2,y3-y2,z3-z2)
call vprod(s,n1(1),n1(2),n1(3),n2(1),n2(2),n2(3))
sign = s(1)*(x2-x1) + s(2)*(y2-y1) + s(3)*(z2-z1)

dihedral = n1(1)*n2(1) + n1(2)*n2(2) + n1(3)*n2(3)
dihedral = 180.*acos(amin1(1.,dihedral))/pi

if(sign < 0) dihedral = -1. * abs(dihedral)
if(sign > 0) dihedral = abs(dihedral)

return
end

!
! Compute a product vector and normalize
!

subroutine vprod(v,x1,x2,x3,y1,y2,y3)

implicit none
real :: v(3), vnorm, x1, x2, x3, y1, y2, y3, xnorm

v(1) = x2*y3 - x3*y2
v(2) = - x1*y3 + x3*y1
v(3) = x1*y2 - x2*y1
vnorm = xnorm(v(1),v(2),v(3))
v(1) = v(1) / vnorm
v(2) = v(2) / vnorm
v(3) = v(3) / vnorm

return
end
                  





 
