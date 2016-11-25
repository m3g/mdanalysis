!
! distance: A program for analysing distances between selected groups
!           in molecular dynamics simulations in NAMD DCD format.
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
!       group2: number of atoms of group2
!       eps, sig, q, e, s, class, segat, resat, typeat, classat: total
!               number of atoms of the system
!

! Static variables
 
use charsize
implicit none
integer, parameter :: memory=15000000
integer :: natom, ngroup1, ngroup2,&
           narg, length, firstframe, lastframe, stride, nclass,&
           nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
           j, iframe, icycle, nfrcycle, iatom, k, ii, jj, catom,&
           status, keystatus, ndim, iargc, lastatom
double precision :: readsidesx, readsidesy, readsidesz, t
real :: side(memory,3), sij,&
        mindist, maxdist, cmdist, mass1, mass2,&
        cm1x, cm1y, cm1z, cm2x, cm2y, cm2z
real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
        x1, y1, z1, x2, y2, z2, time0, etime, tarray(2),&
        scaletime, axis(3)
character(len=200) :: groupfile, line, record, value, keyword,&
                      dcdfile, inputfile, output, psffile
character(len=4) :: dummyc
logical :: periodic, readfromdcd, dcdaxis, centeratom      

! Allocatable arrays

integer, allocatable :: group1(:), group2(:), resid(:)
real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:)
character(len=charsize1), allocatable :: class(:), segat(:), resat(:),&
                                         typeat(:), classat(:)

! Compute time

time0 = etime(tarray)

! Output title

write(*,"(/,' ####################################################',&
          &/,/,&
          & '   DISTANCE: Compute distances from DCD files  ',&
          &/,/,&
          & ' ####################################################',/)")    

call version

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
  write(*,*) ' Run with: ./distance input.inp '
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
  else if(keyword(record) == 'lastframe') then
    line = value(record)
    if(line(1:length(line)) /= 'last') then
      read(line,*,iostat=keystatus) lastframe
    end if
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'scaletime') then
    line = value(record)
    read(line,*,iostat=keystatus) scaletime
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
  else if(keyword(record) == 'group1' .or. &
          keyword(record) == 'group2') then
    write(*,"(a,/,a)") ' ERROR: The options group1 and group2 must be used ',&
                       '        with the distance.sh script, not directly. '
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

! Output some information if not set

write(*,*) ' First frame to be considered: ', firstframe
if(lastframe == 0) then
  write(*,*) ' Last frame to be considered: last '
else
  write(*,*) ' Last frame to be considered: ', lastframe
end if
write(*,*) ' Stride (will jump frames): ', stride
write(*,*) ' Time (x-axis) will be scaled by: ', scaletime

! Read PSF file

write(*,*) ' Reading PSF file: ', psffile(1:length(psffile))
call readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
             resid,classat,typeat,q,e,s,mass,.false.)
write(*,*) ' Number of atoms in PSF file: ', natom
write(*,*)

! Read group information from group file
! First reading the size of the groups to allocate arrays

open(10,file=groupfile,action='read')
ngroup1 = 0
ngroup2 = 0
do 
  read(10,"( a200 )",iostat=status) line
  if(status /= 0) exit
  if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
    if(line(63:66) == '1.00') then      
      ngroup1 = ngroup1 + 1        
    else if(line(63:66) == '2.00') then
      ngroup2 = ngroup2 + 1        
    end if
  end if     
end do
allocate ( group1(ngroup1), group2(ngroup2) )
close(10)

! Now reading reading the group atoms

open(10,file=groupfile,action='read')
ngroup1 = 0
ngroup2 = 0
mass1 = 0.d0
mass2 = 0.d0
iatom = 0
do 
  read(10,"( a200 )",iostat=status) line
  if(status /= 0) exit
  if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
    iatom = iatom + 1
    if(line(63:66) == '1.00') then      

! Read atoms belonging to group1

      ngroup1 = ngroup1 + 1        
      group1(ngroup1) = iatom
      mass1 = mass1 + mass(iatom)

    else if(line(63:66) == '2.00') then

! Read atoms belonging to group2

      ngroup2 = ngroup2 + 1        
      group2(ngroup2) = iatom
      mass2 = mass2 + mass(iatom)

    end if
  end if     
end do
close(10)
lastatom = max0(group1(ngroup1),group2(ngroup2))
lastatom = max0(lastatom,catom)

! Output some group properties for testing purposes

write(*,*) ' Number of atoms of group 1: ', ngroup1 
write(*,*) ' First atom of group 1: ', group1(1)
write(*,*) ' Last atom of group 1: ', group1(ngroup1)
write(*,*) ' Mass of group 1: ', mass1
write(*,*) ' Number of atoms of group 2: ', ngroup2 
write(*,*) ' First atom of group 2: ', group2(1)
write(*,*) ' Last atom of group 2: ', group2(ngroup2)
write(*,*) ' Mass of group 2: ', mass2
 
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

! Stopping if some group is empty

if(ngroup1 == 0 .or. ngroup2 == 0) then
  write(*,*) ' ERROR: Groups must contain at least one atom.'
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
call getnframes(10,nframes,dcdaxis)
write(*,*) ' Total number of frames in this dcd file: ', nframes
if(nframes < lastframe) then
  write(*,*) ' ERROR: lastframe greater than the number of '
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
           &'# Output of distance.f:',/,&
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
           &'# Number of atoms and mass of group 2: ',i6,f12.3,/,&
           &'# First and last atoms of group 2: ',i6,tr1,i6,/,&
           &'#',/,&
           &'# Distances of group 1 with group 2',/,&
           &'#  TIME ',t12,'MINIMUM DISTANCE',&
           &           t36,'CM DISTANCE',&
           &           t50,'MAXIMUM DISTANCE')")& 
           &inputfile(1:length(inputfile)),&
           &dcdfile(1:length(dcdfile)),&
           &groupfile(1:length(groupfile)),&
           &psffile(1:length(psffile)),&
           &firstframe, lastframe, stride,&
           &periodic, readfromdcd, centeratom, catom,&
           &ngroup1, mass1, group1(1), group1(ngroup1),& 
           &ngroup2, mass2, group2(1), group2(ngroup2) 

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

! Reading dcd file and computing distances
 
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

! Computing the distances

  iatom = 0
  do k = 1, nfrcycle
    iframe = iframe + 1

! Reseting the variables that contains the distances in each frame

    if(mod(iframe,stride) == 0 .and. iframe >= firstframe) then
      mindist = 1.e30
      maxdist = 0.d0
      cm1x = 0.d0
      cm1y = 0.d0
      cm1z = 0.d0
      cm2x = 0.d0
      cm2y = 0.d0
      cm2z = 0.d0
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

! Computing center of mass of group 1

        cm1x = cm1x + x1*mass(group1(i))
        cm1y = cm1y + y1*mass(group1(i))
        cm1z = cm1z + z1*mass(group1(i))

        do j = 1, ngroup2
          jj = iatom + group2(j)

! Move atom jj according to reference (either atom ii or centeratom)

          if(periodic) then
            if(centeratom) then
              x2 = xdcd(jj) - xdcd(iatom+catom)
              y2 = ydcd(jj) - ydcd(iatom+catom)
              z2 = zdcd(jj) - zdcd(iatom+catom)
              if(readfromdcd) then
                call image(x2,y2,z2,side(k,1),side(k,2),side(k,3))
              else
                call image(x2,y2,z2,axis(1),axis(2),axis(3))
              end if
            else
              x2 = xdcd(jj) - x1
              y2 = ydcd(jj) - y1
              z2 = zdcd(jj) - z1
              if(readfromdcd) then
                call image(x2,y2,z2,side(k,1),side(k,2),side(k,3))
              else
                call image(x2,y2,z2,axis(1),axis(2),axis(3))
              end if
              x2 = x2 + x1
              y2 = y2 + y1
              z2 = z2 + z1
            end if
          else
            x2 = xdcd(jj)
            y2 = ydcd(jj)
            z2 = zdcd(jj)
          end if

! Computing the distance between ii and jj (atom jj is in the origin)

          sij = sqrt(( x2 - x1 )**2 + ( y2 - y1 )**2 + ( z2 - z1 )**2)
          if(sij < mindist) mindist = sij
          if(sij > maxdist) maxdist = sij

! Computing the center of mass of group 2 (add only if i = 1)
          
          if(i.eq.1) then
            cm2x = cm2x + x2*mass(group2(j))
            cm2y = cm2y + y2*mass(group2(j))
            cm2z = cm2z + z2*mass(group2(j))
          end if

        end do
      end do

! The distance between the centers of mass

      cm1x = cm1x / mass1      
      cm1y = cm1y / mass1      
      cm1z = cm1z / mass1      
      cm2x = cm2x / mass2      
      cm2y = cm2y / mass2      
      cm2z = cm2z / mass2      
      cmdist = sqrt(( cm2x - cm1x )**2 + ( cm2y - cm1y )**2 + ( cm2z - cm1z )**2)
 
! Printing to output file in this frame 
 
      write(51,"( e12.6,3(tr2,f17.6) )") scaletime*iframe, mindist, cmdist, maxdist

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
write(*,*) ' Which contains the minimum distance and the distance'
write(*,*) ' between the centers of mass of each group as a '
write(*,*) ' function of simulation frames. '
write(*,*) 
write(*,*) ' Running time: ', time0
write(*,*) '####################################################'
write(*,*) 
write(*,*) '  END: Normal termination.  '
write(*,*) 
write(*,*) '####################################################'
write(*,*)        

end

