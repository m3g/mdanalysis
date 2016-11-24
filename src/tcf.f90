!
! time_correlation: Compute the reorientational dynamics of a vector, defined
!                   as the average vector of a set of vectors defined as
!                   differences in positions between atoms in the simulation.
!                   (if only one pair of atoms is set, then the result will
!                    be the reorientational dynamics of this vector).
!
! L. Martinez, Jul 31, 2009.
!

program time_correlation

! Static variables
 
implicit none
integer :: narg, length, firstframe, lastframe, &
           dummyi, i, ntotat, nvectors, &
           j, catom, ntotframes, idcd, ndcd, iframe, ifcount, ifcount2,&
           status, keystatus, iargc, lastatom, navs, nave, nevs, neve
double precision :: t, readsidesx, readsidesy, readsidesz
real :: axis(3)
real :: dummyr,& 
        xt, yt, zt, r(3,3), ax, ay, az,&
        xt_avs, yt_avs, zt_avs,&
        xt_ave, yt_ave, zt_ave,&
        xt_evs, yt_evs, zt_evs,&
        xt_eve, yt_eve, zt_eve,&
        vnorm, int_prod, theta, cost, sint, r0,&
        scaletime
character(len=200) :: line, record, value, keyword, inputfile, output, record2
character(len=4) :: dummyc
logical :: periodic, readfromdcd, dcdaxis, centeratom, error

! Allocatable arrays

integer, allocatable :: iavs(:), iave(:), ievs(:), ieve(:), nframes(:)
real, allocatable :: xabs(:), yabs(:), zabs(:),&
                     xemi(:), yemi(:), zemi(:),&
                     tcf(:), legendre(:),&
                     xdcd(:), ydcd(:), zdcd(:)
real, allocatable :: side(:,:)
character(len=200), allocatable :: dcdfile(:)

! Output title

write(*,"(/,' ####################################################',&
          &/,/,&
          & '   ORIENTATION: Time-dependent reorientation ',&
          &/,/,&
          & ' ####################################################',/)")    

call version()

! Some default parameters

nvectors = 0
firstframe = 1
lastframe = 0
periodic = .false.
readfromdcd = .true.
centeratom = .false.
catom = 0   
theta = 0.d0
r0 = 1.0
scaletime = 1.

! Open input file and read parameters

narg = iargc()
if(narg == 0) then
  write(*,*) ' Run with: ./tcf input.inp '
  stop
end if   
call getarg(1,record)

inputfile = record(1:length(record))
open(99,file=inputfile,action='read')

! Count the number of dcd files

ndcd = 0
do
  read(99,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(keyword(record) == 'dcd') ndcd = ndcd + 1
end do
allocate(dcdfile(ndcd),nframes(ndcd)) ; ndcd = 0
rewind(99)

! Read all parameters from the input file

do 
  read(99,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(keyword(record) == 'dcd') then
    ndcd = ndcd + 1
    dcdfile(ndcd) = value(record)
    record2 = dcdfile(ndcd)
    write(*,*) ' DCD file name: ', record2(1:length(dcdfile(ndcd))) 
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
  else if(keyword(record) == 'theta') then
    line = value(record)
    read(line,*,iostat=keystatus) theta
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'r0') then
    line = value(record)
    read(line,*,iostat=keystatus) r0
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'scaletime') then
    line = value(record)
    read(line,*,iostat=keystatus) scaletime
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
  else if(keyword(record) == 'abs_vector_start') then
    read(record,*,iostat=keystatus) line, navs
    if(keystatus /= 0) exit 
    allocate(iavs(navs))
    read(record,*,iostat=keystatus) line, navs, (iavs(i), i = 1, navs)
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'abs_vector_end') then
    read(record,*,iostat=keystatus) line, nave
    if(keystatus /= 0) exit 
    allocate(iave(nave))
    read(record,*,iostat=keystatus) line, nave, (iave(i), i = 1, nave)
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'emi_vector_start') then
    read(record,*,iostat=keystatus) line, nevs
    if(keystatus /= 0) exit 
    allocate(ievs(nevs))
    read(record,*,iostat=keystatus) line, nevs, (ievs(i), i = 1, nevs)
    if(keystatus /= 0) exit 
  else if(keyword(record) == 'emi_vector_end') then
    read(record,*,iostat=keystatus) line, neve
    if(keystatus /= 0) exit 
    allocate(ieve(neve))
    read(record,*,iostat=keystatus) line, neve, (ieve(i), i = 1, neve)
    if(keystatus /= 0) exit 
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

! Find last atom to speed up reading

lastatom = 0
do i = 1, navs
  lastatom = max(lastatom,iavs(i))
end do
do i = 1, nave
  lastatom = max(lastatom,iave(i))
end do
do i = 1, nevs
  lastatom = max(lastatom,ievs(i))
end do
do i = 1, neve
  lastatom = max(lastatom,ieve(i))
end do

write(*,*) ' Absorption vector starts at average position of atoms: '
write(*,*) (iavs(i),i=1,navs)
write(*,*) ' Absorption vector ends at average position of atoms: '
write(*,*) (iave(i),i=1,nave)
write(*,*) ' Emission vector starts at average position of atoms: '
write(*,*) (ievs(i),i=1,nevs)
write(*,*) ' Emission vector ends at average position of atoms: '
write(*,*) (ieve(i),i=1,neve)

! Check for simple input errors

if(lastframe < firstframe .and. lastframe /= 0) then
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
write(*,*) ' Angle between absorption and emission dipoles (theta): ', theta
write(*,*) ' Intrisic (r0) anisotropy: ', r0
write(*,*) ' Time (x-axis) will be scaled by: ', scaletime

! Computing sint and cost for rotation angle correction

if ( theta > 1.d-7 ) then

  ! This only makes sense if the starting point of the two vectors are the same

  error = .false.
  if ( navs == nevs ) then 
    do i = 1, navs
      if ( iavs(i) /= ievs(i) ) then 
        error = .true.
      end if
    end do
  else
    error = .true.
  end if
  if (error) then
    write(*,*) ' Error: If you want to use an angle (theta) to define the emission '
    write(*,*) '        vector, the initial point of the emission vector ' 
    write(*,*) '        must be the same as the initial point of the absorption vector. '
    stop
  end if

  ! And also, this only makes sense if the vectors are not identical

  error = .true.
  if ( navs == nevs .and. nave == neve ) then
    do i = 1, navs
      if ( iavs(i) /= ievs(i) ) error = .false.
    end do
    do i = 1, nave 
      if ( iave(i) /= ieve(i) ) error = .false.
    end do
  else
    error = .false.
  end if
  if (error) then
    write(*,*) ' Error: If you want to compute the emission vector from a rotation theta '
    write(*,*) '        of the absorption vector, the defined absorption and emission '
    write(*,*) '        vectors cannot be equal, because they must define the plane '
    write(*,*) '        of rotation (actually, this is the only function of the definition '
    write(*,*) '        of the emission vector in this case).'
    stop
  end if

  ! Converting angle

  theta = 3.141592654 * theta / 180. ! theta in radians
  cost = cos(theta)
  sint = sin(theta)

end if

! Checking if dcd file contains periodic cell information

record2 = dcdfile(1)
write(*,"( /,tr2,52('-') )")
write(*,*) ' Periodic cell data: Read carefully. '
call chkperiod(record2,dcdaxis,readfromdcd) 
if(.not.readfromdcd.and.periodic) then
  write(*,*) ' User provided periodic cell dimensions: '
  write(*,*) axis(1), axis(2), axis(3)
end if
if(centeratom) then
  write(*,*) ' All interactions will be computed by wrapping' 
  write(*,*) ' coordinates relative to atom ', catom
  write(*,*) ' This means that this is NOT a usual minimum image computation.'
end if
if(.not.periodic.and.centeratom) then
  write(*,*) ' ERROR: centeratom option can only be used '
  write(*,*) ' if some periodic condition is to be applied. '
  stop
end if   

! Count the total number of frames summing up all dcd files

ntotframes = 0
do idcd = 1, ndcd
  write(*,"(tr2,52('-'))")
  open(10,file=dcdfile(idcd),action='read',form='unformatted')
  read(10) dummyc, nframes(idcd), (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
  read(10) dummyi, dummyr
  read(10) ntotat
  ntotframes = ntotframes + nframes(idcd)
  record2 = dcdfile(idcd)
  write(*,*) ' DCDfile: ', record2(1:length(record2))
  write(*,*) ' Number of atoms as specified in the dcd file: ',ntotat     
  write(*,*) ' Total number of frames in this dcd file: ', nframes(idcd)
  close(10)
end do

if(ntotframes < lastframe) then
  write(*,*) ' ERROR: lastframe greater than the number of frames dcd files. '
  stop
end if
if(lastframe == 0) lastframe = ntotframes

! Total number of frames that will actually be used

ntotframes = lastframe - firstframe + 1

! Allocate arrays that use lastatom positions

allocate(  xdcd(lastatom), ydcd(lastatom), zdcd(lastatom) )
 
! Allocate vectors for time average vector and time correlation function

allocate( xabs(ntotframes), yabs(ntotframes), zabs(ntotframes),&
          xemi(ntotframes), yemi(ntotframes), zemi(ntotframes),&
          tcf(ntotframes), &
          side(ntotframes,3), legendre(ntotframes) )

! Reading coordinates and computing 

ifcount = 0
iframe = 0
do idcd = 1, ndcd
  if ( iframe > ntotframes ) exit
  ifcount2 = 1
       
  open(10,file=dcdfile(idcd),action='read',form='unformatted')
  read(10) 
  read(10) 
  read(10) 

! Reading dcd file until the first frame, if this is the first dcd file

  if ( ifcount < firstframe - 1 ) then
    if ( ifcount + nframes(idcd) <= firstframe - 1 ) then
      ifcount = ifcount + nframes(idcd)
      close(10)
      cycle
    else
      do while( ifcount < firstframe - 1 )
        ifcount = ifcount + 1
        ifcount2 = ifcount2 + 1
        if(dcdaxis) read(10)
        read(10)
        read(10)
        read(10)
      end do
    end if
  end if

! Reading dcd file and the average vector of each frame

  do i = ifcount2, nframes(idcd)
    iframe = iframe + 1
    if ( iframe > ntotframes ) exit
  
    if(dcdaxis) then
      read(10) readsidesx, t, readsidesy, t, t, readsidesz
      side(iframe,1) = sngl(readsidesx)
      side(iframe,2) = sngl(readsidesy)
      side(iframe,3) = sngl(readsidesz)
    end if
    read(10) (xdcd(j), j = 1, lastatom)
    read(10) (ydcd(j), j = 1, lastatom)            
    read(10) (zdcd(j), j = 1, lastatom)           
  
  ! Computing absorption vector start position
   
    xt_avs = 0.
    yt_avs = 0.
    zt_avs = 0.
    do j = 1, navs
      xt_avs = xt_avs + xdcd(iavs(j))
      yt_avs = yt_avs + ydcd(iavs(j))
      zt_avs = zt_avs + zdcd(iavs(j))
    end do
    xt_avs = xt_avs / navs
    yt_avs = yt_avs / navs
    zt_avs = zt_avs / navs
  
  ! Computing the absorption vector end position
  
    xt_ave = 0.
    yt_ave = 0.
    zt_ave = 0.
    do j = 1, nave
      xt_ave = xt_ave + xdcd(iave(j))
      yt_ave = yt_ave + ydcd(iave(j))
      zt_ave = zt_ave + zdcd(iave(j))
    end do
    xt_ave = xt_ave / nave
    yt_ave = yt_ave / nave
    zt_ave = zt_ave / nave
  
  ! Thus, computing the absorption vector at this frame
  
    xabs(iframe) = xt_ave - xt_avs
    yabs(iframe) = yt_ave - yt_avs
    zabs(iframe) = zt_ave - zt_avs
    vnorm = sqrt( xabs(iframe)**2 + yabs(iframe)**2 + zabs(iframe)**2 )
    xabs(iframe) = xabs(iframe) / vnorm
    yabs(iframe) = yabs(iframe) / vnorm
    zabs(iframe) = zabs(iframe) / vnorm
  
  ! Computing the emission vector start position
  
    xt_evs = 0.
    yt_evs = 0.
    zt_evs = 0.
    do j = 1, nevs
      xt_evs = xt_evs + xdcd(ievs(j))
      yt_evs = yt_evs + ydcd(ievs(j))
      zt_evs = zt_evs + zdcd(ievs(j))
    end do
    xt_evs = xt_evs / nevs
    yt_evs = yt_evs / nevs
    zt_evs = zt_evs / nevs
  
  ! Computing the emission vector end position
  
    xt_eve = 0.
    yt_eve = 0.
    zt_eve = 0.
    do j = 1, neve
      xt_eve = xt_eve + xdcd(ieve(j))
      yt_eve = yt_eve + ydcd(ieve(j))
      zt_eve = zt_eve + zdcd(ieve(j))
    end do
    xt_eve = xt_eve / neve
    yt_eve = yt_eve / neve
    zt_eve = zt_eve / neve
  
  ! Thus, computing the emission vector at this frame
  
    xemi(iframe) = xt_eve - xt_evs
    yemi(iframe) = yt_eve - yt_evs
    zemi(iframe) = zt_eve - zt_evs
    vnorm = sqrt( xemi(iframe)**2 + yemi(iframe)**2 + zemi(iframe)**2 ) 
    xemi(iframe) = xemi(iframe) / vnorm
    yemi(iframe) = yemi(iframe) / vnorm
    zemi(iframe) = zemi(iframe) / vnorm
  
  ! Computing emission vector from a rotation of the absoprtion vector, if desired (theta>0)
  
    if ( theta > 1.d-7 ) then
  
      ! Rotation axis is perpendicular to both vectors (external product)
  
      ax = yabs(iframe)*zemi(iframe) - zabs(iframe)*yemi(iframe)
      ay = zabs(iframe)*xemi(iframe) - xabs(iframe)*zemi(iframe)
      az = xabs(iframe)*yemi(iframe) - yabs(iframe)*xemi(iframe)
      vnorm = sqrt( ax**2 + ay**2 + az**2 )
      ax = ax / vnorm
      ay = ay / vnorm
      az = az / vnorm
  
      ! Rotation matrix (quaternion derived) 
  
      xt = 1. - cost
      r(1,1) = cost + ax*ax*xt ; r(1,2) = ax*ay*xt-az*sint; r(1,3) = ax*az*xt+ay*sint
      r(2,1) = ay*ax*xt+az*sint; r(2,2) = cost+ay*ay*xt   ; r(2,3) = ay*az*xt-ax*sint
      r(3,1) = az*ax*xt-ay*sint; r(3,2) = az*ay*xt+ax*sint; r(3,3) = cost + az*az*xt
  
      ! Applying rotation matrix
  
      xt = xabs(iframe)
      yt = yabs(iframe)
      zt = zabs(iframe)
      xemi(iframe) = r(1,1)*xt + r(1,2)*yt + r(1,3)*zt
      yemi(iframe) = r(2,1)*xt + r(2,2)*yt + r(2,3)*zt
      zemi(iframe) = r(3,1)*xt + r(3,2)*yt + r(3,3)*zt
  
    end if
  
    if ( iframe == 1 ) then
      write(*,*)
      write(*,"( tr2,52('-') )")
      write(*,*) ' At first frame, the normalized absorption vector is: '
      write(*,"('  {',3(f8.3),'} {',3(f8.3),' }')") xt_avs, yt_avs, zt_avs,&
                   xt_avs+xabs(1), yt_avs+yabs(1), zt_avs+zabs(1)
      write(*,*) ' At first frame, the normalized emission vector is: '
      write(*,"('  {',3(f8.3),'} {',3(f8.3),' }')") xt_evs, yt_evs, zt_evs,&
                   xt_avs+xemi(1), yt_avs+yemi(1), zt_avs+zemi(1)
      write(*,*)
      write(*,*) ' You can use these points to draw lines in vmd with: '
      write(*,*) ' draw line { 0 0 0 } { 1 1 1 } '
      write(*,"( tr2,52('-') )")
      write(*,*)
    end if

  ! Print progress
  
    if ( i == ifcount2 ) then
      write(*,"( /,'  Reading DCD file ',i4,' and computing vectors: ',&
                  f6.2,'%' )", advance="no") idcd, 0.
    else
      write(*,"( 7a,f6.2,'%' )",advance='no')&
           (char(8),j=1,7), 100.*float(iframe)/float(ntotframes)
    end if
  
  end do
  close(10)
end do

! Computing the time-dependent correlation function

do i = 1, ntotframes
  tcf(i) = 0.
  legendre(i) = 0.
end do

do i = 1, ntotframes
  do j = i,  ntotframes

! Computing the internal product of absoprtion and emission vectors

    int_prod = xabs(i)*xemi(j) + yabs(i)*yemi(j) + zabs(i)*zemi(j)
   
! Computing the time correlation function and its second legendre polynomial

    tcf(j-i+1) = tcf(j-i+1) + int_prod
    legendre(j-i+1) = legendre(j-i+1) + 0.5*( 3. * int_prod**2 - 1. )

  end do
end do

! Open output file and writes all information of this run

open(20,file=output(1:length(output)))
write(20,"( '#',/,&
           &'# Output of tcf.f90:',/,&
           &'# Input file: ',a )") inputfile(1:length(inputfile))
do idcd = 1, ndcd
  record2 = dcdfile(idcd)
  write(20,"( '# DCD file: ',i5,a )") idcd, record2(1:length(record2)) 
end do
write(20,"('# First frame: ',i5,' Last frame: ',i5,/,&
           &'#',/,&
           &'# Periodic conditions: ',/,&
           &'# Periodic: ',l1,' Read from DCD: ',l1,&
           &   ' Center atom: ',l1,' ',i6,/,&
           &'#',/,'#   FRAME DELTA          2ND LEG POL',&
           &'      INTER PROD')")& 
           &firstframe, lastframe,&
           &periodic, readfromdcd, centeratom, catom

! Rescaling the tcf data and printing results  

do i = 1, ntotframes - 1
  tcf(i) = tcf(i) / float( ntotframes - i + 1 )
  legendre(i) = legendre(i) / float( ntotframes - i + 1 )
  write(20,"( t4, f12.6, t26, f8.5, t43, f8.5 )") scaletime*(i-1), r0*legendre(i), tcf(i)
end do
close(20)

! Write final messages with names of output files and their content

write(*,*)
write(*,"( tr2,52('-') )")
write(*,*)
write(*,*) ' OUTPUT FILES: ' 
write(*,*)
write(*,*) ' Wrote output file: ', output(1:length(output))
write(*,*)
write(*,*) ' Containing the time correlation function of the '
write(*,*) ' average vector. '
write(*,*) 
write(*,*) '####################################################'
write(*,*) 
write(*,*) '  END: Normal termination.  '
write(*,*) 
write(*,*) '####################################################'
write(*,*)        

end

