!
! hbonds: Computes the number of hydrogen bonds between two groups
!         of atoms (does not consider the nature of the groups, 
!         that is, donnors and acceptors, just geometrical properties
!         of the interaction and the presence of the hydrogen atoms).
!
! Auxiliar dimensions:
!          memory: Controls the amount of memory used for reading dcd
!                  files. If the program ends with segmentation fault
!                  without any other printing, decrease the value of
!                  this parameter.  
!
! L. Martinez, Mar 17, 2008. (modified Aug, 29, 2014 for linked cells)
!
! Dynamically allocatable arrays and their dimensions:
!
!       group1: number of atoms of group1
!       group2: number of atoms of group2
!       eps, sig, q, e, s, class, segat, resat, typeat, classat: total
!               number of atoms of the system
!

! Modules

module hbonds_linkedcells

  use charsize
  implicit none
  integer, parameter :: memory = 15000000 
  real, parameter :: pi=3.141592654
  real :: distance, xdcd(memory), ydcd(memory), zdcd(memory), coshbond, axis(3)
  integer :: nboxes(3), ngroup1, ngroup2
  integer, allocatable :: iatfg1(:,:,:), iatfg2(:,:,:), iatnextg1(:), &
                          iatnextg2(:), group1(:), group2(:), h_bond_first(:), &
                          h_bond_next(:)
  logical :: lifetime
  logical, allocatable :: is_hydrogen(:), g12(:)
  character(len=200) :: lifetime_file
  character(len=charsize1), allocatable :: typeat(:)

end module hbonds_linkedcells

program hbonds

  ! Static variables
  
  use hbonds_linkedcells
  use charsize
  implicit none
  integer :: natom, nh1, nh2, nbonds,&
             narg, length, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
             j, iframe1, iframe2, icycle, nfrcycle, iatom, k, ii, &
             status, keystatus, ndim, iargc, lastatom, nhtot, kframe,&
             ibox, jbox, kbox, ig1, nhparc, nbdim(3), ib1, ib2
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: side(memory,3)
  real :: dummyr, time0, etime, tarray(2), angle, x1, y1, z1, &
          scaletime, xmin(3), xmax(3), &
          boxlength, dbox_x, dbox_y, dbox_z
  character(len=200) :: groupfile, line, record, value, keyword, dcdfile,&
                        inputfile, output, psffile
  character(len=4) :: dummyc
  logical :: periodic, readfromdcd, dcdaxis
  
  ! Allocatable arrays
  
  integer, allocatable :: resid(:)
  real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:)
  character(len=charsize1), allocatable :: class(:), segat(:), resat(:), classat(:)
  character(len=charsize1) :: type1, type2
  
  ! Compute time
  
  time0 = etime(tarray)
  
  ! Output title
  
  write(*,"(/,' ####################################################',&
            &/,/,&
            & '   HBONDS: Compute H-bonds from DCD files.           ',&
            &/,/,&
            & ' ####################################################',/)")    
  
  call version()
  
  ! Some default parameters
  
  firstframe = 1
  lastframe = 0
  stride = 1
  periodic = .true.
  readfromdcd = .true.
  angle = 20.
  distance = 3.
  scaletime = 1.
  lifetime = .false.
  
  ! Open input file and read parameters
  
  narg = iargc()
  if(narg == 0) then
    write(*,*) ' Run with: ./hbonds input.inp '
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
    else if(keyword(record) == 'lifetime') then
      lifetime = .true.
      lifetime_file = value(record)
      write(*,*) ' Will write lifetime file: ', psffile(1:length(psffile))
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
    else if(keyword(record) == 'distance') then
      line = value(record)
      read(line,*,iostat=keystatus) distance
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'angle') then
      line = value(record)
      read(line,*,iostat=keystatus) angle
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'output') then
      output = value(record)
      write(*,*) ' Output file name: ', output(1:length(output))
    else if(keyword(record) == 'group1' .or. &
            keyword(record) == 'group2') then
      write(*,"(a,/,a)") ' ERROR: The options group1 and group2 must be used ',&
                         '        with the hbonds.sh script, not directly. '
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
  
  ! File for lifetime data computing
  
  if ( lifetime ) then
    open(111,file=lifetime_file)
    write(111,"( '# System: ', a )") psffile(1:length(psffile))
    write(111,"( '# Life-time data for DCD file: ', a )") dcdfile(1:length(dcdfile))  
  end if
  
  ! Reading the header of psf file and parameter files to get dimensions
  
  call getdim(psffile,inputfile,ndim)
  allocate( eps(ndim), sig(ndim), q(ndim), e(ndim), s(ndim),&
            segat(ndim), resat(ndim), classat(ndim), typeat(ndim),&
            class(ndim), mass(ndim), resid(ndim) )
         
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
  
  ! Read PSF file to get charges and assign parameters
  
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
      if(line(57:60) == '1.00') ngroup1 = ngroup1 + 1        
      if(line(63:66) == '2.00') ngroup2 = ngroup2 + 1        
    end if     
  end do
  allocate ( group1(ngroup1), group2(ngroup2), iatnextg1(ngroup1), iatnextg2(ngroup2),&
             h_bond_first(natom), h_bond_next(natom), is_hydrogen(natom), g12(ngroup1) )
  close(10)
  
  ! Now reading reading the group atoms
  
  open(10,file=groupfile,action='read')
  ngroup1 = 0
  ngroup2 = 0
  iatom = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      iatom = iatom + 1
  
  ! Read atoms belonging to group1
  
      if(line(57:60) == '1.00') then      
        ngroup1 = ngroup1 + 1        
        group1(ngroup1) = iatom
      end if
  
  ! Read atoms belonging to group2
  
     if(line(63:66) == '2.00') then
        ngroup2 = ngroup2 + 1        
        group2(ngroup2) = iatom
      end if
  
    end if     
  end do
  close(10)
  lastatom = max0(group1(ngroup1),group2(ngroup2))

  ! Check if atoms belong to both groups, to escape from repeted counting

  do i = 1, ngroup1
    g12(i) = .false.
    do j = 1, ngroup2
      if ( group1(i) == group2(j) ) then
        g12(i) = .true.
      end if
    end do
  end do
  
  ! Check which atoms of these groups are bound to hydrogen atoms
  
  do i = 1, natom
    h_bond_first(i) = 0 
    is_hydrogen(i) = .false.
    type1 = typeat(i)
    if ( type1(1:1) == 'H' ) is_hydrogen(i) = .true.
  end do
  open(10,file=psffile,status='old',action='read')
  do
    read(10,*,iostat=status) nbonds, line
    if ( status /= 0 ) cycle
    if ( line == '!NBOND:' .or. line == '!NBONDS' ) then
      i = 0
      do while( i < nbonds )
        read(10,"( a200 )",iostat=status) record
        j = 1
        read(record,*,iostat=status) (ib1, ib2, k = 1, j)
        do while( status == 0 ) 
          type1 = typeat(ib1)
          type2 = typeat(ib2)
          if ( type1(1:1) == 'H' .and. type2(1:1) /= 'H' ) then
            h_bond_next(ib1) = h_bond_first(ib2)
            h_bond_first(ib2) = ib1
          end if
          if ( type2(1:1) == 'H' .and. type1(1:1) /= 'H' ) then
            h_bond_next(ib2) = h_bond_first(ib1)
            h_bond_first(ib1) = ib2
          end if
          j = j + 1
          read(record,*,iostat=status) (ib1, ib2, k = 1, j)
        end do
        i = i + j - 1
      end do
    else
      cycle
    end if
    exit
  end do
  close(10)
  
  ! Output some group properties for testing purposes
  
  write(*,*) ' Number of atoms of group 1: ', ngroup1 
  write(*,*) ' First atom of group 1: ', group1(1)
  write(*,*) ' Last atom of group 1: ', group1(ngroup1)
  write(*,*) ' Number of atoms of group 2: ', ngroup2 
  write(*,*) ' First atom of group 2: ', group2(1)
  write(*,*) ' Last atom of group 2: ', group2(ngroup2)
  if(ngroup1 == 0 .or. ngroup2 == 0) then
    write(*,*) ' ERROR: Groups must contain at least one atom.'
    stop
  end if
   
  ! Counting the number of hydrogen atoms of each group
  
  nh1 = 0
  do i = 1, ngroup1
    record = typeat(group1(i))
    if(record(1:1).eq.'H') nh1 = nh1 + 1
  end do
  nh2 = 0
  do i = 1, ngroup2
    record = typeat(group2(i))
    if(record(1:1).eq.'H') nh2 = nh2 + 1
  end do
  write(*,*) ' Number of hydrogen atoms found in group1: ', nh1
  write(*,*) ' Number of hydrogen atoms found in group2: ', nh2
  if(nh1 == 0 .and. nh2 == 0) then
    write(*,*) ' ERROR: Could not find H atoms in any of the groups. '
    stop
  end if
   
  ! Checking if dcd file contains periodic cell information
  
  write(*,"( /,tr2,52('-') )")
  write(*,*) ' Periodic cell data: Read carefully. '
  call chkperiod(dcdfile,dcdaxis,readfromdcd) 
  if(.not.readfromdcd.and.periodic) then
    write(*,*) ' User provided periodic cell dimensions: '
    write(*,*) axis(1), axis(2), axis(3)
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
             &'# Output of hbonds.f:',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ',a,/,& 
             &'# First frame: ',i5,' Last frame: ',i5,' Stride: ',i5,/,&
             &'#',/,&
             &'# Periodic conditions: ',/,&
             &'# Periodic: ',l1,' Read from DCD: ',l1,&
             &'#',/,&
             &'# Number of atoms of group 1: ', i6,/,&
             &'# First and last atoms of group 1: ',i6,tr1,i6,/,&
             &'# Number of atoms of group 2: ', i6,/,&
             &'# First and last atoms of group 2: ',i6,tr1,i6,/,&
             &'#',/,&
             &'# Number of hydrogens in group 1: ', i6,/,&
             &'# Number of hydrogens in group 2: ', i6,/,&
             &'#',/,&
             &'# Angle cutoff criterium (degrees): ', f8.3,/,&
             &'# Distance cutoff criterium (A): ', f8.3,/,&
             &'#',/,&
             &'# Number of hydrogen bonds formed between the groups:',/,&
             &'#  TIME ',t13,'NUMBER')")& 
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, &
             &ngroup1, group1(1), group1(ngroup1),& 
             &ngroup2, group2(1), group2(ngroup2),&
             &nh1,nh2,angle,distance
   
  ! Convert angle to its cosine
  
  coshbond = cos(pi*angle/180)
  
  ! Now going to read the dcd file
  
  memframes = memory / ntotat
  ncycles = lastframe / memframes + 1
  memlast = lastframe - memframes * ( ncycles - 1 )
  write(*,*) ' Will read and store in memory at most', memframes,&
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
  
  do i = 1, 3
    nbdim(i) = 0
  end do
  allocate( iatfg1(1,1,1), iatfg2(1,1,1) )
  
  ! Reading dcd file and computing h-bonds
  
  iframe1 = 0
  iframe2 = 0
  do icycle = 1, ncycles 
  
    write(*,"( t3,i4,' Read...',t17,'|' )",advance='no') icycle
  
    ! Each cycle fills the memory as specified by the memory parameter 
  
    if(icycle == ncycles) then
      nfrcycle = memlast
    else
      nfrcycle = memframes
    end if
  
    iatom = 0
    do kframe = 1, nfrcycle    
      iframe1 = iframe1 + 1
      if(dcdaxis) then
        read(10) readsidesx, t, readsidesy, t, t, readsidesz
        side(kframe,1) = sngl(readsidesx)
        side(kframe,2) = sngl(readsidesy)
        side(kframe,3) = sngl(readsidesz)
      end if
      read(10) (xdcd(k), k = iatom + 1, iatom + lastatom)
      read(10) (ydcd(k), k = iatom + 1, iatom + lastatom)            
      read(10) (zdcd(k), k = iatom + 1, iatom + lastatom)           
  
      if(mod(iframe1,stride) /= 0 .or. iframe1 < firstframe) then
        iatom = iatom + ntotat
        cycle
      end if
  
      ! Putting the atoms in their minimum image coordinates if using periodic
      ! boundary conditions
    
      if ( periodic ) then
        if ( readfromdcd ) then
          axis(1) = side(kframe,1)
          axis(2) = side(kframe,2)
          axis(3) = side(kframe,3)
        end if
        do i = 1, ngroup1
          ii = iatom + group1(i)
          x1 = xdcd(ii)
          y1 = ydcd(ii)
          z1 = zdcd(ii)
          call image(x1,y1,z1,axis(1),axis(2),axis(3))
          xdcd(ii) = x1
          ydcd(ii) = y1
          zdcd(ii) = z1
        end do
        do i = 1, ngroup2
          ii = iatom + group2(i)
          x1 = xdcd(ii)
          y1 = ydcd(ii)
          z1 = zdcd(ii)
          call image(x1,y1,z1,axis(1),axis(2),axis(3))
          xdcd(ii) = x1
          ydcd(ii) = y1
          zdcd(ii) = z1
        end do
      end if
  
      ! Updating counter 
  
      iatom = iatom + ntotat
      if(mod(kframe,max0(1,nfrcycle/35)) == 0) write(*,"( '*' )",advance='no') 
    end do
    write(*,"( '|' )")   
    write(*,"( t3,'Computing... ',t17,'|' )",advance='no') 
  
    ! Computing the hydrogen bonds
  
    iatom = 0
    do kframe = 1, nfrcycle
      iframe2 = iframe2 + 1
  
      if(mod(iframe2,stride) /= 0 .or. iframe2 < firstframe) then
        iatom = iatom + ntotat
        cycle
      end if
  
      ! Prepare the linked cells
  
      if ( periodic ) then
        if ( readfromdcd ) then
          axis(1) = side(kframe,1)
          axis(2) = side(kframe,2)
          axis(3) = side(kframe,3)
        end if
        xmin(1) = -axis(1)/2.
        xmin(2) = -axis(2)/2.
        xmin(3) = -axis(3)/2.
        xmax(1) = axis(1)/2.
        xmax(2) = axis(2)/2.  
        xmax(3) = axis(3)/2.  
      else 
        do i = 1, 3
          xmin(i) = 1.e30
          xmax(i) = -1.e30
        end do
        do i = 1, ngroup1
          ii = iatom + group1(i)
          xmin(1) = min(xmin(1),xdcd(ii))
          xmin(2) = min(xmin(2),ydcd(ii))
          xmin(3) = min(xmin(3),zdcd(ii))
          xmax(1) = max(xmax(1),xdcd(ii))
          xmax(2) = max(xmax(2),ydcd(ii))
          xmax(3) = max(xmax(3),zdcd(ii))
        end do
        do i = 1, ngroup2
          ii = iatom + group2(i)
          xmin(1) = min(xmin(1),xdcd(ii))
          xmin(2) = min(xmin(2),ydcd(ii))
          xmin(3) = min(xmin(3),zdcd(ii))
          xmax(1) = max(xmax(1),xdcd(ii))
          xmax(2) = max(xmax(2),ydcd(ii))
          xmax(3) = max(xmax(3),zdcd(ii))
        end do
      end if
  
      ! The side of the linked cells must be an exact divisor of the
      ! box side, because of the phantom replicas
  
      boxlength = xmax(1) - xmin(1)
      nboxes(1) = max(1,int(boxlength/distance))
      dbox_x = boxlength / nboxes(1)
  
      boxlength = xmax(2) - xmin(2)
      nboxes(2) = max(1,int(boxlength/distance))
      dbox_y = boxlength / nboxes(2)
  
      boxlength = xmax(3) - xmin(3)
      nboxes(3) = max(1,int(boxlength/distance))
      dbox_z = boxlength / nboxes(3)
  
      ! If the number of cells changed from the previous frame, update dimensions
  
      if ( nboxes(1) > nbdim(1) .or. &
           nboxes(2) > nbdim(2) .or. &
           nboxes(3) > nbdim(3) ) then
        deallocate( iatfg1, iatfg2 )
        nbdim(1) = nboxes(1)
        nbdim(2) = nboxes(2)
        nbdim(3) = nboxes(3)
        allocate( iatfg1(0:nbdim(1)+1,0:nbdim(2)+1,0:nbdim(3)+1), &
                  iatfg2(0:nbdim(1)+1,0:nbdim(2)+1,0:nbdim(3)+1) )
      end if
  
      ! Reset the linked list initial atom vector
  
      do i = 0, nboxes(1) + 1
        do j = 0, nboxes(2) + 1
          do k = 0, nboxes(3) + 1
            iatfg1(i,j,k) = 0
            iatfg2(i,j,k) = 0
          end do
        end do
      end do
  
      ! Putting the atoms in their boxes
  
      do i = 1, ngroup1
        if ( is_hydrogen(group1(i)) ) cycle
        ii = iatom + group1(i)
        ibox = int( (xdcd(ii) - xmin(1)) / dbox_x ) + 1
        jbox = int( (ydcd(ii) - xmin(2)) / dbox_y ) + 1
        kbox = int( (zdcd(ii) - xmin(3)) / dbox_z ) + 1
        if ( ibox == nbdim(1)+1 ) ibox = nbdim(1)
        if ( jbox == nbdim(2)+1 ) jbox = nbdim(2)
        if ( kbox == nbdim(3)+1 ) kbox = nbdim(3)
        iatnextg1(i) = iatfg1(ibox,jbox,kbox)
        iatfg1(ibox,jbox,kbox) = i
      end do
      do i = 1, ngroup2
        if ( is_hydrogen(group2(i)) ) cycle
        ii = iatom + group2(i)
        ibox = int( (xdcd(ii) - xmin(1)) / dbox_x ) + 1
        jbox = int( (ydcd(ii) - xmin(2)) / dbox_y ) + 1
        kbox = int( (zdcd(ii) - xmin(3)) / dbox_z ) + 1
        if ( ibox == nbdim(1)+1 ) ibox = nbdim(1)
        if ( jbox == nbdim(2)+1 ) jbox = nbdim(2)
        if ( kbox == nbdim(3)+1 ) kbox = nbdim(3)
        iatnextg2(i) = iatfg2(ibox,jbox,kbox)
        iatfg2(ibox,jbox,kbox) = i
      end do
  
      ! Fill phantom boxes for periodic boundary conditions with group2
  
      if ( periodic ) call phantomcells(nboxes,iatfg2, nbdim)
  
      ! Now the loop is over boxes

      nhtot = 0
      do i = 1, nboxes(1)
        do j = 1, nboxes(2) 
          do k = 1, nboxes(3)
   
            ig1 = iatfg1(i,j,k)
            do while( ig1 /= 0 ) 
  
              ! Interactions inside box
  
              nhtot = nhtot + nhparc(ig1,i,j,k,iatom,iframe2)
  
              ! Interactions of boxes that share faces
  
              nhtot = nhtot + nhparc(ig1,i+1,j,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j+1,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j,k+1,iatom,iframe2)
  
              nhtot = nhtot + nhparc(ig1,i-1,j,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j-1,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j,k-1,iatom,iframe2)
  
              ! Interactions of boxes that share axes
  
              nhtot = nhtot + nhparc(ig1,i+1,j+1,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i+1,j,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i+1,j-1,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i+1,j,k-1,iatom,iframe2)
  
              nhtot = nhtot + nhparc(ig1,i,j+1,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j+1,k-1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j-1,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i,j-1,k-1,iatom,iframe2)
  
              nhtot = nhtot + nhparc(ig1,i-1,j+1,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i-1,j,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i-1,j-1,k,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i-1,j,k-1,iatom,iframe2)
  
              ! Interactions of boxes that share vertices
  
              nhtot = nhtot + nhparc(ig1,i+1,j+1,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i+1,j+1,k-1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i+1,j-1,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i+1,j-1,k-1,iatom,iframe2)
  
              nhtot = nhtot + nhparc(ig1,i-1,j+1,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i-1,j+1,k-1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i-1,j-1,k+1,iatom,iframe2)
              nhtot = nhtot + nhparc(ig1,i-1,j-1,k-1,iatom,iframe2)
  
              ig1 = iatnextg1(ig1)
            end do
          end do
        end do
      end do
  
      ! Printing to output file in this frame 
   
      write(51,"( e12.6,tr2,i8 )") scaletime*iframe2, nhtot
  
      ! Printing status bar
  
      if(mod(kframe,max0(1,nfrcycle/35)) == 0) write(*,"( '*' )",advance='no') 
  
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
  write(*,*) ' Which contains the total interaction energies of '
  write(*,*) ' the groups as a function of time.'
  if ( lifetime ) then
    close(111)
    write(*,*) 
    write(*,*) ' Wrote lifetime data file: ', lifetime_file(1:length(lifetime_file))
  end if
  write(*,*) 
  write(*,*) ' Running time: ', time0
  write(*,*) '####################################################'
  write(*,*) 
  write(*,*) '  END: Normal termination.  '
  write(*,*) 
  write(*,*) '####################################################'
  write(*,*)        

end program hbonds

!
! Function that checks if there is a hydrogen bond between the atom
! of group1 which is input and the atoms of group2 of defined by the
! first atom of group2 input
!

integer function nhparc(ig1,ibox,jbox,kbox,iatom,iframe)

  use charsize
  use hbonds_linkedcells
  implicit none
  real :: cosangle, sij, xh, yh, zh, xnorm, cosang, x1, y1, z1, x2, y2, z2
  integer :: ig1, ig2, ii, jj, iatom, iframe, jh, ibox, jbox, kbox

  nhparc = 0
  ig2 = iatfg2(ibox,jbox,kbox)
  ii = iatom + group1(ig1)
  x1 = xdcd(ii)
  y1 = ydcd(ii)
  z1 = zdcd(ii)
  ig2do : do while( ig2 /= 0 )

    ! If none of the atoms is bound to any hydrogen, cycle

    if( h_bond_first(group1(ig1)) == 0 .and. h_bond_first(group2(ig2)) == 0 ) then 
      ig2 = iatnextg2(ig2)
      cycle ig2do
    end if

    ! If the atoms belong to both groups, cycle if the computation is not
    ! needed

    if ( g12(ig1) ) then
      if ( group1(ig1) >= group2(ig2) ) then
        ig2 = iatnextg2(ig2)
        cycle ig2do
      end if
    end if

    ! Computing the distance between ii and jj (atom jj is in the origin)

    jj = iatom + group2(ig2)
    call movephantomcoor(xdcd(jj),ydcd(jj),zdcd(jj),x2,y2,z2,&
                         nboxes,ibox,jbox,kbox,axis)

    sij = xnorm(x1-x2,y1-y2,z1-z2)

    ! Cycle if the distance does not satisfy the h-bond criterium

    if ( sij >= distance ) then
      ig2 = iatnextg2(ig2)
      cycle ig2do
    end if

    ! Check if there is a hydrogen bond with the hydrogens bound to ig1

    jh = h_bond_first(group1(ig1)) 
    do while( jh /= 0 ) 

      ! Fix the coordinate of the hydrogen, which may be on the 
      ! the other side of the box
      
      xh = xdcd(iatom+jh) - x1
      yh = ydcd(iatom+jh) - y1
      zh = zdcd(iatom+jh) - z1
      call image(xh,yh,zh,axis(1),axis(2),axis(3))
      xh = xh + x1
      yh = yh + y1
      zh = zh + z1

      ! Compute angle and count h-bond if the criterium is satified

      cosangle = cosang(x2-xh,y2-yh,z2-zh,xh-x1,yh-y1,zh-z1)
      if ( cosangle > coshbond ) then
        nhparc = nhparc + 1
        if ( lifetime ) write(111,*) iframe, group1(ig1), group2(ig2)
        ig2 = iatnextg2(ig2)
        cycle ig2do
      end if
      jh = h_bond_next(jh)
    end do

    ! Check if there is a hydrogen bond with the hydrogens bound to ig1

    jh = h_bond_first(group2(ig2)) 
    do while( jh /= 0 ) 

      ! Fix the coordinate of the hydrogen, which may be on the 
      ! the other side of the box
      
      xh = xdcd(iatom+jh) - x2
      yh = ydcd(iatom+jh) - y2
      zh = zdcd(iatom+jh) - z2
      call image(xh,yh,zh,axis(1),axis(2),axis(3))
      xh = xh + x2
      yh = yh + y2
      zh = zh + z2

      ! Compute angle and count h-bond if the criterium is satified

      cosangle = cosang(x2-xh,y2-yh,z2-zh,xh-x1,yh-y1,zh-z1)
      if ( cosangle > coshbond ) then
        nhparc = nhparc + 1
        if ( lifetime ) write(111,*) iframe, group1(ig1), group2(ig2)
        ig2 = iatnextg2(ig2)
        cycle ig2do
      end if
      jh = h_bond_next(jh)
    end do

    ig2 = iatnextg2(ig2)
  end do ig2do

end function nhparc
