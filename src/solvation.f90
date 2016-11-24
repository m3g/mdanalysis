!
! solvation: A program for computing the solvation of selected groups
!            in molecular dynamics simulations in NAMD DCD format.
!
! Auxiliar dimensions:
!          memory: Controls the amount of memory used for reading dcd
!                  files. If the program ends with segmentation fault
!                  without any other printing, decrease the value of
!                  this parameter.  
!
! L. Martinez, Apr 22, 2008 (modified Sep 3, 2014 for linked cells).
!
! Dynamically allocatable arrays and their dimensions:
!
!       solute: number of atoms of solute
!       solvent: number of atoms of solvent
!       eps, sig, q, e, s, class, segat, resat, typeat, classat: total
!               number of atoms of the system
!

! Module used for linked cell method

module linkedcells_solvation

  integer, parameter :: memory = 15000000
  real :: xdcd(memory), ydcd(memory), zdcd(memory)
  real :: axis(3)

  integer :: nsolvlist, nboxes(3)
  integer, allocatable :: iatomfirst(:,:,:), iatomnext(:), resid(:),&
                          solute(:), solvent(:)

  real :: shellsize2
  logical :: periodic

end module linkedcells_solvation

! Main program

program solvation

  ! Static variables
 
  use charsize
  use linkedcells_solvation
  implicit none
  integer :: natom, narg, length, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
             j, iframe, icycle, nfrcycle, iatom, k, ii, &
             status, keystatus, ndim, iargc, lastatom, nrsolvent, &
             nrsolute, kframe, nsolv, nsolute, nsolvent, nrsolv, nres
  integer :: nbdim(3), ibox, jbox, kbox, isolute
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: mass1, mass2, side(memory,3)
  real :: dummyr, x1, y1, z1, time0, etime, tarray(2),&
          scaletime, xmin(3), xmax(3), shellsize,&
          boxlength, dbox_x, dbox_y, dbox_z
  character(len=200) :: groupfile, line, record, value, keyword,&
                        dcdfile, inputfile, output, psffile, perresidue_output
  character(len=4) :: dummyc
  logical :: readfromdcd, dcdaxis, perresidue
  
  ! Allocatable arrays
  
  integer, allocatable :: solv_list(:,:)
  integer, allocatable :: natres(:), fatres(:), fatrsolute(:), fatrsolvent(:),&
                          solute_box(:,:), solv_per_residue(:), rseqsolute(:)
  real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:)
  character(len=charsize1), allocatable :: class(:), segat(:), resat(:),&
                                           typeat(:), classat(:)
  logical, allocatable :: tcount(:), counted(:), soluteatom(:), solventatom(:), &
                          rsolv(:)
  
  ! Compute time
  
  time0 = etime(tarray)
  
  ! Output title
  
  write(*,"(/,' ####################################################',&
            &/,/,&
            & '   SOLVATION: Compute solvation from DCD files  ',&
            &/,/,&
            & ' ####################################################',/)")    
  
  call version()
  
  ! Some default parameters
  
  firstframe = 1
  lastframe = 0
  stride = 1
  periodic = .true.
  readfromdcd = .true.
  shellsize = 0.d0
  scaletime = 1.
  nsolvlist = 0
  perresidue = .false.
  
  ! Open input file and read parameters
  
  narg = iargc()
  if(narg == 0) then
    write(*,*) ' Run with: ./solvation input.inp '
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
    else if(keyword(record) == 'nsolvlist') then
      line = value(record)
      read(line,*,iostat=keystatus) nsolvlist
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'shellsize') then
      line = value(record)
      read(line,*,iostat=keystatus) shellsize
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
    else if(keyword(record) == 'output') then
      output = value(record)
      write(*,*) ' Output file name: ', output(1:length(output))
    else if(keyword(record) == 'perresidue') then
      perresidue = .true.
      perresidue_output = value(record)
      write(*,*) ' Output file for solvation per residue: ',&
                 perresidue_output(1:length(perresidue_output))
    else if(keyword(record) == 'solute' .or. &
            keyword(record) == 'solvent') then
      write(*,"(a,/,a)") ' ERROR: The options solute and solvent must be used ',&
                         '        with the solvation.sh script, not directly. '
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
  
  ! Reading the header of psf file
  
  call getdim(psffile,inputfile,ndim)
  allocate( eps(ndim), sig(ndim), q(ndim), e(ndim), s(ndim), mass(ndim),&
            segat(ndim), resat(ndim), classat(ndim), typeat(ndim), class(ndim),&
            resid(ndim), natres(ndim), fatres(ndim), fatrsolute(ndim), rseqsolute(ndim),&
            fatrsolvent(ndim), tcount(ndim), counted(ndim), &
            soluteatom(ndim), solventatom(ndim))

  ! Check for simple input errors
  
  if(stride < 1) then
    write(*,*) ' ERROR: stride cannot be less than 1. ' 
    stop
  end if
  if(lastframe < firstframe.and.lastframe /= 0) then
    write(*,*) ' ERROR: lastframe must be greater or equal to firstframe. '
    stop
  end if
  if(shellsize == 0.d0) then
    write(*,*) ' ERROR: shellsize not defined or zero. '
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
  nres = resid(natom)
  write(*,*) ' Number of atoms in PSF file: ', natom
  write(*,*) ' Number of residues in PSF file: ', nres
  write(*,*)
   
  ! Read solute and solvent information from file
  ! First reading the size of the groups to allocate arrays
  
  open(10,file=groupfile,action='read',status='old')
  nsolute = 0
  nsolvent = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      if(line(63:66) == '1.00') then      
        nsolute = nsolute + 1        
      else if(line(63:66) == '2.00') then
        nsolvent = nsolvent + 1        
      end if
    end if     
  end do
  if( nsolute == 0 ) then
    write(*,*) ' ERROR: Could not find solute atoms. '
    stop
  end if
  if( nsolvent == 0 ) then
    write(*,*) ' ERROR: Could not find solvent atoms. '
    stop
  end if
  allocate ( solute(nsolute), solvent(nsolvent), solute_box(nsolute,3) )
  
  ! Now reading the group atoms
  
  rewind(10)
  nsolute = 0
  nsolvent = 0
  mass1 = 0.d0
  mass2 = 0.d0
  iatom = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      iatom = iatom + 1
      soluteatom(iatom) = .false.
      solventatom(iatom) = .false.
      if(line(63:66) == '1.00') then      
  
        ! Read atoms belonging to solute
  
        nsolute = nsolute + 1        
        solute(nsolute) = iatom
        soluteatom(iatom) = .true.
        mass1 = mass1 + mass(iatom)
  
      else if(line(63:66) == '2.00') then
  
        ! Read atoms belonging to solvent
  
        nsolvent = nsolvent + 1        
        solvent(nsolvent) = iatom
        solventatom(iatom) = .true.
        mass2 = mass2 + mass(iatom)
  
      end if
    end if     
  end do
  close(10)
  lastatom = max0(solute(nsolute),solvent(nsolvent))
  
  ! Finding the number of atoms and the first atom of each residue
  
  j = 0
  do i = 1, natom
    if(resid(i).gt.j) then
      fatres(resid(i)) = i
      natres(resid(i)) = 1
      j = resid(i)
    else
      natres(resid(i)) = natres(resid(i)) + 1
    end if
  end do
  
  ! Counting the number of residues of the solute and solvent
  
  j = 0
  nrsolute = 0
  do i = 1, nsolute
    if(resid(solute(i)).gt.j) then
      nrsolute = nrsolute + 1 
      rseqsolute(nrsolute) = resid(solute(i))
      fatrsolute(nrsolute) = fatres(resid(solute(i)))
      j = resid(solute(i))
    end if
  end do
  j = 0
  nrsolvent = 0
  do i = 1, nsolvent
    if(resid(solvent(i)).gt.j) then
      nrsolvent = nrsolvent + 1 
      fatrsolvent(nrsolvent) = fatres(resid(solvent(i)))
      j = resid(solvent(i))
    end if
  end do
  if ( nsolvlist == 0 ) nsolvlist = 2*nres
  allocate(solv_list(nsolvlist,2),rsolv(nres),solv_per_residue(nrsolute))

  ! Output some group properties for testing purposes
  
  write(*,*) ' Number of atoms of solute: ', nsolute 
  write(*,*) ' First atom of solute: ', solute(1)
  write(*,*) ' Last atom of solute: ', solute(nsolute)
  write(*,*) ' Number of residues in solute: ', nrsolute
  write(*,*) ' Mass of solute: ', mass1
  write(*,*) ' Number of atoms of solvent: ', nsolvent 
  write(*,*) ' First atom of solvent: ', solvent(1)
  write(*,*) ' Last atom of solvent: ', solvent(nsolvent)
  write(*,*) ' Number of residues in solvent: ', nrsolvent
  write(*,*) ' Mass of solvent: ', mass2
  if(nsolute == 0 .or. nsolvent == 0) then
    write(*,*) ' ERROR: Both solute and solvent must contain at least one atom.'
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
  
  ! Open output file for total solvation and writes all information of this run
  
  open(20,file=output(1:length(output)))
  write(20,"( '#',/,&
             &'# Output of solvation (TOTAL):',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ',a,/,& 
             &'# First frame: ',i5,' Last frame: ',i5,' Stride: ',i5,/,&
             &'#',/,&
             &'# Periodic conditions: ',/,&
             &'# Periodic: ',l1,' Read from DCD: ',l1,&
             &'#',/,&
             &'# Number of atoms and mass of solute: ',i6,f12.3,/,&
             &'# First and last atoms of solute: ',i6,tr1,i6,/,&
             &'# Number of atoms and mass of solvent: ',i6,f12.3,/,&
             &'# First and last atoms of solvent: ',i6,tr1,i6,/,&
             &'#' )" )&
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, &
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  
  write(20,"( '#',/,&
             &'# Solvation of group 1 by group 2:',/,&
             &'#    TIME ',t16,'TOTAL')")

  ! Open output file for per-residue solvation if required

  if ( perresidue ) then
    open(30,file=perresidue_output(1:length(perresidue_output)))
    write(30,"( '#',/,&
               &'# Output of solvation (PER-RESIDUE):',/,&
               &'# Input file: ',a,/,& 
               &'# DCD file: ',a,/,& 
               &'# Group file: ',a,/,&
               &'# PSF file: ',a,/,& 
               &'# First frame: ',i5,' Last frame: ',i5,' Stride: ',i5,/,&
               &'#',/,&
               &'# Periodic conditions: ',/,&
               &'# Periodic: ',l1,' Read from DCD: ',l1,&
               &'#',/,&
               &'# Number of atoms and mass of solute: ',i6,f12.3,/,&
               &'# First and last atoms of solute: ',i6,tr1,i6,/,&
               &'# Number of atoms and mass of solvent: ',i6,f12.3,/,&
               &'# First and last atoms of solvent: ',i6,tr1,i6,/,&
               &'#' )" )&
               &inputfile(1:length(inputfile)),&
               &dcdfile(1:length(dcdfile)),&
               &groupfile(1:length(groupfile)),&
               &psffile(1:length(psffile)),&
               &firstframe, lastframe, stride,&
               &periodic, readfromdcd, &
               &nsolute, mass1, solute(1), solute(nsolute),& 
               &nsolvent, mass2, solvent(1), solvent(nsolvent)  
    write(30,"( '# COLUMNS CORRESPOND TO: ',/,&
               &'#       1  Time.' )" )
    do i = 1, nrsolute
      write(30,"('#',i8,tr2,a4,tr2,a4,tr2,i4 )" )&
               &i+1, segat(fatrsolute(i)),&
               &resat(fatrsolute(i)),&
               &resid(fatrsolute(i))
    end do
    write(30,"( '#',/,&
               &'# Solvation of group 1 by group 2:',/,&
               &'#    TIME ',t16,'SOLVATION RESIDUE BY RESIDUE')")
  end if
  
  ! Now going to read the dcd file
  
  memframes = memory / ntotat
  ncycles = lastframe / memframes + 1
  memlast = lastframe - memframes * ( ncycles - 1 )
  write(*,*) ' Will read and store in memory at most ', memframes,&
             ' frames per reading cycle. '
  write(*,*) ' There will be ', ncycles, ' cycles of reading. '
  write(*,*) ' Last cycle will read ', memlast,' frames. '
  write(*,*)        
  
  ! Squared shellsize for faster calculation in solvcell routine
  
  shellsize2 = shellsize**2
  
  ! Number of cells of linked cell method in each direction (initialization)
  
  do i = 1, 3
    nbdim(i) = 0
  end do
  allocate( iatomfirst(1,1,1), iatomnext(nsolvent) )
  
  ! Reading dcd file and computing the solvation
   
  iframe = 0
  do icycle = 1, ncycles 
   
    write(*,"( t3,'Cycle',t10,i5,tr2,' Reading: ',f6.2,'%')",&
          advance='no') icycle, 0.d0 
  
    ! Each cycle fills the memory as specified by the memory parameter 

    if(icycle == ncycles) then
      nfrcycle = memlast
    else
      nfrcycle = memframes
    end if

    ! Reading the coordinates
  
    iatom = 0
    do kframe = 1, nfrcycle    
      if(dcdaxis) then
        read(10) readsidesx, t, readsidesy, t, t, readsidesz
        side(kframe,1) = sngl(readsidesx)
        side(kframe,2) = sngl(readsidesy)
        side(kframe,3) = sngl(readsidesz)
      end if
      read(10) (xdcd(k), k = iatom + 1, iatom + lastatom)
      read(10) (ydcd(k), k = iatom + 1, iatom + lastatom)            
      read(10) (zdcd(k), k = iatom + 1, iatom + lastatom)           
      iatom = iatom + ntotat
      write(*,"( 7a,f6.2,'%' )",advance='no')&
           (char(8),i=1,7), 100.d0*float(kframe)/nfrcycle
    end do
  
    ! Computing the solvation
  
    write(*,"(' Computing: ',f6.2,'%')",advance='no') 0.d0
  
    iatom = 0
    do kframe = 1, nfrcycle
      iframe = iframe + 1
  
      if(mod(iframe,stride) == 0 .and. iframe >= firstframe) then
  
        ! Putting the atoms in their minimum image coordinates if using periodic
        ! boundary conditions
  
        if ( periodic ) then
          if ( readfromdcd ) then
            axis(1) = side(kframe,1)
            axis(2) = side(kframe,2)
            axis(3) = side(kframe,3)
          end if
          do i = 1, nsolute
            ii = iatom + solute(i)
            x1 = xdcd(ii)
            y1 = ydcd(ii)
            z1 = zdcd(ii)
            call image(x1,y1,z1,axis(1),axis(2),axis(3))
            xdcd(ii) = x1
            ydcd(ii) = y1
            zdcd(ii) = z1
          end do
          do i = 1, nsolvent
            ii = iatom + solvent(i)
            x1 = xdcd(ii)
            y1 = ydcd(ii)
            z1 = zdcd(ii)
            call image(x1,y1,z1,axis(1),axis(2),axis(3))
            xdcd(ii) = x1
            ydcd(ii) = y1
            zdcd(ii) = z1
          end do
        end if

        ! Prepare the linked cells
  
        if ( periodic ) then
          xmin(1) = -axis(1) / 2.
          xmin(2) = -axis(2) / 2.
          xmin(3) = -axis(3) / 2.
          xmax(1) = axis(1) / 2.
          xmax(2) = axis(2) / 2.
          xmax(3) = axis(3) / 2.
        else
          do i = 1, 3
            xmin(i) = 1.e30
            xmax(i) = -1.e30
          end do
          do i = 1, nsolute
            ii = iatom + solute(i)
            xmin(1) = min(xmin(1),xdcd(ii))
            xmin(2) = min(xmin(2),ydcd(ii))
            xmin(3) = min(xmin(3),zdcd(ii))
            xmax(1) = max(xmax(1),xdcd(ii))
            xmax(2) = max(xmax(2),ydcd(ii))
            xmax(3) = max(xmax(3),zdcd(ii))
          end do
          do i = 1, nsolvent
            ii = iatom + solvent(i)
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
        nboxes(1) = max(1,int(boxlength/shellsize))
        dbox_x = boxlength / nboxes(1)

        boxlength = xmax(2) - xmin(2) 
        nboxes(2) = max(1,int(boxlength/shellsize))
        dbox_y = boxlength / nboxes(2)

        boxlength = xmax(3) - xmin(3) 
        nboxes(3) = max(1,int(boxlength/shellsize))
        dbox_z = boxlength / nboxes(3)

        ! If the number of cells changed from the previous frame, update dimensions

        if ( nboxes(1) > nbdim(1) .or. &
             nboxes(2) > nbdim(2) .or. &
             nboxes(3) > nbdim(3) ) then
          deallocate( iatomfirst )
          nbdim(1) = nboxes(1)
          nbdim(2) = nboxes(2)
          nbdim(3) = nboxes(3)
          allocate( iatomfirst(0:nbdim(1)+1,0:nbdim(2)+1,0:nbdim(3)+1) )
        end if

        ! Reseting the iatomfirst array
  
        do i = 0, nboxes(1) + 1
          do j = 0, nboxes(2) + 1
            do k = 0, nboxes(3) + 1
              iatomfirst(i,j,k) = 0
            end do
          end do
        end do

        ! Putting the atoms in their boxes
  
        do i = 1, nsolute
          ii = iatom + solute(i)
          ibox = int( (xdcd(ii) - xmin(1)) / dbox_x ) + 1
          jbox = int( (ydcd(ii) - xmin(2)) / dbox_y ) + 1
          kbox = int( (zdcd(ii) - xmin(3)) / dbox_z ) + 1
          if ( ibox == nbdim(1)+1 ) ibox = nbdim(1)
          if ( jbox == nbdim(2)+1 ) jbox = nbdim(2)
          if ( kbox == nbdim(3)+1 ) kbox = nbdim(3)
          solute_box(i,1) = ibox
          solute_box(i,2) = jbox
          solute_box(i,3) = kbox
        end do
        do i = 1, nsolvent
          ii = iatom + solvent(i)
          ibox = int( (xdcd(ii) - xmin(1)) / dbox_x ) + 1
          jbox = int( (ydcd(ii) - xmin(2)) / dbox_y ) + 1
          kbox = int( (zdcd(ii) - xmin(3)) / dbox_z ) + 1
          if ( ibox == nbdim(1)+1 ) ibox = nbdim(1)
          if ( jbox == nbdim(2)+1 ) jbox = nbdim(2)
          if ( kbox == nbdim(3)+1 ) kbox = nbdim(3)
          iatomnext(i) = iatomfirst(ibox,jbox,kbox)
          iatomfirst(ibox,jbox,kbox) = i
        end do

        ! Filling up boundaries of periodic cell with phantom copies of the solvent atoms

        if ( periodic ) call phantomcells(nboxes,iatomfirst,nbdim)

        ! Loop over solute atoms

        nsolv = 0
        do isolute = 1, nsolute
  
          i = solute_box(isolute,1)  
          j = solute_box(isolute,2)  
          k = solute_box(isolute,3)  
  
          ! Check on the adjacent boxes if there is an atom of the
          ! solvent which is close enough
  
          ! Inside box
  
          call solvcell(iatom,isolute,i,j,k,solv_list,nsolv)
  
          ! Interactions of boxes that share faces
  
          call solvcell(iatom,isolute,i+1,j,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j+1,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j,k+1,solv_list,nsolv)
  
          call solvcell(iatom,isolute,i-1,j,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j-1,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j,k-1,solv_list,nsolv)
  
          ! Interactions of boxes that share axes
  
          call solvcell(iatom,isolute,i+1,j+1,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i+1,j,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i+1,j-1,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i+1,j,k-1,solv_list,nsolv)
  
          call solvcell(iatom,isolute,i,j+1,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j+1,k-1,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j-1,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i,j-1,k-1,solv_list,nsolv)
  
          call solvcell(iatom,isolute,i-1,j+1,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i-1,j,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i-1,j-1,k,solv_list,nsolv)
          call solvcell(iatom,isolute,i-1,j,k-1,solv_list,nsolv)
  
          ! Interactions of boxes that share vertices
  
          call solvcell(iatom,isolute,i+1,j+1,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i+1,j+1,k-1,solv_list,nsolv)
          call solvcell(iatom,isolute,i+1,j-1,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i+1,j-1,k-1,solv_list,nsolv)
  
          call solvcell(iatom,isolute,i-1,j+1,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i-1,j+1,k-1,solv_list,nsolv)
          call solvcell(iatom,isolute,i-1,j-1,k+1,solv_list,nsolv)
          call solvcell(iatom,isolute,i-1,j-1,k-1,solv_list,nsolv)

        end do

        ! Now analysing the solvation list

        ! Computing the total number of residues of the solvent that solvate the solute

        do i = 1, nres
          rsolv(i) = .false.
        end do
        do i = 1, nsolv
          rsolv(solv_list(i,2)) = .true.
        end do
        nrsolv = count(rsolv)
        write(20,"( e12.6,10000(tr2,i8) )") scaletime*iframe, nrsolv

        ! Now computing the solvation per residue of the solute (this is not
        ! done efficiently here...)

        if ( perresidue ) then
          do i = 1, nrsolute 
            do j = 1, nres
              rsolv(j) = .false.
            end do
            do j = 1, nsolv
              if ( solv_list(j,1) == rseqsolute(i) ) then
                rsolv(solv_list(j,2)) = .true.
              end if
            end do
            solv_per_residue(i) = count(rsolv)
          end do
          write(30,"( e12.6,10000(tr2,i8) )") scaletime*iframe, (solv_per_residue(i),i=1,nrsolute)
        end if

      end if
  
      ! Printing progress
  
      write(*,"( 7a,f6.2,'%' )",advance='no')&
           (char(8),i=1,7), 100.d0*float(kframe)/nfrcycle
  
      iatom = iatom + ntotat
    end do
    write(*,*)
           
  end do
  close(10)
  close(20)
  close(30)
  
  ! Write final messages with names of output files and their content
  
  time0 = etime(tarray) - time0
  write(*,*)
  write(*,"( tr2,52('-') )")
  write(*,*)
  write(*,*) ' OUTPUT FILES: ' 
  write(*,*)
  write(*,*) ' Wrote output file: ', output(1:length(output))
  write(*,*)
  write(*,*) ' Which contains the total and residue-by-residue'
  write(*,*) ' solvation of the solute by the solvent.'
  write(*,*) 
  write(*,*) ' Running time: ', time0
  write(*,*) '####################################################'
  write(*,*) 
  write(*,*) '  END: Normal termination.  '
  write(*,*) 
  write(*,*) '####################################################'
  write(*,*)        

end program solvation

!
! Subroutine that checks the distance crterium and annotates the
! possible satisfaction of the solvation criterium
!

subroutine solvcell(iatom,isolute,ibox,jbox,kbox,solv_list,nsolv)

   use linkedcells_solvation, only : xdcd, ydcd, zdcd, iatomnext, shellsize2, &
                                     solute, solvent, resid, nsolvlist, iatomfirst, &
                                     axis, periodic, nboxes
   implicit none
   real :: d2, x1, y1, z1
   integer :: ii, jj, nsolv, iatom, isolute, isolvent, solv_list(nsolvlist,2), &
              ibox, jbox, kbox

   ii = iatom + solute(isolute)
   isolvent = iatomfirst(ibox,jbox,kbox)
   do while( isolvent /= 0 )

     jj = iatom + solvent(isolvent)

     x1 = xdcd(jj)
     y1 = ydcd(jj)
     z1 = zdcd(jj)
     if ( periodic ) call movephantomcoor(xdcd(jj),ydcd(jj),zdcd(jj),x1,y1,z1,&
                                          nboxes,ibox,jbox,kbox,axis)

     d2 = ( xdcd(ii) - x1 )**2 + &
          ( ydcd(ii) - y1 )**2 + &
          ( zdcd(ii) - z1 )**2

     if ( d2 <= shellsize2 ) then
       nsolv = nsolv + 1
       if ( nsolv > nsolvlist ) then
         write(*,*)
         write(*,*) ' ERROR: Increase the nsolvlist parameter in the input file. Try ', nsolv*2 
         stop
       end if
       solv_list(nsolv,1) = resid(solute(isolute))
       solv_list(nsolv,2) = resid(solvent(isolvent))
     end if

     isolvent = iatomnext(isolvent)
   end do

end subroutine solvcell
