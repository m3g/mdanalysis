!
! gss: A program to compute gss[1] radial distribution functions from
!      molecular dynamics simulations in NAMD DCD format.
!
!      Important: THIS IS NOT THE CLASSICAL RADIAL DISTRIBUTION
!                 FUNCTION. It is the shape-dependent RDF used
!                 for non-spherical species. It will only coincide
!                 with the classical RDF for perfectly spherical
!                 species. The normalization of this distribution 
!                 function is more complicated than the normalization
!                 of the radial distribution function for spherical
!                 solutes. Here, we estimate the volume of the 
!                 solute* by partitioning the space into bins
!                 of side "probeside", and count how many bins contain
!                 atoms of the solute. With this estimate of the
!                 volume of the solute, the number of solvent
!                 that would be present in the simulation box for
!                 the same solvent density, but without the solute,
!                 is calculated. This number of random solvent molecules
!                 is then generated, and the gss relative to the
!                 solute for these molecules is computed for 
!                 normalization. This procedure is performed 
!                 independently for each solute structure in each
!                 frame of the trajectory.
!                 *The volume can be of a different selection than
!                  the solute selection.
!
! Please cite the following reference when using this package:
!
! I. P. de Oliveira, L. Martinez, Molecular basis of co-solvent induced
! stability of lipases. To be published. 
!
! Reference of this distribution function:
! 1. W. Song, R. Biswas and M. Maroncelli. Intermolecular interactions
! and local density augmentation in supercritical solvation: A survey of
! simulation and experimental results. Journal of Physical Chemistry A,
! 104:6924-6939, 2000.  
!
! Auxiliar dimensions:
!          memory: Controls the amount of memory used for reading dcd
!                  files. If the program ends with segmentation fault
!                  without any other printing, decrease the value of
!                  this parameter.  
!
! L. Martinez, Mar 13, 2014.
! Institute of Chemistry, State University of Campinas (UNICAMP)
!

program g_solute_solvent
 
  ! Static variables

  implicit none
  integer, parameter :: memory=15000000
  integer :: natom, nsolute, nsolvent, nexclude, isolute, isolvent, nsolvent2, &
             iexclude, narg, length, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
             j, iframe, icycle, nfrcycle, iatom, k, ii, jj, &
             status, keystatus, ndim, iargc, lastatom, nres, nrsolute, nrexclude,&
             nrsolvent, kframe, irad, gsssum, nslabs, natoms_solvent,&
             gsssum_random, nsmalld, nrsolvent2, ibox, jbox, kbox,&
             noccupied, nbdim(3), nboxes(3), maxsmalld, frames
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: side(memory,3), mass1, mass2, seed,&
          cmx, cmy, cmz, beta, gamma, theta, random, axis(3)
  real, parameter :: twopi = 2.*3.1415925655
  real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
          x1, y1, z1, time0, etime, tarray(2),&
          gssnorm, gssstep, &
          density, dbox_x, dbox_y, dbox_z, cutoff, probeside, exclude_volume,&
          totalvolume, xmin(3), xmax(3), gssmax, kbint, gsslast, radius, bulkdensity
  character(len=200) :: groupfile, line, record, value, keyword,&
                        dcdfile, inputfile, psffile, file,&
                        output
  character(len=4) :: dummyc
  logical :: readfromdcd, dcdaxis, periodic, scalelast
  real :: shellvolume
  
  ! Allocatable arrays
  
  integer, allocatable :: solute(:), solvent(:), resid(:), solute2(:), solvent2(:), &
                          natres(:), fatres(:), fatrsolute(:), fatrsolvent(:), exclude(:),&
                          fatrexclude(:),&
                          nrsolv(:), gss(:), gss_random(:), gss_max(:), gss_max_random(:),&
                          gss_avg(:), gss_avg_random(:), nrsolv2(:),&
                          ismalld(:)
  real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:),&
                       solvent_molecule(:,:), &
                       xref(:), yref(:), zref(:), xrnd(:), yrnd(:), zrnd(:),&
                       x(:), y(:), z(:), mind(:), dsmalld(:)
  character(len=6), allocatable :: class(:), segat(:), resat(:),&
                                   typeat(:), classat(:)
  logical, allocatable :: hasatoms(:,:,:)
  
  ! Compute time
  
  time0 = etime(tarray)
  
  ! Output title
  
  write(*,"(/,' ####################################################',&
            &/,/,&
            & '   GSS: Compute gss distribution from DCD files      ',&
            &/,/,&
            & ' ####################################################',/)")    
  
  call version()
  
  ! Seed for random number generator
  
  seed = 0.48154278727e0
  
  ! Some default parameters
  
  firstframe = 1
  lastframe = 0
  stride = 1
  periodic = .true.
  readfromdcd = .true.
  nslabs = 1000 
  cutoff = 10.
  gssmax = 20.
  density = 1.
  probeside = 2.
  scalelast = .true.
  
  ! Open input file and read parameters
  
  narg = iargc()
  if(narg == 0) then
    write(*,*) ' Run with: ./gss input.inp '
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
    else if(keyword(record) == 'stride') then
      line = value(record)
      read(line,*,iostat=keystatus) stride
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'nslabs') then
      line = value(record)
      read(line,*,iostat=keystatus) nslabs
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'cutoff') then
      line = value(record)
      read(line,*,iostat=keystatus) cutoff
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'gssmax') then
      line = value(record)
      read(line,*,iostat=keystatus) gssmax
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'probeside') then
      line = value(record)
      read(line,*,iostat=keystatus) probeside
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'density') then
      line = value(record)
      read(line,*,iostat=keystatus) density
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
      write(*,*) ' GSS output file name: ', output(1:length(output))
    else if(keyword(record) == 'solute' .or. &
            keyword(record) == 'exclude' .or. &
            keyword(record) == 'solvent') then
      write(*,"(a,/,a)") ' ERROR: The options solute and solvent must be used ',&
                         '        with the gss.sh script, not directly. '
      stop
    else if(keyword(record) == 'scalelast') then
      line = value(record)
      if ( trim(line) == 'yes' ) scalelast = .true.
      if ( trim(line) == 'no' ) scalelast = .false.
    else if(record(1:1) /= '#' .and. & 
            keyword(record) /= 'par' .and. &
            record(1:1) > ' ') then
      write(*,*) ' ERROR: Unrecognized keyword found: ',keyword(record)
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
            resid(ndim), natres(ndim), fatres(ndim), fatrsolute(ndim),&
            fatrsolvent(ndim), nrsolv(ndim), fatrexclude(ndim),&
            nrsolv2(ndim*2), x(ndim*2), y(ndim*2), z(ndim*2), solvent2(ndim*2),&
            solute2(ndim), mind(ndim*2) )
  
  ! Reading parameter files to get the vdW sigmas for the definition of exclusion
  ! zones
  
  nclass = 0
  open(99,file=inputfile,action='read',status='old')
  do while(.true.)
    read(99,"( a200 )",iostat=status) record
    if(status /= 0) exit
    if(keyword(record) == 'par') then
      file = value(record)
      write(*,*) ' Reading parameter file: ', file(1:length(file))
      call readpar(file,nclass,class,eps,sig)
    end if
  end do
  close(99)
  
  ! Allocate gss array according to nslabs

  if ( gssmax > cutoff ) then
    write(*,*) ' Warning: gssmax > cutoff, the actual gssmax will be equal to cutoff. '
    nslabs = int(cutoff/(gssmax/nslabs))
    write(*,*) '          to keep same bin precision, nlsabs = ', nslabs 
    gssmax = cutoff
  end if
  
  write(*,*) ' Number of volume slabs: ', nslabs
  allocate( gss(nslabs), gss_random(nslabs), gss_max(nslabs), gss_max_random(nslabs),&
            gss_avg(nslabs), gss_avg_random(nslabs) )

  gssstep = gssmax / float(nslabs)
         
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
  write(*,*) ' Cutoff for linked cells: ', cutoff
  
  ! Read PSF file
  
  write(*,*) ' Reading PSF file: ', psffile(1:length(psffile))
  call readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
               resid,classat,typeat,q,e,s,mass,.false.)
  nres = resid(natom)
  write(*,*) ' Number of atoms in PSF file: ', natom
  write(*,*) ' Number of residues in PSF file: ', nres
  write(*,*)
  
  ! Read solute and solvent and exclude information from file
  ! First reading the size of the groups to allocate arrays

  open(10,file=groupfile,action='read')
  nsolute = 0
  nsolvent = 0
  nexclude = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      if(line(63:66) == '1.00') then      
        nsolute = nsolute + 1        
      end if
      if(line(63:66) == '2.00') then      
        nsolvent = nsolvent + 1        
      end if
      if(line(57:60) == '1.00') then
        nexclude = nexclude + 1        
      end if
    end if     
  end do
  if ( nsolute < 1 .or. nsolvent < 1 .or. nexclude < 1 ) then
    write(*,*) ' ERROR: No atom selected for solute, solvent, or exclude. '
    write(*,*) '        nsolute = ', nsolute
    write(*,*) '        nsolvent = ', nsolvent
    write(*,*) '        nexclude = ', nexclude
    stop
  end if
  allocate ( solute(nsolute), solvent(nsolvent), exclude(nexclude) )
  close(10)
  
  ! Now reading reading the group atoms
  
  open(10,file=groupfile,action='read')
  isolute = 0
  isolvent = 0
  iexclude = 0
  mass1 = 0.
  mass2 = 0.
  iatom = 0
  do 
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) exit
    if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
      iatom = iatom + 1

      ! Read atoms belonging to solute

      if(line(63:66) == '1.00') then      
        isolute = isolute + 1        
        solute(isolute) = iatom
        mass1 = mass1 + mass(iatom)
      end if
  
      ! Read atoms belonging to solvent

      if(line(63:66) == '2.00') then
        isolvent = isolvent + 1        
        solvent(isolvent) = iatom
        mass2 = mass2 + mass(iatom)
      end if

      ! Read atoms belonging to solvent

      if(line(57:60) == '1.00') then
        iexclude = iexclude + 1        
        exclude(iexclude) = iatom
      end if

    end if     
  end do
  close(10)
  lastatom = max0(solute(nsolute),solvent(nsolvent),exclude(nexclude))

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
    nrsolv(i) = nrsolvent
  end do
  j = 0
  nrexclude = 0
  do i = 1, nexclude
    if(resid(exclude(i)).gt.j) then
      nrexclude = nrexclude + 1 
      fatrexclude(nrexclude) = fatres(resid(exclude(i)))
      j = resid(exclude(i))
    end if
  end do

  ! This is for the initialization of the smalldistances routine

  allocate( ismalld(1), dsmalld(1) )

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
  write(*,*) 
  write(*,*) ' Excluded volume information: ' 
  write(*,*) ' Probe side for solute volume estimation: ', probeside
  write(*,*) ' Number of atoms of excluded volume: ', nexclude
  write(*,*) ' First atom of excluded volume: ', exclude(1)
  write(*,*) ' Last atom of excluded volume: ', exclude(nexclude)
  write(*,*) ' Number of residues in exclude: ', nrexclude

  ! Check if the solvent atoms have obvious reading problems
  
  if ( mod(nsolvent,nrsolvent) /= 0 ) then
    write(*,*) ' ERROR: Incorrect count of solvent atoms or residues. '
    stop
  end if

  natoms_solvent = nsolvent / nrsolvent 
  write(*,*)  ' Number of atoms of each solvent molecule: ', natoms_solvent
  
  ! Allocate solvent molecule (this will be used to generate random coordinates
  ! for each solvent molecule, one at a time, later)
  
  allocate( solvent_molecule(natoms_solvent,3), &
            xref(natoms_solvent), yref(natoms_solvent), zref(natoms_solvent),&
            xrnd(natoms_solvent), yrnd(natoms_solvent), zrnd(natoms_solvent) )
  
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
  
  ! Number of frames (used for output only)
  
  frames=(lastframe-firstframe+1)/stride
  write(*,*) ' Number of frames to read: ', frames

  ! Now going to read the dcd file
  
  memframes = memory / ntotat
  ncycles = lastframe / memframes + 1
  memlast = lastframe - memframes * ( ncycles - 1 )
  write(*,*) ' Will read and store in memory at most ', memframes,&
             ' frames per reading cycle. '
  write(*,*) ' There will be ', ncycles, ' cycles of reading. '
  write(*,*) ' Last cycle will read ', memlast,' frames. '
  write(*,*)        
  
  ! Reseting the gss distribution function
  
  do i = 1, nslabs
    gss(i) = 0
    gss_random(i) = 0
    gss_max(i) = 0
    gss_max_random(i) = 0
    gss_avg(i) = 0
    gss_avg_random(i) = 0
  end do

  ! Initializing hasatoms array

  nbdim(1) = 1
  nbdim(2) = 1
  nbdim(3) = 1
  allocate( hasatoms(0:nbdim(1)+1,0:nbdim(2)+1,0:nbdim(3)+1))

  ! Reading dcd file and computing the gss function
   
  iframe = 0
  do icycle = 1, ncycles 
   
    write(*,"( t3,'Cycle',t10,i5,tr2,' Reading: ',f6.2,'%')",&
          advance='no') icycle, 0. 
  
    ! Each cycle fills the memory as specified by the memory parameter 
  
    if(icycle == ncycles) then
      nfrcycle = memlast
    else
      nfrcycle = memframes
    end if
  
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
           (char(8),i=1,7), 100.*float(kframe)/nfrcycle
    end do
    write(*,"(' Computing: ',f6.2,'%')",advance='no') 0.
  
    ! Computing the gss function
  
    iatom = 0
    do kframe = 1, nfrcycle
      iframe = iframe + 1

      if(mod(iframe-firstframe,stride) /= 0 .or. iframe < firstframe ) then
        iatom = iatom + ntotat
        cycle
      end if

      ! Sides of the periodic cell in this frame
  
      axis(1) = side(kframe,1) 
      axis(2) = side(kframe,2) 
      axis(3) = side(kframe,3) 

      ! Normalization: creates a box with random coordinates for the solvent molecules, excluding
      ! the excluded selection volume. The gss is computed for this random box independently, and this is 
      ! the normalization for the gss that allows for the computation of the potential of mean
      ! force.

      ! 
      ! Estimate the volume of the solute in this frame
      !

      call getmaxmin(nexclude,exclude,iatom,xdcd,ydcd,zdcd,xmin,xmax,.true.)
      nboxes(1) = max(1,int((xmax(1)-xmin(1))/probeside))
      dbox_x = max((xmax(1)-xmin(1)) / nboxes(1),probeside)
      nboxes(2) = max(1,int((xmax(2)-xmin(2))/probeside))
      dbox_y = max((xmax(2)-xmin(2)) / nboxes(2),probeside)
      nboxes(3) = max(1,int((xmax(3)-xmin(3))/probeside))
      dbox_z = max((xmax(3)-xmin(3)) / nboxes(3),probeside)
      if ( nboxes(1) > nbdim(1) .or. &
           nboxes(2) > nbdim(2) .or. &
           nboxes(3) > nbdim(3) ) then
         deallocate( hasatoms )
         nbdim(1) = nboxes(1)
         nbdim(2) = nboxes(2)
         nbdim(3) = nboxes(3)
         allocate( hasatoms(0:nbdim(1)+1,0:nbdim(2)+1,0:nbdim(3)+1) )
      end if
      do i = 0, nboxes(1)+1
        do j = 0, nboxes(2)+1
          do k = 0, nboxes(3)+1
            hasatoms(i,j,k) = .false.
          end do
        end do
      end do
      do i = 1, nexclude
        ii = iatom + exclude(i)
        x1 = xdcd(ii)
        y1 = ydcd(ii)
        z1 = zdcd(ii)
        ibox = int( (x1 - xmin(1)) / dbox_x ) + 1
        jbox = int( (y1 - xmin(2)) / dbox_y ) + 1
        kbox = int( (z1 - xmin(3)) / dbox_z ) + 1
        if ( ibox == nbdim(1)+1 ) ibox = nbdim(1)
        if ( jbox == nbdim(2)+1 ) jbox = nbdim(2)
        if ( kbox == nbdim(3)+1 ) kbox = nbdim(3)
        hasatoms(ibox,jbox,kbox) = .true.
      end do
      noccupied = 0
      do i = 0, nboxes(1)+1
        do j = 0, nboxes(2)+1 
          do k = 0, nboxes(3)+1 
            if( hasatoms(i,j,k) ) noccupied = noccupied + 1
          end do
        end do
      end do
      exclude_volume = noccupied * ( probeside**3 ) 

      ! Computing the number of random solvent molecules that has to be
      ! generated

      totalvolume = axis(1)*axis(2)*axis(3)
      nrsolvent2 = nrsolvent * int(totalvolume / ( totalvolume - exclude_volume ))
      nsolvent2 = nrsolvent2*natoms_solvent

      ! Create random coordinates for 'nrsolvent2' solvent molecules, inside
      ! the box

      do isolvent = 1, nrsolvent2

        ! First, pick randomly a solvent molecule from the box, to minimize conformational
        ! biases of the solvent molecules
    
        ii = (nrsolvent-1)*int(random(seed)) + 1
    
        ! Save the coordinates of this molecule in this frame in the solvent_molecule array
    
        jj = iatom + solvent(1) + natoms_solvent*(ii-1)
        do i = 1, natoms_solvent
          solvent_molecule(i,1) = xdcd(jj+i-1)
          solvent_molecule(i,2) = ydcd(jj+i-1)
          solvent_molecule(i,3) = zdcd(jj+i-1)
        end do
  
        ! Put molecule in its center of mass for it to be the reference coordinate for the 
        ! random coordinates that will be generated
  
        cmx = 0.
        cmy = 0.
        cmz = 0.
        do i = 1, natoms_solvent
          cmx = cmx + solvent_molecule(i,1)
          cmy = cmy + solvent_molecule(i,2)
          cmz = cmz + solvent_molecule(i,3)
        end do
        cmx = cmx / float(natoms_solvent)
        cmy = cmy / float(natoms_solvent)
        cmz = cmz / float(natoms_solvent)
        do i = 1, natoms_solvent
          xref(i) = solvent_molecule(i,1) - cmx
          yref(i) = solvent_molecule(i,2) - cmy
          zref(i) = solvent_molecule(i,3) - cmz
        end do
  
        ! Generate a random position for this molecule
  
        cmx = -axis(1)/2. + random(seed)*axis(1) 
        cmy = -axis(2)/2. + random(seed)*axis(2) 
        cmz = -axis(3)/2. + random(seed)*axis(3) 
        beta = random(seed)*twopi
        gamma = random(seed)*twopi
        theta = random(seed)*twopi
        call compcart(natoms_solvent,xref,yref,zref,xrnd,yrnd,zrnd,&
                      cmx,cmy,cmz,beta,gamma,theta)

        ! Add this molecule to x, y, z arrays

        do i = 1, natoms_solvent
          call image(xrnd(i),yrnd(i),zrnd(i),side(kframe,1),side(kframe,2),side(kframe,3))
          ii = nsolute + (isolvent-1)*natoms_solvent + i
          x(ii) = xrnd(i)
          y(ii) = yrnd(i)
          z(ii) = zrnd(i)
          solvent2(ii-nsolute) = ii
          ! Annotate to which molecule this atom pertains
          nrsolv2(ii-nsolute) = isolvent
        end do
         
      end do

      ! Add solute coordinates to x, y and z arrays

      do isolute = 1, nsolute
        ii = iatom + solute(isolute)
        x1 = xdcd(ii)
        y1 = ydcd(ii)
        z1 = zdcd(ii)
        if ( periodic ) call image(x1,y1,z1,axis(1),axis(2),axis(3)) 
        x(isolute) = x1
        y(isolute) = y1
        z(isolute) = z1
        solute2(isolute) = isolute
      end do

      ! Now compute all the distances that are smaller than the desired cutoff

      maxsmalld = nsolute*nsolvent
      if ( maxsmalld > size(dsmalld) ) then
        deallocate( ismalld, dsmalld )
        allocate( ismalld(maxsmalld), dsmalld(maxsmalld) )
      end if
      call smalldistances(nsolute,solute2,nsolvent2,solvent2,x,y,z,cutoff,&
                          nsmalld,ismalld,dsmalld,axis,maxsmalld)
      
      !
      ! Computing the gss functions from distance data
      !

      ! For each solvent residue, get the MINIMUM, MAXIMUM and AVERAGE distances to the solute
    
      do i = 1, nrsolvent2
        mind(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        isolvent = nrsolv2(ismalld(i))
        if ( dsmalld(i) < mind(isolvent) ) then
          mind(isolvent) = dsmalld(i)
        end if
      end do

      ! Summing up current data to the gss histogram

      do i = 1, nrsolvent2
        irad = int(float(nslabs)*mind(i)/gssmax)+1
        if ( irad <= nslabs ) then
          gss_random(irad) = gss_random(irad) + 1
        end if
      end do

      !
      ! Now computing the GSS data from the simulation itself
      !

      do i = 1, nsolute
        ii = iatom + solute(i)
        x(solute(i)) = xdcd(ii)
        y(solute(i)) = ydcd(ii)
        z(solute(i)) = zdcd(ii)
      end do
      do i = 1, nsolvent
        ii = iatom + solvent(i)
        x(solvent(i)) = xdcd(ii)
        y(solvent(i)) = ydcd(ii)
        z(solvent(i)) = zdcd(ii)
      end do

      ! Compute all distances that are smaller than the cutoff

      maxsmalld = nsolute*nsolvent
      if ( maxsmalld > size(dsmalld) ) then
        deallocate( ismalld, dsmalld )
        allocate( ismalld(maxsmalld), dsmalld(maxsmalld) )
      end if
      call smalldistances(nsolute,solute,nsolvent,solvent,x,y,z,cutoff,&
                          nsmalld,ismalld,dsmalld,axis,maxsmalld)

      !
      ! Computing the gss functions from distance data
      !

      ! For each solvent residue, get the MINIMUM distance to the solute
    
      do i = 1, nrsolvent
        mind(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        isolvent = nrsolv(ismalld(i))
        if ( dsmalld(i) < mind(isolvent) ) then
          mind(isolvent) = dsmalld(i)
        end if
      end do

      ! Summing up current data to the gss histogram

      do i = 1, nrsolvent
        irad = int(float(nslabs)*mind(i)/gssmax)+1
        if( irad <= nslabs ) then
          gss(irad) = gss(irad) + 1
        end if
      end do

      ! Print progress

      write(*,"( 7a,f6.2,'%' )",advance='no')&
           (char(8),i=1,7), 100.*float(kframe)/nfrcycle
  
      iatom = iatom + ntotat
    end do
    write(*,*)
  end do
  close(10)
  
  ! Open output file and writes all information of this run

  !
  ! GSS computed with minimum distance
  !
  
  open(20,file=output(1:length(output)))
  write(20,"( '#',/,&
             &'# Output of gss.f90: Using MINIMUM distance to solute.',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ',a,/,& 
             &'# First frame: ',i5,' Last frame: ',i5,' Stride: ',i5,/,&
             &'#',/,&
             &'# Periodic boundary conditions: ',/,&
             &'# Periodic: ',l1,' Read from DCD: ',l1,/,&
             &'#',/,&
             &'# Number of atoms and mass of group 1: ',i6,f12.3,/,&
             &'# First and last atoms of group 1: ',i6,tr1,i6,/,&
             &'# Number of atoms and mass of group 2: ',i6,f12.3,/,&
             &'# First and last atoms of group 2: ',i6,tr1,i6,/,&
             &'#' )" )&
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, &
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  

  ! Compute error of gssrand distribution at large distance (expected to be 1.00)
  
  if ( gss(nslabs) < 1.e-10 ) then
    write(*,*) ' ERROR: Something wrong with random normalization. Contact the developer. '
    stop
  end if
  gsslast = float(gss_random(nslabs))/float(gss(nslabs))
  write(20,"('# Error in random normalization at large distances: ',f12.3,'%' )") (gsslast - 1.d0)*100
  if ( scalelast ) then
    write(20,"('# scalelast is true, so GSS/GSSRND (column 2) is divided by: ',f12.3 )") gsslast
  end if

  radius = nslabs*gssstep - gssstep/2.
  bulkdensity = gss(nslabs) / (frames*shellvolume(radius,gssstep))
  write(20,"('# Solvent density at bulk (largest distance) in sites/vol: ',e12.5 )") bulkdensity
  write(20,"('#')")

  ! Output table

  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
             &'#       1  Minimum distance to solute (dmin)',/,&
             &'#       2  GSS normalized by the GSS RAND distribution. ',/,&
             &'#       3  GSS normalized according to spherical volume of radius dmin.',/,&
             &'#       4  GSS not normalized at all (just site count for each dmin, averaged over frames)',/,&
             &'#       5  Cumulative sum of sites (averaged over the number of frames) ',/,&
             &'#       6  GSS computed from random solvent distribution, not normalized ',/,&
             &'#       7  Cumulative sum of sites for the random distribution, averaged on frames.',/,&
             &'#       8  Kirwood-Buff integral (int gss - 1) ')")
  write(20,"( '#',/,&      
   &'#',t5,'1-DISTANCE',t17,'2-GSS/GSSRND',t32,'3-GSS/SPHER',t52,'4-GSS',t64,'5-CUMUL',&
   &t76,'6-GSS RND',t88,'7-CUMUL RND',t105,'8-KB INT' )" )

  gsssum = 0
  gsssum_random = 0
  kbint = 0.e0
  do i = 1, nslabs
    if ( gss(i) == 0 ) cycle
    radius = i*gssstep - gssstep/2.
    gsssum = gsssum + gss(i)
    gsssum_random = gsssum_random + gss_random(i)
    ! Normalization by spherical shell of this radius
    gssnorm = gss(i) / ( bulkdensity*shellvolume(radius,gssstep) )
    ! Normalization by random distribution of molecules
    x1 = float(gss(i))
    y1 = float(gss_random(i))/gsslast
    if ( y1 > 0. ) then
      z1 = x1 / y1
    else
      z1 = 0.
    end if
    kbint = kbint + z1 - 1.e0
    write(20,"( 8(tr2,f12.7) )")&
    radius, z1, gssnorm/frames, float(gss(i))/frames, float(gsssum)/frames, &
            float(gss_random(i))/frames, float(gsssum_random)/frames, kbint
  end do
  close(20)

  ! Write final messages with names of output files and their content
  
  time0 = etime(tarray) - time0
  write(*,*)
  write(*,"( tr2,52('-') )")
  write(*,*)
  write(*,*) ' OUTPUT FILES: ' 
  write(*,*)
  write(*,*) ' Wrote GSS output file: ', output(1:length(output))
  write(*,*)
  write(*,*) ' Which contains the volume-normalized and'
  write(*,*) ' unnormalized gss functions. '
  write(*,*) 
  write(*,*) ' Running time: ', time0
  write(*,*) '####################################################'
  write(*,*) 
  write(*,*) '  END: Normal termination.  '
  write(*,*) 
  write(*,*) '####################################################'
  write(*,*)        

end program g_solute_solvent

! Computes the volume of the spherical shell of radius
! 'radius', defined within [radius-step/2,radius+step,2].

real function shellvolume(radius,step)

  implicit none
  real :: radius, step
  real, parameter :: fourpi = 4.*3.1415925655

  shellvolume = fourpi * 2 * radius * step

end function shellvolume










