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
  integer :: maxatom, nrandom
  integer :: natom, nsolute, nsolvent, isolute, isolvent, &
             narg, length, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
             j, iframe, icycle, nfrcycle, iatom, k, ii, &
             status, keystatus, iargc, lastatom, nres, nrsolute,&
             nrsolvent, kframe, irad, nbins, natoms_solvent,&
             nsmalld, & 
             maxsmalld, frames
  integer :: nrsolvent_random
  real :: dbulk
  integer :: ibulk, nintegral
  real :: site_sum, convert
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: side(memory,3), mass1, mass2, seed, random, axis(3)
  real, parameter :: twopi = 2.*3.1415925655
  real, parameter :: mole = 6.022140857e23
  real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
          time0, etime, tarray(2),&
          binstep, &
          cutoff,  kbint, kbintsphere
  real :: bulkdensity, totalvolume, bulkvolume, simdensity, solutevolume
  real :: bulkerror, sdbulkerror
  character(len=200) :: groupfile, line, record, value, keyword,&
                        dcdfile, inputfile, psffile, file,&
                        output
  character(len=4) :: dummyc
  logical :: readfromdcd, dcdaxis, periodic 
  real :: shellradius, rshift
  real :: sphericalshellvolume, sphereradiusfromshellvolume

  integer :: jj
  real :: beta, gamma, theta, cmx, cmy, cmz

  ! Allocatable arrays
  
  integer, allocatable :: irandom(:), solute2(:), solvent_random(:)
  integer, allocatable :: solute(:), solvent(:), resid(:), &
                          natres(:), fatres(:), fatrsolute(:), fatrsolvent(:), &
                          irsolv(:), ismalld(:), irsolv_random(:)

  real, allocatable :: site_count(:), site_count_random(:), random_count(:), shellvolume(:)
  real, allocatable :: gss(:)

  real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:),&
                       solvent_molecule(:,:), &
                       xref(:), yref(:), zref(:), xrnd(:), yrnd(:), zrnd(:),&
                       x(:), y(:), z(:), mind(:), dsmalld(:)
  character(len=6), allocatable :: class(:), segat(:), resat(:),&
                                   typeat(:), classat(:)
  
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
  nbins = 1000
  nintegral = 5
  dbulk = 12.
  cutoff = 14.
  binstep = 0.1e0
  
  ! Open input file and read parameters
  
  narg = iargc()
  if(narg == 0) then
    write(*,*) ' Run with: ./gss input.inp '
    stop
  end if   
  call getarg(1,record)
  
  inputfile = record(1:length(record))
  open(10,file=inputfile,action='read',iostat=status,status='old')
  if ( status /= 0 ) then
    write(*,*) ' ERROR: Could not open input file: ', trim(adjustl(inputfile))
    stop
  end if
  do 
    read(10,"( a200 )",iostat=status) record
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
    else if(keyword(record) == 'binstep') then
      line = value(record)
      read(line,*,iostat=keystatus) binstep
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'nint') then
      line = value(record)
      read(line,*,iostat=keystatus) nintegral
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'cutoff') then
      line = value(record)
      read(line,*,iostat=keystatus) cutoff
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'dbulk') then
      line = value(record)
      read(line,*,iostat=keystatus) dbulk
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
            keyword(record) == 'solvent') then
      write(*,"(a,/,a)") ' ERROR: The options solute and solvent must be used ',&
                         '        with the gss.sh script, not directly. '
      stop
    else if(record(1:1) /= '#' .and. & 
            keyword(record) /= 'par' .and. &
            record(1:1) > ' ') then
      write(*,*) ' ERROR: Unrecognized keyword found: ',keyword(record)
      stop
    end if
  end do               
  close(10)

  ! If some error was found in some keyword value, report error and stop
  
  if(keystatus /= 0) then
    line = keyword(record)
    write(*,*) ' ERROR: Could not read value for keyword: ',line(1:length(line))
    stop
  end if
  
  ! Reading the header of psf file
  
  call getdim(psffile,inputfile,natom)
  allocate( eps(natom), sig(natom), q(natom), e(natom), s(natom), mass(natom),&
            segat(natom), resat(natom), classat(natom), typeat(natom), class(natom),&
            resid(natom), natres(natom), fatres(natom), fatrsolute(natom),&
            fatrsolvent(natom), irsolv(natom), solute2(natom) )

  ! These will contain indexes for the atoms of the randomly generated solvent molecules,
  ! which are more than the number of the atoms of the solvent in the actual
  ! simulation. Therefore, 2*natom is an overstimated guess. If this turns out not
  ! to be sufficient, it will be fixed afterwards.

  allocate( solvent_random(2*natom), irsolv_random(2*natom) )

  ! Reading parameter files to get the vdW sigmas for the definition of exclusion
  ! zones
  
  nclass = 0
  open(10,file=inputfile,action='read',status='old')
  do while(.true.)
    read(10,"( a200 )",iostat=status) record
    if(status /= 0) exit
    if(keyword(record) == 'par') then
      file = value(record)
      write(*,*) ' Reading parameter file: ', file(1:length(file))
      call readpar(file,nclass,class,eps,sig)
    end if
  end do
  close(10)
  
  ! Conversion factor for volumes (as KB integrals), from A^3 to cm^3/mol

  convert = mole / 1.e24

  ! compute ibulk from dbulk (distance from which the solvent is considered bulk,
  ! in the estimation of bulk density)

  if ( dbulk >= cutoff ) then
    write(*,*) ' ERROR: The bulk volume is zero (dbulk >= cutoff). '
    stop
  end if
  if ( dbulk-int(dbulk/binstep)*binstep > 1.e-5 ) then
    write(*,*) ' ERROR: dbulk must be a multiple of binstep. '  
    stop
  end if
  if ( (cutoff-dbulk)-int((cutoff-dbulk)/binstep)*binstep > 1.e-5 ) then
    write(*,*) ' ERROR: (cutoff-dbulk) must be a multiple of binstep. '  
    stop
  end if

  nbins = int(cutoff/binstep)
  ibulk = int(dbulk/binstep) + 1

  write(*,*) ' Width of histogram bins: ', binstep
  write(*,*) ' Number of bins of histograms: ', nbins

  ! Allocate gss array according to nbins

  allocate( gss(nbins), site_count(nbins), site_count_random(nbins), shellvolume(nbins), random_count(nbins) )
  
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
  write(*,*) ' Bulk distance: ', dbulk
  write(*,*) ' Multiplying factor for random count: ', nintegral
  
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

  open(10,file=groupfile,action='read')
  nsolute = 0
  nsolvent = 0
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
    end if     
  end do
  if ( nsolute < 1 .or. nsolvent < 1 ) then
    write(*,*) ' ERROR: No atom selected for solute or solvent. '
    write(*,*) '        nsolute = ', nsolute
    write(*,*) '        nsolvent = ', nsolvent
    stop
  end if
  allocate ( solute(nsolute), solvent(nsolvent) )
  close(10)
  
  ! Now reading reading the group atoms
  
  open(10,file=groupfile,action='read')
  isolute = 0
  isolvent = 0
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
    irsolv(i) = nrsolvent
  end do

  ! This is for the initialization of the smalldistances routine

  nrandom = nintegral*nrsolvent
  maxsmalld = nsolute*nrandom
  allocate( ismalld(maxsmalld), dsmalld(maxsmalld), mind(maxsmalld) )
  maxatom = 2*natom
  allocate( x(maxatom), y(maxatom), z(maxatom), irandom(nrandom) )

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
  write(*,*) ' Number of atoms as specified in the dcd file: ', ntotat     
  call getnframes(10,nframes,dcdaxis,lastframe)
  if(ntotat /= natom) then
    write(*,"(a,/,a)") ' ERROR: Number of atoms in the dcd file does not',&
                      &'        match the number of atoms in the psf file'
    stop
  end if
  
  ! Number of frames (used for normalization of counts)
  
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
  
  ! Reseting the counters
  
  do i = 1, nbins
    site_count(i) = 0.e0
    random_count(i) = 0.e0
    site_count_random(i) = 0.e0
    shellvolume(i) = 0.e0
  end do
  bulkdensity = 0.e0
  simdensity = 0.e0

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

      !
      ! Computing the GSS data the simulation
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
        isolvent = irsolv(ismalld(i))
        if ( dsmalld(i) < mind(isolvent) ) then
          mind(isolvent) = dsmalld(i)
        end if
      end do

      ! Summing up current data to the gss histogram

      do i = 1, nrsolvent
        irad = int(float(nbins)*mind(i)/cutoff)+1
        if( irad <= nbins ) then
          site_count(irad) = site_count(irad) + 1.e0
        end if
      end do

      !
      ! Computing volumes
      !

      ! Solute coordinates are put at the first nsolute positions of x,y,z

      do i = 1, nsolute
        ii = iatom + solute(i)
        x(i) = xdcd(ii)
        y(i) = ydcd(ii)
        z(i) = zdcd(ii)
        solute2(i) = i
      end do

      ! Generate nrandom random points

      ii = nsolute
      do i = 1, nrandom 
        ii = ii + 1
        x(ii) = -axis(1)/2. + random(seed)*axis(1) 
        y(ii) = -axis(2)/2. + random(seed)*axis(2) 
        z(ii) = -axis(3)/2. + random(seed)*axis(3) 
        irandom(i) = ii
      end do

      ! Computes distances which are smaller than the cutoff
      
      call smalldistances(nsolute,solute2,nrandom,irandom,x,y,z,cutoff,&
                          nsmalld,ismalld,dsmalld,axis,maxsmalld)

      do i = 1, nrandom
        mind(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        if ( dsmalld(i) < mind(ismalld(i)) ) then
          mind(ismalld(i)) = dsmalld(i)
        end if
      end do 
      
      ! Estimating volume of bins from count of random points 

      do i = 1, nbins
        random_count(i) = 0.e0
      end do
      do i = 1, nrandom
        irad = int(float(nbins)*mind(i)/cutoff)+1
        if ( irad <= nbins ) then
          random_count(irad) = random_count(irad) + 1.e0
        end if
      end do

      ! Converting counts to volume

      totalvolume = axis(1)*axis(2)*axis(3)
      do i = 1, nbins
        shellvolume(i) = shellvolume(i) + random_count(i)*totalvolume / nrandom
      end do

      ! Estimating the bulk density from site count at large distances

      site_sum = 0.e0
      bulkvolume = 0.e0
      do i = ibulk, nbins
        site_sum = site_sum + site_count(i)
        bulkvolume = bulkvolume + shellvolume(i)
      end do
      bulkdensity = bulkdensity + site_sum/bulkvolume
      simdensity = simdensity + nrsolvent/totalvolume

      !
      ! Bulk density estimated, now we generate a random distribution of solvent
      ! molecules with the appropriate density, to compute the gss in the absence
      ! of solute-solvent interactions: site_sum/bulkvolume is the bulkdensity at
      ! this frame, as computed above
      !

      nrsolvent_random = nint(totalvolume*(site_sum/bulkvolume))

      if ( nrsolvent_random > size( solvent_random ) ) then
        deallocate( solvent_random, irsolv_random ) 
        allocate( solvent_random(nrsolvent_random), irsolv_random(nrsolvent_random) )
      end if

      do isolvent = 1, nrsolvent_random

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

        ! Put molecule in its center of coordinates for it to be the reference coordinate for the 
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
          ii = nsolute + (isolvent-1)*natoms_solvent + i
          x(ii) = xrnd(i)
          y(ii) = yrnd(i)
          z(ii) = zrnd(i)
          solvent_random(ii-nsolute) = ii
          ! Annotate to which molecule this atom pertains
          irsolv_random(ii-nsolute) = isolvent
        end do

      end do

      ! The solute atom was already added to the xyz array for computing volumes, so now
      ! we have only to compute the distances

      call smalldistances(nsolute,solute2,nrsolvent_random,solvent_random,x,y,z,cutoff,&
                          nsmalld,ismalld,dsmalld,axis,maxsmalld)

      !
      ! Computing the gss functions from distance data
      !

      do i = 1, nrsolvent_random
        mind(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        isolvent = irsolv_random(ismalld(i))
        if ( dsmalld(i) < mind(isolvent) ) then
          mind(isolvent) = dsmalld(i)
        end if
      end do

      ! Summing up current data to the gss histogram

      do i = 1, nrsolvent_random
        irad = int(float(nbins)*mind(i)/cutoff)+1
        if ( irad <= nbins ) then
          site_count_random(irad) = site_count_random(irad) + 1.e0
        end if
      end do

      ! Write progress

      write(*,"( 7a,f6.2,'%' )",advance='no') (char(8),i=1,7), 100.*float(kframe)/nfrcycle
  
      iatom = iatom + ntotat
    end do
    write(*,*)
  end do
  close(10)

  ! Averaging results on the number of frames

  bulkdensity = bulkdensity / frames
  simdensity = simdensity / frames
  do i = 1, nbins
    shellvolume(i) = shellvolume(i) / frames
    site_count(i) = site_count(i) / frames
    site_count_random(i) = site_count_random(i) / frames
    if ( site_count_random(i) > 0.e0 ) then
      gss(i) = site_count(i) / site_count_random(i)
    else
      gss(i) = 0.e0
    end if
  end do

  write(*,*)
  write(*,"(a,f12.5)") '  Solvent density in simulation box (sites/A^3): ', simdensity
  write(*,"(a,f12.5)") '  Estimated bulk solvent density (sites/A^3): ', bulkdensity
  write(*,*)
  write(*,"(a,f12.5)") '  Molar volume of solvent simulation box (cc/mol): ', convert/simdensity
  write(*,"(a,f12.5)") '  Molar volume of solvent in bulk (cc/mol): ', convert/bulkdensity

  solutevolume = convert*nrsolvent*(1.e0/simdensity-1e0/bulkdensity)
  write(*,*)
  write(*,"(a,f12.5)") '  Solute partial volume (cc/mol): ', solutevolume

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
             &'# Density of solvent in simulation box (sites/A^3): ',f12.5,/,&
             &'# Density of solvent in bulk (estimated) (sites/A^3): ',f12.5,/,&
             &'# Molar volume of solvent in simulation (cc/mol): ',f12.5,/,&
             &'# Molar volume of solvent in bulk (estimated) (cc/mol): ',f12.5,/,&
             &'#',/,&
             &'# Solute partial volume estimate (cc/mol): ',f12.5,/,&
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
             &simdensity, bulkdensity, convert/simdensity, convert/bulkdensity, solutevolume, &
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  

  bulkerror = 0.e0
  do i = ibulk, nbins
    bulkerror = bulkerror + gss(i)
  end do
  bulkerror = bulkerror / ( nbins-ibulk+1 )
  do i = ibulk, nbins
    sdbulkerror = (bulkerror - gss(i))**2
  end do
  sdbulkerror = sqrt(sdbulkerror/(nbins-ibulk+1))
  write(*,*)
  write(*,"('  Average and standard deviation of bulk-gss: ',f12.5,' +/-',f12.5 )") bulkerror, sdbulkerror 
  write(29,"('#')")
  write(20,"('# Average and standard deviation of bulk-gss: ',f12.5,'+/-',f12.5 )") bulkerror, sdbulkerror 

  ! Output table

  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
  &'#       1  Minimum distance to solute (dmin)',/,&
  &'#       2  GSS normalized by the GSS RAND distribution. ',/,&
  &'#       3  GSS normalized according to spherical volume of radius dmin.',/,&
  &'#       4  Site count for each dmin, averaged over frames',/,&
  &'#       5  Cumulative sum of sites, averaged over the number of frames  ',/,&
  &'#       6  Site count computed from random solvent distribution, averaged over frames.',/,&
  &'#       7  Cumulative sum of sites for the random distribution, averaged over frames.',/,&
  &'#       8  Kirwood-Buff integral (cc/mol) computed from column 2 with volume estimated from col 6 (int V(r)*(gss-1) dr ',/,&
  &'#       9  Kirwood-Buff integral (cc/mol) computed from column 2 with spherical shell volume (int 4*pi*r^2*(gss-1) dr ',/,&
  &'#      10  Spherical shifted minimum distance ')")
  write(20,"('#')")
  write(20,"('#   1-DISTANCE         2-GSS  3-COUNT RAND  3-SITE COUNT   4-SHELL VOL  5-SPHERE VOL  6-GSS/SPHERE      7-KB INT&
            &   8-KB SPHERE      9-RSHIFT')")

  kbint = 0.e0
  kbintsphere = 0.e0
  do i = 1, nbins

    ! KB integral using shell volumes computed
    kbint = kbint + convert*(gss(i)-1.e0)*shellvolume(i)

    ! KB integral using spherical shell volume 
    kbintsphere = kbintsphere + convert*(gss(i)-1.e0)*sphericalshellvolume(i,binstep)

    ! Distance transformation
    rshift = sphereradiusfromshellvolume(shellvolume(i),binstep)

    write(20,"( 10(tr2,f12.7) )") &
    shellradius(i,binstep),&                                      !  1-DISTANCE
    gss(i),&                                                      !  2-GSS
    site_count(i),&                                               !  3-SITE COUNT
    site_count_random(i),&                                        !  4-COUNT RAND
    shellvolume(i),&                                              !  5-SHELL VOL
    sphericalshellvolume(i,binstep),&                             !  6-SPHER VOL
    site_count(i)/(bulkdensity*sphericalshellvolume(i,binstep)),& !  7-GSS/SPHER
    kbint,&                                                       !  8-KB INT
    kbintsphere,&                                                 !  9-KB SPHER
    rshift                                                        ! 10-RSHIFT

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

! Computes the volume of the spherical shell 
! defined within [(i-1)*step,i*step]

real function sphericalshellvolume(i,step)

  implicit none
  integer :: i
  real :: step, rmin
  real, parameter :: fourthirdsofpi = (4./3.)*3.1415925655

  rmin = (i-1)*step
  sphericalshellvolume = fourthirdsofpi*( (rmin+step)**3 - rmin**3 )

end function sphericalshellvolume

! Compute the point in which the radius comprises half of the
! volume of the shell

real function shellradius(i,step)

  implicit none
  integer :: i
  real :: step, rmin

  rmin = (i-1)*step
  shellradius = ( 0.5e0*( (rmin+step)**3 + rmin**3 ) )**(1.e0/3.e0)

end function shellradius

! Computes the radius that corresponds to a spherical shell of
! a given volume

real function sphereradiusfromshellvolume(volume,step)
 
  implicit none
  real :: volume, step, rmin
  real, parameter :: pi = 3.1415925655
  real, parameter :: fourthirdsofpi = (4./3.)*3.1415925655
  
  if ( 3*step*volume - pi*step**4 <= 0.d0 ) then
    sphereradiusfromshellvolume = 0.d0
    return
  end if
  rmin = (sqrt(3*pi)*sqrt(3*step*volume-pi*step**4)-3*pi*step**2)/(6*pi*step)
  sphereradiusfromshellvolume = ( 0.5e0*( volume/fourthirdsofpi + 2*rmin**3 ) )**(1.e0/3.e0)

end function sphereradiusfromshellvolume










