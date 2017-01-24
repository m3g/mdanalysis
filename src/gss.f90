!
! gss: A program to compute gss[1] radial distribution functions from
!      molecular dynamics simulations in NAMD DCD format.
!
!      Important: THIS IS NOT THE CLASSICAL RADIAL DISTRIBUTION
!                 FUNCTION. It is the shape-dependent RDF used
!                 for non-spherical solutes. It will only coincide
!                 with the classical RDF for perfectly spherical
!                 solutes, for instance if single atoms are used
!                 to define the solute and the solvent. 
!                 The normalization of this distribution 
!                 function is more complicated than the normalization
!                 of the radial distribution function for spherical
!                 solutes. Here, the bulk density of the solvent
!                 is estimated by the counting at long distances, and
!                 a random distribution of solvent molecules is used
!                 to estimate the volumes corresponding to each 
!                 minimum distance count. The normalization if done
!                 by dividing the actual count of sites by the
!                 expected non-interacting count estimated from this volume
!                 and the estimated bulk density.
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
  real, parameter :: pi = 4.d0*atan(1.e0)
  integer, parameter :: memory=15000000
  integer :: maxatom 
  integer :: natom, nsolute, nsolvent, isolute, isolvent, isolvent_random, &
             narg, length, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
             j, iframe, icycle, nfrcycle, iatom, k, ii, &
             status, keystatus, iargc, lastatom, nres, nrsolute,&
             nrsolvent, kframe, irad, nbins, natoms_solvent,&
             nsmalld, & 
             frames
  integer :: nrsolvent_random, natsolvent_random
  real :: dbulk, density_fix
  integer :: nbulk, ibulk, nintegral
  real :: site_sum, convert, solute_volume
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: side(memory,3), mass1, mass2, seed, random, axis(3)
  real, parameter :: mole = 6.022140857e23
  real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
          time0, etime, tarray(2),&
          binstep, &
          cutoff,  kbint, kbintsphere, bulkdensity_at_frame
  real :: bulkdensity, totalvolume, bulkvolume, simdensity, solutevolume
  real :: bulkerror, sdbulkerror
  character(len=200) :: groupfile, line, record, value, keyword,&
                        dcdfile, inputfile, psffile, file,&
                        output, lineformat, gssperatom
  character(len=4) :: dummyc
  logical :: readfromdcd, dcdaxis, periodic, onscreenprogress
  real :: shellradius, rshift
  real :: sphericalshellvolume, sphereradiusfromshellvolume

  integer :: jj
  real :: beta, gamma, theta, cmx, cmy, cmz

  ! Allocatable arrays
  
  integer, allocatable :: solute2(:), solvent_random(:)
  integer, allocatable :: solute(:), solvent(:), resid(:), &
                          natres(:), fatres(:), fatrsolute(:), fatrsolvent(:), &
                          irsolv(:), ismalld(:), irsolv_random(:)
  integer, allocatable :: imind(:)

  real, allocatable :: site_count_random(:)
  real, allocatable :: site_count(:), shellvolume(:),&
                       site_count_at_frame(:)
  real, allocatable :: gss(:)
  real, allocatable :: gss_per_atom(:,:)

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
  nintegral = 10
  dbulk = 12.
  cutoff = 14.
  binstep = 0.1e0
  onscreenprogress = .false.

  ! Default output file names

  output = "gss.dat"
  gssperatom = "gssperatom.dat"
  
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
    else if(keyword(record) == 'onscreenprogress') then
      onscreenprogress = .true.
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
    else if(keyword(record) == 'gssperatom') then
      gssperatom = value(record)
      write(*,*) ' GSS per-atom output file name: ', gssperatom(1:length(gssperatom))
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

  allocate( gss(nbins), site_count(nbins), site_count_random(nbins), shellvolume(nbins), &
            site_count_at_frame(nbins) )
  
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

  allocate( ismalld(nsolute*(2*nsolvent)), dsmalld(nsolute*(2*nsolvent)), mind(2*nrsolvent), &
            imind(nrsolvent) )
  maxatom = 2*natom
  allocate( x(maxatom), y(maxatom), z(maxatom) )

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
            xrnd(natoms_solvent), yrnd(natoms_solvent), zrnd(natoms_solvent),& 
            gss_per_atom(natoms_solvent,nbins) )
  
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
    site_count_random(i) = 0.e0
    shellvolume(i) = 0.e0
    do j = 1, natoms_solvent
      gss_per_atom(j,i) = 0.e0
    end do
  end do
  bulkdensity = 0.e0
  simdensity = 0.e0

  ! Reading dcd file and computing the gss function
   
  iframe = 0
  do icycle = 1, ncycles 
   
    if ( onscreenprogress ) then
      write(*,"( t3,'Cycle',t10,i5,tr2,' Reading: ',f6.2,'%')",&
            advance='no') icycle, 0. 
    end if
  
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
      if ( onscreenprogress ) then
        write(*,"( 7a,f6.2,'%' )",advance='no')&
             (char(8),i=1,7), 100.*float(kframe)/nfrcycle
      end if
    end do
    if ( onscreenprogress ) then
      write(*,"(' Computing: ',f6.2,'%')",advance='no') 0.
    end if
  
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

      if ( cutoff > axis(1)/2.e0 .or. &
           cutoff  > axis(2)/2.e0 .or. &
           cutoff > axis(3)/2.d0 ) then
        write(*,*)
        write(*,*) " ERROR: cutoff > periodic_dimension/2 "
        stop
      end if

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
                          nsmalld,ismalld,dsmalld,axis)

      !
      ! Computing the gss functions from distance data
      !

      ! For each solvent residue, get the MINIMUM distance to the solute
    
      do i = 1, nrsolvent
        mind(i) = cutoff + 1.e0
        imind(i) = 0
      end do
      do i = 1, nsmalld
        isolvent = irsolv(ismalld(i))
        if ( dsmalld(i) < mind(isolvent) ) then
          mind(isolvent) = dsmalld(i)
          j = mod(ismalld(i),natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          imind(isolvent) = j
        end if
      end do

      ! Summing up current data to the gss histogram

      do i = 1, nbins
        site_count_at_frame(i) = 0.e0
      end do
      do i = 1, nrsolvent
        irad = int(float(nbins)*mind(i)/cutoff)+1
        if( irad <= nbins ) then
          site_count(irad) = site_count(irad) + 1.e0
          site_count_at_frame(irad) = site_count_at_frame(irad) + 1.e0
          if ( imind(i) > 0 ) then
            gss_per_atom(imind(i),irad) = gss_per_atom(imind(i),irad) + 1.e0
          end if
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

      ! Total volume of the box at this frame

      totalvolume = axis(1)*axis(2)*axis(3)

      ! Very rough (over)estimate of the solute partial volume

      solute_volume = nsolute*(4./3.)*pi*9.

      ! First (over)esimate of the number of solvent molecules that must
      ! be used for normalization

      nrsolvent_random = nintegral*nint(nrsolvent*totalvolume/(totalvolume-solute_volume))

      !
      ! Using this large sampling of random solvent positions, we will compute the
      ! minimum-distance distribution function, which will we first be used only
      ! to estimate the (minimum-distance) volume of the bulk region. Aftwerwards,
      ! the actual normalization is computed from a subset of this solvent molecules
      ! with the better estimate of the (minimum-distance) bulk density. 
      !

      natsolvent_random = natoms_solvent*nrsolvent_random

      ! Checking the size of the arrays and updating if necessary

      if ( natsolvent_random > size( solvent_random ) ) then
        deallocate( solvent_random, irsolv_random ) 
        allocate( solvent_random(natsolvent_random), irsolv_random(natsolvent_random) )
      end if
      if ( nsolute+natsolvent_random > size(x) ) then
        deallocate( x, y, z )
        allocate( x(nsolute+natsolvent_random), &
                  y(nsolute+natsolvent_random), &
                  z(nsolute+natsolvent_random) )
      end if
      if ( natsolvent_random*nsolute > size(dsmalld) ) then
        deallocate( dsmalld, ismalld ) 
        allocate( dsmalld(natsolvent_random*nsolute), ismalld(natsolvent_random*nsolute) )
      end if
      if ( nrsolvent_random > size(mind) ) then
        deallocate( mind )
        allocate( mind(nrsolvent_random) )
      end if

      do isolvent_random = 1, nrsolvent_random

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
        beta = random(seed)*2.e0*pi
        gamma = random(seed)*2.e0*pi
        theta = random(seed)*2.e0*pi
        call compcart(natoms_solvent,xref,yref,zref,xrnd,yrnd,zrnd,&
                      cmx,cmy,cmz,beta,gamma,theta)

        ! Add this molecule to x, y, z arrays

        do i = 1, natoms_solvent
          ii = nsolute + (isolvent_random-1)*natoms_solvent + i
          x(ii) = xrnd(i)
          y(ii) = yrnd(i)
          z(ii) = zrnd(i)
          solvent_random(ii-nsolute) = ii
          ! Annotate to which molecule this atom pertains
          irsolv_random(ii-nsolute) = isolvent_random
        end do

      end do

      ! The solute atom was already added to the xyz array for computing volumes, so now
      ! we have only to compute the distances

      call smalldistances(nsolute,solute2,natsolvent_random,solvent_random,x,y,z,cutoff,&
                          nsmalld,ismalld,dsmalld,axis)

      ! Up to now, we are only interested in the site count in the bulk region, as set
      ! by the user 

      do i = 1, nrsolvent_random
        mind(i) = cutoff + 1.e0
      end do
      do i = 1, nsmalld
        isolvent = irsolv_random(ismalld(i))
        if ( dsmalld(i) < mind(isolvent) ) then
          mind(isolvent) = dsmalld(i)
        end if
      end do
      
      ! Compute the fraction of minimum distance sites found at the bulk volume region

      nbulk = 0
      do i = 1, nrsolvent_random
        if ( mind(i) >= dbulk .and. mind(i) <= cutoff ) then
          nbulk = nbulk + 1
        end if
      end do

      ! The minimum-distance volume of the bulk is, then...

      bulkvolume = (float(nbulk)/nrsolvent_random)*totalvolume 

      ! Therefore, since we have already computed the gss at these distances,
      ! we can estimate the minimum-distance bulk density

      site_sum = 0.e0
      do i = ibulk, nbins
        site_sum = site_sum + site_count_at_frame(i)
      end do
      bulkdensity_at_frame = site_sum/bulkvolume

      ! These are averaged at the end for final report:

      simdensity = simdensity + nrsolvent/totalvolume
      bulkdensity = bulkdensity + bulkdensity_at_frame

      ! With the bulk density properly estimated, we now compute the total
      ! number of random molecules that should actually be used for normalization, 
      ! which should be bulkdensity*totalvolume. Therefore, the site count above
      ! is overstimated because the density is greater than the actual desired 
      ! density by

      density_fix = (bulkdensity_at_frame*totalvolume)/nrsolvent_random

      ! So lets count the sites at each bin distance for the ideal gas distribution
      ! rescaled for the correct density

      do i = 1, nrsolvent_random
        irad = int(float(nbins)*mind(i)/cutoff)+1
        if ( irad <= nbins ) then
          site_count_random(irad) = site_count_random(irad) + 1.e0*density_fix
        end if
      end do

      ! Write progress

      if ( onscreenprogress ) then
        write(*,"( 7a,f6.2,'%' )",advance='no') (char(8),i=1,7), 100.*float(kframe)/nfrcycle
      else
        if ( mod(iframe,max(1,(frames/1000))) == 0 ) then
          write(*,"( '  Progress: ',f6.2,'%' )") 100.*float(iframe)/frames
        end if
      end if
  
      iatom = iatom + ntotat
    end do
    write(*,*)
  end do
  close(10)

  !
  ! Averaging results on the number of frames
  !

  bulkdensity = bulkdensity / frames
  simdensity = simdensity / frames
  do i = 1, nbins
    site_count(i) = site_count(i)/frames
    site_count_random(i) = site_count_random(i)/frames
    shellvolume(i) = site_count_random(i)/bulkdensity
    if ( site_count_random(i) > 0.e0 ) then
      gss(i) = site_count(i)/site_count_random(i)
      do j = 1, natoms_solvent
        gss_per_atom(j,i) = gss_per_atom(j,i)/frames  
        gss_per_atom(j,i) = gss_per_atom(j,i)/site_count_random(i)
      end do
    else
      gss(i) = 0.e0
      do j = 1, natoms_solvent
        gss_per_atom(j,i) = 0.e0
      end do
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
  write(20,"('#')")
  write(20,"('# Average and standard deviation of bulk-gss: ',f12.5,'+/-',f12.5 )") bulkerror, sdbulkerror 

  ! Output table

  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
  &'#       1  Minimum distance to solute (dmin)',/,&
  &'#       2  GSS distribution, normalized by ideal gas distribution. ',/,&
  &'#       3  Site count for each dmin (GSS without normalization).',/,&
  &'#       4  Site count for ideal gas distribution.',/,&
  &'#       5  Shell volume associated with each dmin. ',/,&
  &'#       6  Spherical shell volume with radius dmin. ',/,&
  &'#       7  GSS normalized with spherical shell volume. ',/,&
  &'#       8  Kirwood-Buff integral (cc/mol) computed from column 2. ',/,&
  &'#       9  Kirwood-Buff integral (cc/mol) computed from column 2 with spherical shell volume (int 4*pi*r^2*(gss-1) dr ',/,&
  &'#      10  Spherical-shifted minimum distance ')")
  write(20,"('#')")
  write(20,"('#   1-DISTANCE         2-GSS  3-SITE COUNT  4-COUNT RAND   5-SHELL VOL  6-SPHERE VOL  7-GSS/SPHERE      8-KB INT&
            &   9-KB SPHERE     10-RSHIFT')")

  kbint = 0.e0
  kbintsphere = 0.e0
  do i = 1, nbins

    ! KB integral using shell volumes computed
    kbint = kbint + convert*(gss(i)-1.e0)*shellvolume(i)

    ! KB integral using spherical shell volume 
    kbintsphere = kbintsphere + &
                  convert*(site_count(i)/(bulkdensity*sphericalshellvolume(i,binstep))-1.e0)*sphericalshellvolume(i,binstep)

    ! Distance transformation
    rshift = sphereradiusfromshellvolume(shellvolume(i),binstep)

    lineformat = "(10(tr2,f12.7))"
    if ( abs(kbint) > 999.e0 .or. abs(kbintsphere) > 999.e0 ) then
      lineformat = "(7(tr2,f12.7),2(tr2,e12.5),(tr2,f12.7))"
    end if

    write(20,lineformat) &
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

  ! Writting gss per atom contributions 

  open(20,file=gssperatom)
  write(20,"(a)") "# GSS contribution per atom. "
  write(20,"( '#',/,&
             &'# Input file: ',a,/,& 
             &'# DCD file: ',a,/,& 
             &'# Group file: ',a,/,&
             &'# PSF file: ' )")
  write(20,"(a)") "#"
  write(20,"(a)") "# Atoms: "
  do i = 1, natoms_solvent
    write(20,"( '#', i6, 2(tr2,a), tr2,' mass: ',f12.5 )") i, typeat(solvent(i)), classat(solvent(i)), mass(solvent(i))
  end do
  write(20,"(a)") "#"
  write(lineformat,*) "('#',t7,'DISTANCE     GSS TOTAL',",natoms_solvent,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,natoms_solvent)
  write(lineformat,*) "(",natoms_solvent+2,"(tr2,f12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), gss(i), (gss_per_atom(j,i),j=1,natoms_solvent)
  end do
  close(20)

  ! Write final messages with names of output files and their content
  
  time0 = etime(tarray) - time0
  write(*,*)
  write(*,"( tr2,52('-') )")
  write(*,*)
  write(*,*) ' OUTPUT FILES: ' 
  write(*,*)
  write(*,*) ' Wrote GSS output file: ', trim(adjustl(output))
  write(*,*) ' Wrote GSS per atom output file: ', trim(adjustl(gssperatom))
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
  real, parameter :: pi = 4.d0*atan(1.e0)
  real, parameter :: fourthirdsofpi = (4./3.)*pi

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
  real, parameter :: pi = 4.d0*atan(1.e0)
  real, parameter :: fourthirdsofpi = (4./3.)*pi
  
  if ( 3*step*volume - pi*step**4 <= 0.d0 ) then
    sphereradiusfromshellvolume = 0.d0
    return
  end if
  rmin = (sqrt(3*pi)*sqrt(3*step*volume-pi*step**4)-3*pi*step**2)/(6*pi*step)
  sphereradiusfromshellvolume = ( 0.5e0*( volume/fourthirdsofpi + 2*rmin**3 ) )**(1.e0/3.e0)

end function sphereradiusfromshellvolume










