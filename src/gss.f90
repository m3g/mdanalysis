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
! L. Martinez, S. Shimizu, Minimum distance distribution functions for
! the analysis of the solvation of complex solutes and solvents.
! to be published.
!
! Auxiliar dimensions:
!          memory: Controls the amount of memory used for reading dcd
!                  files. If the program ends with segmentation fault
!                  without any other printing, decrease the value of
!                  this parameter.  
!
! L. Martinez, Mar 13, 2014. (first version)
! Institute of Chemistry, State University of Campinas (UNICAMP)
!
! L. Martinez, Mar 1, 2017 (with KB integrals)
! Institute of Chemistry, State University of Campinas (UNICAMP)
!
! http://leandro.iqm.unicamp.br/mdanalysis
! http://github.com/leandromaratinez98/mdanalysis
!

program g_solute_solvent
 
  ! Static variables

  use file_operations
  implicit none

  real, parameter :: pi = 4.d0*atan(1.e0)
  real, parameter :: mole = 6.022140857e23
  integer, parameter :: memory=15000000

  integer :: i, j, k, ii, jj
  integer :: maxatom
  integer :: natom, nsolute, nsolvent, isolute, isolvent, isolvent_random, &
             narg, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, ntotat, memframes, ncycles, memlast,&
             iframe, icycle, nfrcycle, iatom, &
             status, keystatus, iargc, lastatom, nres, nrsolute,&
             nrsolvent, kframe, ibin, nbins, natoms_solvent,&
             nsmalld, & 
             frames
  integer :: nrsolvent_random, natsolvent_random
  integer :: maxsmalld
  logical :: memerror
  real :: dbulk, density_fix 
  real :: md_int, md_rand_int, volume
  integer :: nbulk, ibulk, nint, nbulk_random, notbulk
  real :: convert
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: side(memory,3), mass1, mass2, random, axis(3)
  real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
          time0, etime, tarray(2),&
          binstep, kbintsphere
  real :: bulkdensity, totalvolume, bulkvolume, simdensity, solutevolume, av_totalvolume, &
          domaindensity
  real :: bulkerror, sdbulkerror
  character(len=200) :: groupfile, line, record, value, keyword,&
                        dcdfile, inputfile, psffile, file,&
                        lineformat
  character(len=200) :: output, output_atom_phi,  &
                        output_atom_kb, output_atom_gss_contrib, output_atom_phi_contrib
  character(len=4) :: dummyc
  logical :: readfromdcd, dcdaxis, periodic, onscreenprogress
  real :: shellradius, rshift
  real :: sphericalshellvolume, sphereradiusfromshellvolume

  real :: beta, gamma, theta, cmx, cmy, cmz

  !
  ! Allocatable arrays
  !
  
  integer, allocatable :: solute2(:), solvent_random(:)
  integer, allocatable :: solute(:), solvent(:), resid(:), &
                          irsolv(:), ismalld(:), irsolv_random(:)
  integer, allocatable :: imind(:)

  ! Shell volume, estimated from atom count
  real, allocatable :: shellvolume(:)

  ! These are the global (whole-solvent-molecule) counts of minimum-distances
  real, allocatable :: md_count(:)
  real, allocatable :: md_count_random(:)
  real, allocatable :: md_atom_contribution(:,:)

  ! This is the resulting phi (md_count/md_count_random)
  real, allocatable :: phi(:) 

  ! This is the resulting gss (md_count/bulkdensity)
  real, allocatable :: gss(:), gss_phantom(:)

  ! KB integral 
  real, allocatable :: kb(:)

  ! These are the counts to compute the phi per atom
  real, allocatable :: site_count_atom(:,:)
  real, allocatable :: site_count_atom_random(:,:)
  real, allocatable :: phi_atom(:,:), kb_atom(:,:)

  ! This is to compute the atomic contributions to phi (md_atom_contribution/md_count_random)
  real, allocatable :: phi_atom_contribution(:,:)

  ! This is to compute the atomic contributions to gss (md_atom_contribution/bulkdensity)
  real, allocatable :: gss_atom_contribution(:,:)

  ! Data read from the psf file, not necessarily used here
  real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:)
  character(len=6), allocatable :: class(:), segat(:), resat(:),&
                                   typeat(:), classat(:)

  ! Array that will contain the coordinates of a single solvent molecule to
  ! be used in random generation
  real, allocatable :: solvent_molecule(:,:)

  ! Coordinates and distances
  real, allocatable :: x(:), y(:), z(:), dsmalld(:), &
                       xref(:), yref(:), zref(:), xrnd(:), yrnd(:), zrnd(:)

  real, allocatable :: mind_mol(:), mind_atom(:)
  
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
  
  call init_random_number(654321)
  
  ! Some default parameters
  
  firstframe = 1
  lastframe = 0
  stride = 1
  periodic = .true.
  readfromdcd = .true.
  nbins = 1000
  nint = 10
  dbulk = 12.
  binstep = 0.1e0
  onscreenprogress = .false.

  ! Default output file names

  output = "gss.dat"

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
      read(line,*,iostat=keystatus) nint
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

  ! Names of auxiliary output files
  
  output_atom_gss_contrib = output(1:length(remove_extension(output)))//"-GSS_ATOM_CONTRIB."//file_extension(output)
  output_atom_phi_contrib = output(1:length(remove_extension(output)))//"-PHI_ATOM_CONTRIB."//file_extension(output)
  output_atom_phi = output(1:length(remove_extension(output)))//"-PHI_PERATOM."//file_extension(output)
  output_atom_kb = output(1:length(remove_extension(output)))//"-KB_PERATOM."//file_extension(output)

  ! Reading the header of psf file
  
  call getdim(psffile,inputfile,natom)
  allocate( eps(natom), sig(natom), q(natom), e(natom), s(natom), mass(natom),&
            segat(natom), resat(natom), classat(natom), typeat(natom), class(natom),&
            resid(natom), &
            irsolv(natom), solute2(natom) )

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

  ! Check compatibility of dbulk and binstep

  if ( dbulk-int(dbulk/binstep)*binstep > 1.e-5 ) then
    write(*,*) ' ERROR: dbulk must be a multiple of binstep. '  
    stop
  end if

  nbins = int(dbulk/binstep)

  write(*,*) ' Width of histogram bins: ', binstep
  write(*,*) ' Number of bins of histograms: ', nbins

  ! Allocate phi array according to nbins

  allocate( phi(nbins), kb(nbins), &
            gss(nbins), gss_phantom(nbins), &
            md_count(nbins), md_count_random(nbins), shellvolume(nbins) )
  
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
  write(*,*) ' Bulk distance: ', dbulk
  write(*,*) ' Multiplying factor for random count: ', nint
  
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

  ! Counting the number of residues of the solute and solvent
  
  j = 0
  nrsolute = 0
  do i = 1, nsolute
    if(resid(solute(i)).gt.j) then
      nrsolute = nrsolute + 1 
      j = resid(solute(i))
    end if
  end do
  j = 0
  nrsolvent = 0
  do i = 1, nsolvent
    if(resid(solvent(i)).gt.j) then
      nrsolvent = nrsolvent + 1 
      j = resid(solvent(i))
    end if
    irsolv(i) = nrsolvent
  end do

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
  
  ! The number of random molecules for numerical normalization 

  nrsolvent_random = nint*nrsolvent
  natsolvent_random = natoms_solvent*nrsolvent_random

  ! Initialization of the smalldistances routine arrays
 
  maxsmalld = nrsolvent_random
  allocate( ismalld(maxsmalld), dsmalld(maxsmalld) )

  ! Allocate xyz and minimum-distance count arrays

  maxatom = max(natom,nsolute+natsolvent_random)
  allocate( x(maxatom), y(maxatom), z(maxatom) )
  allocate( imind(nrsolvent_random), mind_mol(nrsolvent_random), mind_atom(natsolvent_random) )

  ! Allocate solvent molecule (this will be used to generate random coordinates
  ! for each solvent molecule, one at a time, later)
  
  allocate( solvent_molecule(natoms_solvent,3), &
            xref(natoms_solvent), yref(natoms_solvent), zref(natoms_solvent),&
            xrnd(natoms_solvent), yrnd(natoms_solvent), zrnd(natoms_solvent) )
  allocate( site_count_atom(natoms_solvent,nbins), &
            site_count_atom_random(natoms_solvent,nbins), &
            kb_atom(natoms_solvent,nbins), &
            phi_atom(natoms_solvent,nbins) )
  allocate( phi_atom_contribution(natoms_solvent,nbins), & 
            gss_atom_contribution(natoms_solvent,nbins), &
            md_atom_contribution(natoms_solvent,nbins) )

  ! These will contain indexes for the atoms of the randomly generated solvent molecules,
  ! which are more than the number of the atoms of the solvent in the actual
  ! simulation. 

  allocate( solvent_random(natsolvent_random), irsolv_random(natsolvent_random) )

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
  ncycles = (lastframe-1) / memframes + 1
  memlast = lastframe - memframes * ( ncycles - 1 )
  write(*,*) ' Will read and store in memory at most ', memframes,&
             ' frames per reading cycle. '
  write(*,*) ' There will be ', ncycles, ' cycles of reading. '
  write(*,*) ' Last cycle will read ', memlast,' frames. '
  write(*,*)        
  
  ! Reseting the counters
  
  do i = 1, nbins
    md_count(i) = 0.e0
    md_count_random(i) = 0.e0
    shellvolume(i) = 0.e0
    do j = 1, natoms_solvent
      md_atom_contribution(j,i) = 0.e0
      site_count_atom(j,i) = 0.e0
      site_count_atom_random(j,i) = 0.e0
    end do
  end do
  domaindensity = 0.e0
  bulkdensity = 0.e0
  simdensity = 0.e0
  av_totalvolume = 0.e0

  ! Reading dcd file and computing the phi function
   
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
  
    ! Computing the phi function
  
    iatom = 0
    do kframe = 1, nfrcycle
      iframe = iframe + 1

      ! Write progress

      if ( iframe >= firstframe ) then
        if ( onscreenprogress ) then
          write(*,"( 7a,f6.2,'%' )",advance='no') (char(8),i=1,7), 100.*float(kframe)/nfrcycle
        else
          if ( mod(iframe-firstframe,max(1,(frames/1000))) == 0 .and. &
               mod(iframe-firstframe,stride) == 0 ) then
            write(*,"( '  Progress: ',f6.2,'%' )") 100.*float((iframe-firstframe)/stride)/frames
          end if
        end if
      end if

      ! Skip if stride determines so
  
      if(mod(iframe-firstframe,stride) /= 0 .or. iframe < firstframe ) then
        iatom = iatom + ntotat
        cycle
      end if

      ! Sides of the periodic cell in this frame
  
      axis(1) = side(kframe,1) 
      axis(2) = side(kframe,2) 
      axis(3) = side(kframe,3) 

      if ( dbulk > axis(1)/2.e0 .or. &
           dbulk > axis(2)/2.e0 .or. &
           dbulk > axis(3)/2.d0 ) then
        write(*,*)
        write(*,*) " ERROR: dbulk > periodic_dimension/2 "
        stop
      end if

      !
      ! Computing the PHI data the simulation
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

      ! Compute all distances that are smaller than dbulk

      memerror = .true.
      do while ( memerror ) 
        memerror = .false.
        call smalldistances(nsolute,solute,nsolvent,solvent,x,y,z,dbulk,&
                            nsmalld,ismalld,dsmalld,axis,maxsmalld,memerror)
        if ( memerror ) then
          deallocate( ismalld, dsmalld )
          maxsmalld = int(1.5*nsmalld)
          allocate( ismalld(maxsmalld), dsmalld(maxsmalld) )
        end if
      end do

      !
      ! Computing the phi functions from distance data
      !

      ! For each solvent residue, get the MINIMUM distance to the solute
    
      do i = 1, nrsolvent
        mind_mol(i) = dbulk + 1.e0
        imind(i) = 0
      end do
      do i = 1, nsolvent
        mind_atom(i) = dbulk + 1.e0
      end do
      do i = 1, nsmalld
        ! Counting for computing the whole-molecule phi 
        isolvent = irsolv(ismalld(i))
        if ( dsmalld(i) < mind_mol(isolvent) ) then
          mind_mol(isolvent) = dsmalld(i)
          j = mod(ismalld(i),natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          imind(isolvent) = j
        end if
        ! Counting for computing atom-specific phi
        if ( dsmalld(i) < mind_atom(ismalld(i)) ) then
          mind_atom(ismalld(i)) = dsmalld(i)
        end if
      end do

      ! Summing up current data to the phi histogram

      do i = 1, nrsolvent
        ibin = int(float(nbins)*mind_mol(i)/dbulk)+1
        if( ibin <= nbins ) then
          md_count(ibin) = md_count(ibin) + 1.e0
          if ( imind(i) > 0 ) then
            md_atom_contribution(imind(i),ibin) = md_atom_contribution(imind(i),ibin) + 1.e0
          end if
        end if
      end do

      ! Site count at frame, to estimate the bulk density, is performed for a
      ! single solvent reference site, which is taken as atom of type 1 of the solvent

      notbulk = 0
      do i = 1, nsolvent
        ibin = int(float(nbins)*mind_atom(i)/dbulk)+1
        if( ibin <= nbins ) then
          j = mod(i,natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          site_count_atom(j,ibin) = site_count_atom(j,ibin) + 1.e0
          if ( j == 1 ) notbulk = notbulk + 1
        end if
      end do
      nbulk = nrsolvent - notbulk

      ! Total volume of the box at this frame

      totalvolume = axis(1)*axis(2)*axis(3)
      av_totalvolume = av_totalvolume + totalvolume

      ! This is the average density of the solvent in the simulation box, that will
      ! be averaged at the end 

      simdensity = simdensity + nrsolvent/totalvolume

      !
      ! Computing random counts
      !

      ! Solute coordinates are put at the first nsolute positions of x,y,z

      do i = 1, nsolute
        ii = iatom + solute(i)
        x(i) = xdcd(ii)
        y(i) = ydcd(ii)
        z(i) = zdcd(ii)
        solute2(i) = i
      end do

      !
      ! Generating random distribution of solvent molecules in box
      !

      do isolvent_random = 1, nrsolvent_random

        ! First, pick randomly a solvent molecule from the bulk (mind_mol was just computed above
        ! for the actual simulation)
    
        ii = int((nrsolvent-1)*random()) + 1
        do while( mind_mol(ii) < dbulk ) 
          ii = int((nrsolvent-1)*random()) + 1
        end do

        ! Save the coordinates of this molecule in this frame in the solvent_molecule array
    
        jj = iatom + solvent(1) + natoms_solvent*(ii-1)
        do i = 1, natoms_solvent
          solvent_molecule(i,1) = xdcd(jj+i-1)
          solvent_molecule(i,2) = ydcd(jj+i-1)
          solvent_molecule(i,3) = zdcd(jj+i-1)
        end do

        ! Put atom 1 in the center of rotation of the molecule

        do i = 1, natoms_solvent
          xref(i) = solvent_molecule(i,1) - solvent_molecule(1,1)
          yref(i) = solvent_molecule(i,2) - solvent_molecule(1,2)
          zref(i) = solvent_molecule(i,3) - solvent_molecule(1,3)
        end do

        ! Generate a random position for this molecule
  
        cmx = -axis(1)/2. + random()*axis(1) 
        cmy = -axis(2)/2. + random()*axis(2) 
        cmz = -axis(3)/2. + random()*axis(3) 
        beta = random()*2.e0*pi
        gamma = random()*2.e0*pi
        theta = random()*2.e0*pi
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

      memerror = .true.
      do while ( memerror ) 
        memerror = .false.
        call smalldistances(nsolute,solute2,natsolvent_random,solvent_random,x,y,z,dbulk,&
                            nsmalld,ismalld,dsmalld,axis,maxsmalld,memerror)
        if ( memerror ) then
          deallocate( ismalld, dsmalld )
          maxsmalld = int(1.5*nsmalld)
          allocate( ismalld(maxsmalld), dsmalld(maxsmalld) )
        end if
      end do

      ! Counting the number or random molecules with minimum distances
      ! in each region

      do i = 1, nrsolvent_random
        mind_mol(i) = dbulk + 1.e0
        imind(i) = 0
      end do
      do i = 1, nsmalld
        isolvent = irsolv_random(ismalld(i))
        if ( dsmalld(i) < mind_mol(isolvent) ) then
          mind_mol(isolvent) = dsmalld(i)
          j = mod(ismalld(i),natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          imind(isolvent) = j
        end if
      end do

      ! So lets count the sites at each bin distance for the non-interacting distribution 

      do i = 1, nrsolvent_random
        ibin = int(float(nbins)*mind_mol(i)/dbulk)+1
        if ( ibin <= nbins ) then
          md_count_random(ibin) = md_count_random(ibin) + 1.e0
        end if
      end do

      ! Accumulate site count for each atom

      do i = 1, natsolvent_random
        mind_atom(i) = dbulk + 1.e0
      end do
      do i = 1, nsmalld
        if ( dsmalld(i) < mind_atom(ismalld(i)) ) then
          mind_atom(ismalld(i)) = dsmalld(i)
        end if
      end do
      notbulk = 0
      do i = 1, natsolvent_random
        ibin = int(float(nbins)*mind_atom(i)/dbulk)+1
        if ( ibin <= nbins ) then
          j = mod(i,natoms_solvent) 
          if ( j == 0 ) j = natoms_solvent
          site_count_atom_random(j,ibin) = site_count_atom_random(j,ibin) + 1.e0

          ! The counting of single-sites at the bulk region will be used to estimate
          ! the volumes of spherical shells of radius ibin

          if ( j == 1 ) then
            shellvolume(ibin) = shellvolume(ibin) + 1.e0
            notbulk = notbulk + 1
          end if

        end if
      end do
      nbulk_random = nrsolvent_random - notbulk
      if ( nbulk_random == 0 ) then
        write(*,*) 
        write(*,*) ' ERROR: zero volume estimated for bulk region. Either the region is '
        write(*,*) '        too thin, or there is a numerical error. '
        write(*,*) ' frame = ', kframe
        stop
      end if

      ! We have just counted the number of times an atom of type 1 was found
      ! at the bulk region. The minimum-distance volume of the bulk is, then...

      bulkvolume = totalvolume*(float(nbulk_random)/nrsolvent_random)

      ! These are averaged at the end for final report:

      bulkdensity = bulkdensity + float(nbulk)/bulkvolume
      domaindensity = domaindensity + float(nbulk)/bulkvolume

      iatom = iatom + ntotat
    end do
    write(*,*)
  end do
  close(10)

  !
  ! Averaging results on the number of frames
  !

  simdensity = simdensity / frames
  bulkdensity = bulkdensity / frames
  domaindensity = domaindensity / frames
  av_totalvolume = av_totalvolume / frames
  density_fix = (bulkdensity*av_totalvolume)/nrsolvent_random

  write(*,*)
  write(*,"(a,e12.5)") '  Solvent density in simulation box (sites/A^3): ', simdensity
  write(*,"(a,e12.5)") '  Estimated bulk solvent density (sites/A^3): ', bulkdensity
  write(*,"(a,e12.5)") '  Estimated solvent density on solute domain (sites/A^3): ', domaindensity
  write(*,*)
  write(*,"(a,e12.5)") '  Molar volume of solvent in simulation box (cc/mol): ', convert/simdensity
  write(*,"(a,e12.5)") '  Molar volume of solvent in bulk (cc/mol): ', convert/bulkdensity
  write(*,*)
  write(*,"(a,f12.5)") '  Density scaling factor for numerical integration: ', density_fix

  solutevolume = convert*(bulkdensity*av_totalvolume - nrsolvent)/bulkdensity
  write(*,*)
  write(*,"(a,e12.5)") '  Solute partial volume (cc/mol): ', solutevolume
   
  do i = 1, nbins

    md_count(i) = md_count(i)/frames
    md_count_random(i) = density_fix*md_count_random(i)/frames
    do j = 1, natoms_solvent
      md_atom_contribution(j,i) = md_atom_contribution(j,i)/frames
    end do
    shellvolume(i) = ((shellvolume(i)/nrsolvent_random)*av_totalvolume)/frames

    ! PHI distributions

    if ( md_count_random(i) > 0.e0 ) then
      phi(i) = md_count(i)/md_count_random(i)
      do j = 1, natoms_solvent
        phi_atom_contribution(j,i) = md_atom_contribution(j,i)/md_count_random(i)
      end do
    else
      phi(i) = 0.e0
      do j = 1, natoms_solvent
        phi_atom_contribution(j,i) = 0.e0
      end do
    end if

    ! GSS distributions

    if ( shellvolume(i) > 0.e0 ) then
      gss(i) = md_count(i)/(bulkdensity*shellvolume(i))
      do j = 1, natoms_solvent
        gss_atom_contribution(j,i) = md_atom_contribution(j,i)/(bulkdensity*shellvolume(i))
      end do
      gss_phantom(i) = md_count_random(i)/(bulkdensity*shellvolume(i))
    else 
      gss(i) = 0.e0
      do j = 1, natoms_solvent
        gss_atom_contribution(j,i) = 0.e0
      end do
      gss_phantom(i) = natoms_solvent
    end if

    ! Additional distribution required for computing atomic parameters
    
    do j = 1, natoms_solvent
      site_count_atom(j,i) = site_count_atom(j,i) / frames
      site_count_atom_random(j,i) = density_fix*site_count_atom_random(j,i) / frames
      if ( site_count_atom_random(j,i) > 0.e0 ) then
        phi_atom(j,i) = site_count_atom(j,i) / site_count_atom_random(j,i)
      else
        phi_atom(j,i) = 0.e0
      end if
    end do

  end do

  ! Open output file and writes all information of this run

  !
  ! PHI computed with minimum distance
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
             &'# Density of solvent in simulation box (sites/A^3): ',e12.5,/,&
             &'# Density of solvent in bulk (estimated) (sites/A^3): ',e12.5,/,&
             &'# Molar volume of solvent in simulation (cc/mol): ',e12.5,/,&
             &'# Molar volume of solvent in bulk (estimated) (cc/mol): ',e12.5,/,&
             &'#',/,&
             &'# Solute partial volume estimate (cc/mol): ',e12.5,/,&
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

  ! Check the error of the last Angstrom of the computed phi relative to bulk density

  ibulk = int((dbulk-1.e0)/binstep) + 1
  bulkerror = 0.e0
  do i = ibulk, nbins
    bulkerror = bulkerror + phi(i)
  end do
  bulkerror = bulkerror / ( nbins-ibulk+1 )
  do i = ibulk, nbins
    sdbulkerror = (bulkerror - phi(i))**2
  end do
  sdbulkerror = sqrt(sdbulkerror/(nbins-ibulk+1))
  write(*,*)
  write(*,"('  Average and standard deviation long-range bulk-phi: ',f12.5,' +/-',f12.5 )") bulkerror, sdbulkerror 
  write(20,"('# Average and standard deviation of long-range bulk-phi: ',f12.5,' +/-',f12.5 )") bulkerror, sdbulkerror 
  write(20,"('#')")

  ! Check the error of the last Angstrom of the computed gss relative to bulk density

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
  write(*,"('  Average and standard deviation long-range bulk-gss: ',f12.5,' +/-',f12.5 )") bulkerror, sdbulkerror 
  write(20,"('# Average and standard deviation of long-range bulk-gss: ',f12.5,'+/-',f12.5 )") bulkerror, sdbulkerror 
  write(20,"('#')")

  ! Output table

  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
  &'#       1  Minimum distance to solute (dmin)',/,&
  &'#       2  PHI distribution, normalized by phantom-solute distribution. ',/,&
  &'#       3  GSS distribution, normalized by bulk density ',/,&
  &'#       4  GSS PHANTOM distribution, normalized by bulk density ',/,&
  &'#       5  Site count for each dmin (PHI without normalization).',/,&
  &'#       6  Site count for ideal gas distribution.',/,&
  &'#       7  Shell volume associated with each dmin. ',/,&
  &'#       8  Spherical shell volume with radius dmin. ',/,&
  &'#       9  PHI normalized with spherical shell volume. ',/,&
  &'#      10  Kirwood-Buff integral (cc/mol) computed from columns 3 and 4. ',/,&
  &'#      11  Kirwood-Buff integral (cc/mol) computed from column 3 with spherical shell volume (int 4*pi*r^2*(gss-1) dr ',/,&
  &'#      12  Spherical-shifted minimum distance ')")
  write(20,"('#')")
  write(20,"('#   1-DISTANCE         2-PHI         3-GSS 4-GSS PHANTOM  5-SITE COUNT  6-COUNT RAND&
            &   7-SHELL VOL  8-SPHERE VOL  9-PHI/SPHERE     10-KB INT  11-KB SPHERE     12-RSHIFT')")

  do i = 1, nbins
    kb(i) = 0.e0
    do k = 1, natoms_solvent
      kb_atom(k,i) = 0.e0
    end do
  end do
  kbintsphere = 0.e0
  volume = 0.e0
  md_int = 0.e0
  md_rand_int = 0.e0
  do i = 1, nbins

    ! KB integrals 

    md_int = md_int + md_count(i)
    md_rand_int = md_rand_int + md_count_random(i)
    volume = volume + shellvolume(i)

    !kb(i) = convert*(1.e0/bulkdensity)*(md_int - md_rand_int) 
    if ( md_rand_int > 0 ) then
      kb(i) = convert*(volume)*(md_int/md_rand_int - 1.e0) 
    else
      kb(i) = 0.e0
    end if

    do j = i, nbins
      if ( md_count_random(i) > 0.e0 ) then
        kb(j) = kb(j) + convert*(1.e0/bulkdensity)*(md_count(i)-md_count_random(i))
      end if
      do k = 1, natoms_solvent
        if ( site_count_atom_random(k,i) > 0.e0 ) then
          kb_atom(k,j) = kb_atom(k,j) &
                       + convert*(1.e0/bulkdensity)*(site_count_atom(k,i)-site_count_atom_random(k,i))
        end if
      end do
    end do

    ! KB integral using spherical shell volume 
    kbintsphere = kbintsphere + &
                  convert*(md_count(i)/(bulkdensity*sphericalshellvolume(i,binstep))-1.e0)*sphericalshellvolume(i,binstep)

    ! Distance transformation
    rshift = sphereradiusfromshellvolume(shellvolume(i),binstep)

    lineformat = "(12(tr2,f12.7))"
    if ( abs(kb(i)) > 999.e0 .or. abs(kbintsphere) > 999.e0 ) then
      lineformat = "(9(tr2,f12.7),2(tr2,e12.5),(tr2,f12.7))"
    end if

    write(20,lineformat) &
    shellradius(i,binstep),&                                      !  1-DISTANCE
    phi(i),&                                                      !  2-PHI
    gss(i),&                                                      !  3-GSS
    gss_phantom(i),&                                              !  4-GSS PHANTOM
    md_count(i),&                                                 !  5-SITE COUNT
    md_count_random(i),&                                          !  6-COUNT RAND
    shellvolume(i),&                                              !  7-SHELL VOL
    sphericalshellvolume(i,binstep),&                             !  8-SPHER VOL
    md_count(i)/(bulkdensity*sphericalshellvolume(i,binstep)),&   !  9-PHI/SPHER
    kb(i),&                                                       ! 10-KB INT
    kbintsphere,&                                                 ! 11-KB SPHER
    rshift                                                        ! 12-RSHIFT

  end do
  close(20)

  ! Writting gss per atom contributions 

  open(20,file=output_atom_gss_contrib)
  write(20,"(a)") "# Total GSS contribution per atom. "
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
    write(20,lineformat) shellradius(i,binstep), gss(i), (gss_atom_contribution(j,i),j=1,natoms_solvent)
  end do
  close(20)

  ! Writting phi per atom contributions 

  open(20,file=output_atom_phi_contrib)
  write(20,"(a)") "# Total PHI contribution per atom. "
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
  write(lineformat,*) "('#',t7,'DISTANCE     PHI TOTAL',",natoms_solvent,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,natoms_solvent)
  write(lineformat,*) "(",natoms_solvent+2,"(tr2,f12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), phi(i), (phi_atom_contribution(j,i),j=1,natoms_solvent)
  end do
  close(20)

  ! Writting independent atomic PHIs 

  open(20,file=output_atom_phi)
  write(20,"(a)") "# Indepdendent PHI for each atom "
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
  write(lineformat,*) "('#',t7,'DISTANCE     PHI TOTAL',",natoms_solvent,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,natoms_solvent)
  write(lineformat,*) "(",natoms_solvent+2,"(tr2,f12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), phi(i), (phi_atom(j,i),j=1,natoms_solvent)
  end do
  close(20)

  ! Writting independent atomic KB integrals

  open(20,file=output_atom_kb)
  write(20,"(a)") "# Indepdendent KB integrals computed from atom distributions "
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
  write(lineformat,*) "('#',t7,'DISTANCE      KB TOTAL',",natoms_solvent,"(tr2,i12) )"
  write(20,lineformat) (i,i=1,natoms_solvent)
  write(lineformat,*) "(",natoms_solvent+2,"(tr2,e12.5))"
  do i = 1, nbins
    write(20,lineformat) shellradius(i,binstep), kb(i), (kb_atom(j,i),j=1,natoms_solvent)
  end do
  close(20)

  ! Write final messages with names of output files and their content
  
  time0 = etime(tarray) - time0
  write(*,*)
  write(*,"( tr2,52('-') )")
  write(*,*)
  write(*,*) ' OUTPUT FILES: ' 
  write(*,*)
  write(*,*) ' Wrote atomic PHI to file: ', trim(adjustl(output_atom_phi)) 
  write(*,*) ' Wrote atomic KB to file: ', trim(adjustl(output_atom_kb)) 
  write(*,*) ' Wrote atomic GSS contributions to file: ', trim(adjustl(output_atom_gss_contrib))
  write(*,*) ' Wrote atomic PHI contributions to file: ', trim(adjustl(output_atom_phi_contrib))
  write(*,*)
  write(*,*) ' Wrote main output file: ', trim(adjustl(output))
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










