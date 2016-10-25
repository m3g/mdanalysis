!
! This version is the remaining version of GSS 14.255
!
! It computes the dmin_gss; dmax_gss; and davg_gss, but is slow.
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
!                 solutes. Here, a random distribution of solvent
!                 molecules is created, in such a way that the number
!                 of molecules not overlaping the solute coincides 
!                 with the number of molecules of the solvent in the
!                 simulation. The critical parameter is the criterium
!                 for deciding whether two molecules overlap or not.
!                 We use a default overlap distance computed as half
!                 the sum of the vdW radii of the two atoms involved.
!                 This might be controled by the scale_sigma input
!                 parameter, which multiplies that criterium, and
!                 has a default value of scale_sigma = 1.
!                 The way in which we generate the random distribution
!                 of molecules, if we overestimate the number of 
!                 overlaps (if this distance is too large), the solvent
!                 will end up confined in a smaller region of space,
!                 and thus will be denser than in the simulation even
!                 for large distances. On the other side, if the 
!                 overlaps are understimated, the volume available for
!                 solvent molecules will be overestimated and the 
!                 solvent will be less dense than in the simulation. 
!                 Therefore, a good choice for the overlap parameter
!                 is that that provides a unitary gss for large
!                 distances, as we are used. Deviations from unity
!                 are probably caused by poor choices of the overlap
!                 distance.
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
! Version 16.146
!

program g_solute_solvent_old
 
  ! Static variables

  implicit none
  integer, parameter :: memory=15000000
  integer :: natom, nsolute, nsolvent, isolute, isolvent,&
             narg, length, firstframe, lastframe, stride, nclass,&
             nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
             j, iframe, icycle, nfrcycle, iatom, k, ii, jj, catom,&
             status, keystatus, ndim, iargc, lastatom, nres, nrsolute,&
             nrsolvent, kframe, irad, gsssum, nslabs, natoms_solvent,&
             gsssum_random, non_overlap
  double precision :: readsidesx, readsidesy, readsidesz, t
  real :: side(memory,3), sij, d, axis(3), mass1, mass2, seed,&
          cmx, cmy, cmz, beta, gamma, theta, random
  real, parameter :: twopi = 2.*3.1415925655
  real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
          x1, y1, z1, x2, y2, z2, time0, etime, tarray(2),&
          dmin, dmax, davg, gssnorm, gssstep, gssmax, frames,&
          density, scale_sigma, dmin_temp
  character(len=200) :: groupfile, line, record, value, keyword,&
                        dcdfile, inputfile, output, psffile, file
  character(len=4) :: dummyc
  logical :: periodic, readfromdcd, dcdaxis, centeratom, overlap
  
  ! Allocatable arrays
  
  integer, allocatable :: solute(:), solvent(:), resid(:),&
                          natres(:), fatres(:), fatrsolute(:), fatrsolvent(:),&
                          nrsolv(:), gss(:), gss_random(:), gss_max(:), gss_max_random(:),&
                          gss_avg(:), gss_avg_random(:), resnum(:)
  real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:), &
                       solvent_molecule(:,:), &
                       xref(:), yref(:), zref(:), xrnd(:), yrnd(:), zrnd(:),&
                       xsolute(:), ysolute(:), zsolute(:),&
                       sigsolute(:), sigsolvent(:),  &
                       gss_residue(:,:), gss_residue_random(:,:)
  character(len=6), allocatable :: class(:), segat(:), resat(:),&
                                   typeat(:), classat(:)
  logical, allocatable :: tcount(:), counted(:), soluteatom(:), solventatom(:)
  
  ! Compute time
  
  time0 = etime(tarray)
  
  ! Output title
  
  write(*,"(/,' #####################################################' ,&
            &/,/,&
            & '   GSS_MAXMIN: Compute gss distribution from DCD files',&
            &/,/,&
            & ' #####################################################',/)")    
  
  call version
  
  ! Seed for random number generator
  
  seed = 0.48154278727d0
  
  ! Some default parameters
  
  firstframe = 1
  lastframe = 0
  stride = 1
  periodic = .true.
  readfromdcd = .true.
  centeratom = .false.
  catom = 0  
  nslabs = 1000 
  gssmax = 20.
  density = 1.
  scale_sigma = 1.
  
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
    else if(keyword(record) == 'gssmax') then
      line = value(record)
      read(line,*,iostat=keystatus) gssmax
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'scale_sigma') then
      line = value(record)
      read(line,*,iostat=keystatus) scale_sigma
      if(keystatus /= 0) exit 
    else if(keyword(record) == 'density') then
      line = value(record)
      read(line,*,iostat=keystatus) density
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
            fatrsolvent(ndim), tcount(ndim), counted(ndim), nrsolv(ndim),&
            soluteatom(ndim), solventatom(ndim) )
  
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
  
  ! Read PSF file
  
  write(*,*) ' Reading PSF file: ', psffile(1:length(psffile))
  call readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,&
               resid,classat,typeat,q,e,s,mass,.true.)
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
      else if(line(63:66) == '2.00') then
        nsolvent = nsolvent + 1        
      end if
    end if     
  end do
  allocate ( solute(nsolute), solvent(nsolvent), sigsolute(nsolute), resnum(nsolute) )
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
      soluteatom(iatom) = .false.
      solventatom(iatom) = .false.
      if(line(63:66) == '1.00') then      
  
  ! Read atoms belonging to solute
  
        isolute = isolute + 1        
        solute(isolute) = iatom
        soluteatom(iatom) = .true.
        mass1 = mass1 + mass(iatom)
  
      else if(line(63:66) == '2.00') then
  
  ! Read atoms belonging to solvent
  
        isolvent = isolvent + 1        
        solvent(isolvent) = iatom
        solventatom(iatom) = .true.
        mass2 = mass2 + mass(iatom)
  
      end if
    end if     
  end do
  close(10)
  lastatom = max0(solute(nsolute),solvent(nsolvent))
  lastatom = max0(lastatom,catom)

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
    resnum(i) = nrsolute
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

  ! Alocate gss per-residue array

  allocate( gss_residue(nrsolute,nslabs), gss_residue_random(nrsolute,nslabs) )

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

  ! Check if the solute and solvent atoms have obvious reading problems
  
  if(nsolute == 0 .or. nsolvent == 0) then
    write(*,*) ' ERROR: Both solute and solvent must contain at least one atom.'
    stop
  end if
  if ( mod(nsolvent,nrsolvent) /= 0 ) then
    write(*,*) ' ERROR: Incorrect count of solvent atoms or residues. '
    stop
  end if

  natoms_solvent = nsolvent / nrsolvent 
  write(*,*)  ' Number of atoms of each solvent molecule: ', natoms_solvent
  
  ! Allocate solvent molecule (this will be used to generate random coordinates
  ! for each solvent molecule, one at a time, later)
  
  allocate( solvent_molecule(natoms_solvent,3), sigsolvent(natoms_solvent), &
            xref(natoms_solvent), yref(natoms_solvent), zref(natoms_solvent),&
            xrnd(natoms_solvent), yrnd(natoms_solvent), zrnd(natoms_solvent),&
            xsolute(nsolute), ysolute(nsolute), zsolute(nsolute) )
  
  ! Assign the sigma values to sigsolvent and sigsolute vectors
  
  do i = 1, nsolute
    sigsolute(i) = s(solute(i))
  end do
  do i = 1, natoms_solvent
    sigsolvent(i) = s(solvent(i))
  end do
  
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
    do j = 1, nrsolute
      gss_residue(j,i) = 0
      gss_residue_random(j,i) = 0
    end do
  end do
  
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
      write(*,"( 7a,f6.2,'%' )",advance='no')&
           (char(8),i=1,7), 100.*float(j)/nfrcycle
    end do
    write(*,"(' Computing: ',f6.2,'%')",advance='no') 0.
  
    ! Computing the gss function
  
    iatom = 0
    do kframe = 1, nfrcycle
      iframe = iframe + 1
  
      if(mod(iframe,stride) == 0 .and. iframe >= firstframe) then
  
        ! Normalization: creates a box with random coordinates for the solvent molecules, excluding
        ! the solute volume. The gss is computed for this random box independently, and this is 
        ! the normalization for the gss that allows for the computation of the potential of mean
        ! force.
  
        ! Put solute center of coordinates in the origin
  
        cmx = 0.
        cmy = 0.
        cmz = 0.
        do ii = iatom + solute(1), iatom + solute(nsolute)
          if(.not.soluteatom(ii-iatom)) cycle
          cmx = cmx + xdcd(ii)
          cmy = cmy + ydcd(ii)
          cmz = cmz + zdcd(ii)
        end do
        cmx = cmx / float(nsolute)
        cmy = cmy / float(nsolute)
        cmz = cmz / float(nsolute)
        i = 0
        do ii = iatom + solute(1), iatom + solute(nsolute)
          if(.not.soluteatom(ii-iatom)) cycle
          i = i + 1
          xsolute(i) = xdcd(ii) - cmx
          ysolute(i) = ydcd(ii) - cmy
          zsolute(i) = zdcd(ii) - cmz
        end do
  
        ! Create a random box of solvent. Solvent molecules will placed at random
        ! inside the box, until the number of molecules that DO NOT overlap with
        ! the solute is the same as the number of solvent molecules in the actual
        ! simulation. The molecules that overlap are kept. We hope that, in this
        ! way, the number of molecules that would occupy the fraction of the volume 
        ! occupied by the solute will be reasonably estimated.
  
        non_overlap = 0
        do while( non_overlap < nrsolvent )
  
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
  
          ! Generate a random position for this molecule, until none of its atoms overlap with
          ! the solute
  
          cmx = -side(kframe,1)/2. + random(seed)*side(kframe,1) 
          cmy = -side(kframe,2)/2. + random(seed)*side(kframe,2) 
          cmz = -side(kframe,3)/2. + random(seed)*side(kframe,3) 
          beta = random(seed)*twopi
          gamma = random(seed)*twopi
          theta = random(seed)*twopi
          call compcart(natoms_solvent,xref,yref,zref,xrnd,yrnd,zrnd,&
                        cmx,cmy,cmz,beta,gamma,theta)
  
          ! Compute distance between the atoms of this molecule and the solute molecules
  
          overlap = .false.
          dmin = 1.e30
          dmax = 0.
          davg = 0.
          do i = 1, natoms_solvent 
            dmin_temp = 1.e30
            do j = 1, nsolute
              x1 = xrnd(i) - xsolute(j)
              y1 = yrnd(i) - ysolute(j)
              z1 = zrnd(i) - zsolute(j)
              call image(x1,y1,z1,side(kframe,1),side(kframe,2),side(kframe,3))
              d = sqrt( x1**2 + y1**2 + z1**2 )
              dmin_temp = amin1(dmin_temp,d)
              if ( d < ( scale_sigma*(sigsolute(j)+sigsolvent(i)) / 2. ) ) overlap = .true.
            end do
            dmin = amin1(dmin,dmin_temp) 
            dmax = amax1(dmax,dmin_temp)
            davg = davg + dmin_temp
          end do
          davg = davg / natoms_solvent
  
          ! If this molecule is distant from the solute, increase non_overlap
  
          if ( .not. overlap ) non_overlap = non_overlap + 1
  
          ! Count this minimum distance in the gss_random histogram
  
          if( dmin < gssmax ) then
            irad = int(float(nslabs)*dmin/gssmax)
            gss_random(irad) = gss_random(irad) + 1
          end if
          if( dmax < gssmax ) then
            irad = int(float(nslabs)*dmax/gssmax)
            gss_max_random(irad) = gss_max_random(irad) + 1
          end if
          if( davg < gssmax ) then
            irad = int(float(nslabs)*davg/gssmax)
            gss_avg_random(irad) = gss_avg_random(irad) + 1
          end if
  
        end do
  
        ! Now the actual GSS from the simulation will be computed
  
        ! Cycle over solvent atoms
  
        jj = iatom + solvent(1)
        solvent_atom: do while( jj < iatom + solvent(nsolvent) ) 
  
          ! If this atom is not part of solvent selection, cycle
  
          if(.not.solventatom(jj-iatom)) cycle solvent_atom
  
          dmin = 1.e30
          dmax = 0.
          davg = 0.
          do j = 1, natoms_solvent
            dmin_temp = 1.e30
            solute_atom: do ii = iatom + solute(1), iatom + solute(nsolute)
              
              ! If this atom is not part of solute selection, cycle
  
              if(.not.soluteatom(ii-iatom)) cycle solute_atom
  
              ! Move atom ii if centeratom option is set 
  
              if(centeratom) then
                x1 = xdcd(ii) - xdcd(iatom+catom)
                y1 = ydcd(ii) - ydcd(iatom+catom)
                z1 = zdcd(ii) - zdcd(iatom+catom)
                if(readfromdcd) then
                  call image(x1,y1,z1,side(kframe,1),side(kframe,2),side(kframe,3))
                else
                  call image(x1,y1,z1,axis(1),axis(2),axis(3))
                end if
              else
                x1 = xdcd(ii)
                y1 = ydcd(ii)
                z1 = zdcd(ii)
              end if
   
              ! Move atom jj according to reference (either atom ii or centeratom)
  
              if(periodic) then
                if(centeratom) then
                  x2 = xdcd(jj) - xdcd(iatom+catom)
                  y2 = ydcd(jj) - ydcd(iatom+catom)
                  z2 = zdcd(jj) - zdcd(iatom+catom)
                  if(readfromdcd) then
                    call image(x2,y2,z2,side(kframe,1),side(kframe,2),side(kframe,3))
                  else
                    call image(x2,y2,z2,axis(1),axis(2),axis(3))
                  end if
                else
                  x2 = xdcd(jj) - x1
                  y2 = ydcd(jj) - y1
                  z2 = zdcd(jj) - z1
                  if(readfromdcd) then
                    call image(x2,y2,z2,side(kframe,1),side(kframe,2),side(kframe,3))
                  else
                    call image(x2,y2,z2,axis(1),axis(2),axis(3))
                  end if
                end if
                x2 = x2 + x1
                y2 = y2 + y1
                z2 = z2 + z1
              else
                x2 = xdcd(jj)
                y2 = ydcd(jj)
                z2 = zdcd(jj)
              end if
  
              ! Computing distance between atom ii and jj
              
              sij = sqrt( ( x2 - x1 )**2 +&
                          ( y2 - y1 )**2 +&
                          ( z2 - z1 )**2 )

              dmin_temp = amin1(dmin_temp,sij)

            end do solute_atom

            dmin = amin1(dmin_temp,dmin)
            dmax = amax1(dmin_temp,dmax)
            davg = davg + dmin_temp

            jj = jj + 1
          end do
          davg = davg / natoms_solvent
  
          if( dmin < gssmax ) then
            irad = int(float(nslabs)*dmin/gssmax)
            gss(irad) = gss(irad) + 1
          end if
          if( dmax < gssmax ) then
            irad = int(float(nslabs)*dmax/gssmax)
            gss_max(irad) = gss_max(irad) + 1
          end if
          if ( davg < gssmax ) then
            irad = int(float(nslabs)*davg/gssmax)
            gss_avg(irad) = gss_avg(irad) + 1
          end if
  
        end do solvent_atom
      end if
  
      ! Printing progress
  
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
  
  open(20,file="dmin_"//output(1:length(output)))
  write(20,"( '#',/,&
             &'# Output of gss.f90: Using MINIMUM distance to solute.',/,&
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
             &'#' )" )&
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, centeratom, catom,&
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  
  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
             &'#       1  Minimum distance to solute (dmin)',/,&
             &'#       2  GSS normalized by the GSS RAND distribution. ',/,&
             &'#       3  GSS normalized according to spherical volume of radius dmin.',/,&
             &'#       4  GSS not normalized at all (just site count for each dmin)',/,&
             &'#       5  Cumulative sum of sites ',/,&
             &'#       6  GSS computed from random solvent distribution, not normalized ',/,&
             &'#       7  Cumulative sum of sites for the random distribution. ')")
  write(20,"( '#',/,&      
    '#',t5,'1-DISTANCE',t17,'2-GSS/GSSRND',t32,'3-GSS/SPHER',t52,'4-GSS',t64,'5-CUMUL',&
        t76,'6-GSS RND',t88,'7-CUMUL RND' )" )
  
  frames=float(lastframe-firstframe+1)/float(stride)
  gsssum = 0
  gsssum_random = 0
  do i = 1, nslabs
    gsssum = gsssum + gss(i)
    gsssum_random = gsssum_random + gss_random(i)
    if(i > 1) then
      gssnorm = gss(i)/(frames*2.*twopi*((i-1)*gssstep)**2*gssstep)
    else
      gssnorm = 0.
    end if
    x1 = float(gss(i))/frames
    y1 = float(gss_random(i))/frames
    if ( y1 > 0. ) then
      z1 = x1 / y1
    else
      z1 = 0.
    end if
    write(20,"( 4(tr2,f12.7),tr2,i12,tr2,f12.7,tr2,i12 )")&
    i*gssstep-gssstep/2., z1, gssnorm, float(gss(i))/frames, int(gsssum/frames), &
                          float(gss_random(i))/frames, int(gsssum_random/frames)
  end do
  close(20)
  
  ! 
  ! GSS computed with maximum distance
  !
  
  open(20,file="dmax_"//output(1:length(output)))
  write(20,"( '#',/,&
             &'# Output of gss.f90: Using MAXIMUM distance to solute.',/,&
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
             &'#' )" )&
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, centeratom, catom,&
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  
  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
             &'#       1  Maximum distance to solute (dmax)',/,&
             &'#       2  GSS normalized by the GSS RAND distribution. ',/,&
             &'#       3  GSS normalized according to spherical volume of radius dmax.',/,&
             &'#       4  GSS not normalized at all (just site count for each dmax)',/,&
             &'#       5  Cumulative sum of sites ',/,&
             &'#       6  GSS computed from random solvent distribution, not normalized ',/,&
             &'#       7  Cumulative sum of sites for the random distribution. ')")
  write(20,"( '#',/,&      
    '#',t5,'1-DISTANCE',t17,'2-GSS/GSSRND',t32,'3-GSS/SPHER',t52,'4-GSS',t64,'5-CUMUL',&
        t76,'6-GSS RND',t88,'7-CUMUL RND' )" )
  
  frames=float(lastframe-firstframe+1)/float(stride)
  gsssum = 0
  gsssum_random = 0
  do i = 1, nslabs
    gsssum = gsssum + gss_max(i)
    gsssum_random = gsssum_random + gss_max_random(i)
    if(i > 1) then
      gssnorm = gss_max(i)/(frames*2.*twopi*((i-1)*gssstep)**2*gssstep)
    else
      gssnorm = 0.
    end if
    x1 = float(gss_max(i))/frames
    y1 = float(gss_max_random(i))/frames
    if ( y1 > 0. ) then
      z1 = x1 / y1
    else
      z1 = 0.
    end if
    write(20,"( 4(tr2,f12.7),tr2,i12,tr2,f12.7,tr2,i12 )")&
    i*gssstep-gssstep/2., z1, gssnorm, float(gss_max(i))/frames, int(gsssum/frames), &
                          float(gss_max_random(i))/frames, int(gsssum_random/frames)
  end do
  close(20)

  !
  ! GSS computed from average distance 
  !
  
  open(20,file="davg_"//output(1:length(output)))
  write(20,"( '#',/,&
             &'# Output of gss.f90: Using AVERAGE distance to solute.',/,&
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
             &'#' )" )&
             &inputfile(1:length(inputfile)),&
             &dcdfile(1:length(dcdfile)),&
             &groupfile(1:length(groupfile)),&
             &psffile(1:length(psffile)),&
             &firstframe, lastframe, stride,&
             &periodic, readfromdcd, centeratom, catom,&
             &nsolute, mass1, solute(1), solute(nsolute),& 
             &nsolvent, mass2, solvent(1), solvent(nsolvent)  
  write(20,"( '# COLUMNS CORRESPOND TO: ',/,&
             &'#       1  Average distance to solute (davg)',/,&
             &'#       2  GSS normalized by the GSS RAND distribution. ',/,&
             &'#       3  GSS normalized according to spherical volume of radius davg.',/,&
             &'#       4  GSS not normalized at all (just site count for each davg)',/,&
             &'#       5  Cumulative sum of sites ',/,&
             &'#       6  GSS computed from random solvent distribution, not normalized ',/,&
             &'#       7  Cumulative sum of sites for the random distribution. ')")
  write(20,"( '#',/,&      
    '#',t5,'1-DISTANCE',t17,'2-GSS/GSSRND',t32,'3-GSS/SPHER',t52,'4-GSS',t64,'5-CUMUL',&
        t76,'6-GSS RND',t88,'7-CUMUL RND' )" )
  
  frames=float(lastframe-firstframe+1)/float(stride)
  gsssum = 0
  gsssum_random = 0
  do i = 1, nslabs
    gsssum = gsssum + gss_avg(i)
    gsssum_random = gsssum_random + gss_avg_random(i)
    if(i > 1) then
      gssnorm = gss_avg(i)/(frames*2.*twopi*((i-1)*gssstep)**2*gssstep)
    else
      gssnorm = 0.
    end if
    x1 = float(gss_avg(i))/frames
    y1 = float(gss_avg_random(i))/frames
    if ( y1 > 0. ) then
      z1 = x1 / y1
    else
      z1 = 0.
    end if
    write(20,"( 4(tr2,f12.7),tr2,i12,tr2,f12.7,tr2,i12 )")&
    i*gssstep-gssstep/2., z1, gssnorm, float(gss_avg(i))/frames, int(gsssum/frames), &
                          float(gss_avg_random(i))/frames, int(gsssum_random/frames)
  end do
  close(20)

  ! Write final messages with names of output files and their content
  
  time0 = etime(tarray) - time0
  write(*,*)
  write(*,"( tr2,52('-') )")
  write(*,*)
  write(*,*) ' OUTPUT FILES: ' 
  write(*,*)
  write(*,*) ' Wrote output files: ' 
  write(*,*) ' dmin_'//output(1:length(output))
  write(*,*) ' dmax_'//output(1:length(output))
  write(*,*) ' davg_'//output(1:length(output))
  write(*,*)
  write(*,*) ' Which contain gss functions computed with the '
  write(*,*) ' 1-MINIMUM distance between solute and solvent (standard) '
  write(*,*) ' 2-The maximum-minimum distance. '
  write(*,*) ' 3-The average-minimum distance. ' 
  write(*,*) 
  write(*,*) ' Running time: ', time0
  write(*,*) '####################################################'
  write(*,*) 
  write(*,*) '  END: Normal termination.  '
  write(*,*) 
  write(*,*) '####################################################'
  write(*,*)        

end program g_solute_solvent_old

