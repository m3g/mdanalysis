!
! align: A program for analysing mobility between selected groups
!        in molecular dynamics simulations in NAMD DCD format.
!
! Auxiliar dimensions:
!          memory: Controls the ammount of memory used for reading dcd
!                  files. If the program ends with segmentation fault
!                  without any other printing, decrease the value of
!                  this parameter.  
!
! L. Martinez, Mar 17, 2008.
!
! Dynamically allocatable arrays and their dimensions:
!
!       dcdfile: number of dcdfiles
!       nfdcd: number of frames per dcd file
!       i_reference, xref1: number of atoms of reference
!       i_rmsd, xref2: number of atoms of "compute_rmsd_of"
!       eps, sig, q, e, s, class, segat, resat, typeat, classat, mass,
!               xref, yref, zref : total number of atoms of the system
!       xm ... zp: Auxiliar vectors for the align subroutine: maximum
!                  between n_reference and n_rmsd
!       xfirst...zfirst: Auxiliar vectors for computing
!                        the average structure, only allocated
!                        in the case of using it.
!       rmsd_residue: Array that contains the RMSD per residue.
!       i_residue: Array containing the sequential numbers of the 
!                  different residues
!

! Static variables
 
use charsize
implicit none
integer, parameter :: memory=15000000
integer :: natom, n_reference, n_rmsd,&
           narg, length, firstframe, lastframe, stride, nclass,&
           nframes, dummyi, i, ntotat, memframes, ncycles, memlast,&
           j, iframe, icycle, nfrcycle, iatom, k, ii, ic,&
           status, keystatus, ndim, iargc, alignto, lastatom,&
           number_of_residues, last_rmsd_atom, ndcd, idcd, ifcount
double precision :: readsidesx, readsidesy, readsidesz, t
real :: side(memory,3), cm1x, cm1y, cm1z, cm2x, cm2y, cm2z,&
        cm1rx, cm1ry, cm1rz,&
        cm2rx, cm2ry, cm2rz,&
        mass1, mass2, rmsdint, rmsd, rmsdnow,&
        xmove, ymove, zmove, u(3,3), max_rmsd, min_rmsd
real :: dummyr, xdcd(memory), ydcd(memory), zdcd(memory),&
        time0, etime, tarray(2), relative_rmsd, scaletime
character(len=200) :: groupfile, line, record, value, keyword,&
                      inputfile, output, psffile, record2,&
                      pdboutinter, pdboutrel, pdbfile, write_per_residue_file
character(len=4) :: dummyc
logical :: dcdaxis, writeinter, writerel, write_per_residue

! Allocatable arrays

integer, allocatable :: i_reference(:), i_rmsd(:), resid(:), i_residue(:),&
                        natoms_per_residue(:), nfdcd(:)
real, allocatable :: eps(:), sig(:), q(:), e(:), s(:), mass(:),&
                     xm(:), ym(:), zm(:), xp(:), yp(:), zp(:)
real, allocatable :: xref1(:), yref1(:), zref1(:),&
                     xref2(:), yref2(:), zref2(:),&
                     xrmsd_ini(:), yrmsd_ini(:), zrmsd_ini(:),&
                     xfirst1(:), yfirst1(:), zfirst1(:),&
                     xfirst2(:), yfirst2(:), zfirst2(:),&
                     rmsd_residue(:), rmsd_residue_frame(:)
character(len=charsize1), allocatable :: class(:), segat(:), resat(:),&
                                        typeat(:), classat(:)
character(len=200), allocatable :: dcdfile(:)

! Compute time

time0 = etime(tarray)

! Output title

write(*,"(/,' ####################################################',&
          &/,/,&
          & '   ALIGN: Compute RMSD from DCD files  ',&
          &/,/,&
          & ' ####################################################',/)")    

call version()

! Some default parameters

firstframe = 1
lastframe = 0
stride = 1
alignto = 1
writerel = .false.
writeinter = .false.
write_per_residue = .false.
scaletime = 1.

! Open input file and read parameters

narg = iargc()
if(narg == 0) then
  write(*,*) ' Run with: ./align input.inp '
  stop
end if   
call getarg(1,record)

inputfile = record(1:length(record))
open(99,file=inputfile,action='read')

! Read the number of dcd files set

ndcd = 0
do 
  read(99,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(keyword(record) == 'dcd') ndcd = ndcd + 1
end do
if ( ndcd < 1 ) then
  write(*,*) ' ERROR: No DCD trajectory file set. '
  stop
end if
allocate(dcdfile(ndcd),nfdcd(ndcd))
ndcd = 0
rewind(99)

! Now, reading all parameters 

do 
  read(99,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(keyword(record) == 'dcd') then
    ndcd = ndcd + 1
    dcdfile(ndcd) = value(record)
    record2 = dcdfile(ndcd)
    write(*,*) ' DCD file ', ndcd,' name: ', record2(1:length(record2))
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
  else if(keyword(record) == 'writeinter') then
    writeinter = .true.
    pdboutinter = value(record)
    write(*,*) ' PDB output file of internal motions: ',&
               pdboutinter(1:length(pdboutinter))
  else if(keyword(record) == 'writerel') then
    writerel = .true.
    pdboutrel = value(record)
    write(*,*) ' PDB output file of relative motions: ',&
               pdboutrel(1:length(pdboutrel))
  else if(keyword(record) == 'write_per_residue') then
    write_per_residue = .true.
    write_per_residue_file = value(record)
    write(*,*) ' File with RMSD per residue: ',&
               write_per_residue_file(1:length(write_per_residue_file))
  else if(keyword(record) == 'alignto') then
    if(value(record) == 'dcdfirst') then
      alignto = 1
    else if(value(record) == 'dcdlast') then
      alignto = 0
    else if(value(record) == 'dcdaverage') then
      alignto = -1
    else
      line = value(record)
      read(line,"( i20 )",iostat=keystatus) alignto
      if(keystatus /= 0) then
        alignto = -2
        pdbfile=value(record)
        write(*,*) ' PDB reference file: ', pdbfile(1:length(pdbfile))
      end if
    end if
  else if(keyword(record) == 'output') then
    output = value(record)
    write(*,*) ' Output file name: ', output(1:length(output))
  else if(keyword(record) == 'compute_rmsd_of' .or. &
          keyword(record) == 'reference') then
    write(*,"(a,/,a)") ' ERROR: The options reference and compute_rmsd_of',&
                       '        must be used with the align.sh script, not directly. '
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
          resid(ndim), segat(ndim), resat(ndim), classat(ndim),&
          typeat(ndim), class(ndim) )
       
! Check for simple input errors

if(stride < 1) then
  write(*,*) ' ERROR: stride cannot be less than 1. ' 
  stop
end if
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
write(*,*) ' Stride (will jump frames): ', stride
write(*,*) ' Time (x-axis) will be scaled by: ', scaletime

! Read PSF file

write(*,*) ' Reading PSF file: ', psffile(1:length(psffile))
call readpsf(psffile,nclass,class,eps,sig,natom,segat,resat,resid,&
             classat,typeat,q,e,s,mass,.false.)
write(*,*) ' Number of atoms in PSF file: ', natom
write(*,*)

! Read group information from group file
! First reading the size of the groups to allocate arrays

open(10,file=groupfile,action='read')
n_reference = 0
n_rmsd = 0
number_of_residues = 1
iatom = 0
last_rmsd_atom = 1
do 
  read(10,"( a200 )",iostat=status) line
  if(status /= 0) exit
  if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
    iatom = iatom + 1
    if(line(57:60) == '1.00') then      
      n_reference = n_reference + 1        
    end if
    if(line(63:66) == '2.00') then
      n_rmsd = n_rmsd + 1        

! Computing the number of residues of the structure

      if(n_rmsd > 1) then
        if(resid(iatom) /= resid(last_rmsd_atom)) then
          number_of_residues = number_of_residues + 1
          last_rmsd_atom = iatom
        end if
      end if

    end if
  end if     
end do
close(10)

allocate ( i_reference(max0(n_reference,n_rmsd)),&
           i_rmsd(max0(n_reference,n_rmsd)),&
           xref1(max0(n_reference,n_rmsd)),&
           yref1(max0(n_reference,n_rmsd)),&
           zref1(max0(n_reference,n_rmsd)),&
           xref2(max0(n_reference,n_rmsd)),&
           yref2(max0(n_reference,n_rmsd)),&
           zref2(max0(n_reference,n_rmsd)),&
           xrmsd_ini(max0(n_reference,n_rmsd)),&
           yrmsd_ini(max0(n_reference,n_rmsd)),&
           zrmsd_ini(max0(n_reference,n_rmsd)),&
           xm(max0(n_reference,n_rmsd)),& 
           ym(max0(n_reference,n_rmsd)),&
           zm(max0(n_reference,n_rmsd)),& 
           xp(max0(n_reference,n_rmsd)),& 
           yp(max0(n_reference,n_rmsd)),&
           zp(max0(n_reference,n_rmsd)),&
           rmsd_residue(number_of_residues),&
           rmsd_residue_frame(number_of_residues),&
           natoms_per_residue(number_of_residues),&
           i_residue(n_rmsd))

! Now reading reading the group atoms

open(10,file=groupfile,action='read')
n_reference = 0
n_rmsd = 0
mass1 = 0.d0
mass2 = 0.d0
iatom = 0
do i = 1, number_of_residues
  natoms_per_residue(i) = 0
end do
number_of_residues = 1

do 
  read(10,"( a200 )",iostat=status) line
  if(status /= 0) exit
  if(line(1:4) == 'ATOM' .or. line(1:6) == 'HETATM') then 
    iatom = iatom + 1

  ! Read atoms belonging to reference

    if(line(57:60) == '1.00') then      
      n_reference = n_reference + 1        
      i_reference(n_reference) = iatom
      mass1 = mass1 + mass(iatom)
    end if

    ! Read atoms belonging to compute_rmsd_of

    if(line(63:66) == '2.00') then
      n_rmsd = n_rmsd + 1        
      i_rmsd(n_rmsd) = iatom
      mass2 = mass2 + mass(iatom)
      if ( n_rmsd > 1 ) then
        if (resid(i_rmsd(n_rmsd)) /= resid(i_rmsd(n_rmsd-1))) then
          number_of_residues = number_of_residues + 1 
          natoms_per_residue(number_of_residues) = natoms_per_residue(number_of_residues) + 1 
        else
          natoms_per_residue(number_of_residues) = natoms_per_residue(number_of_residues) + 1 
        end if
      else
        natoms_per_residue(1) = natoms_per_residue(1) + 1 
      end if
      i_residue(n_rmsd) = number_of_residues
    end if

  end if     
end do
close(10)
lastatom = max0(i_reference(n_reference),i_rmsd(n_rmsd))

! Output some group properties for testing purposes

write(*,*) ' Number of atoms of group 1 (reference): ', n_reference 
write(*,*) ' First atom of group 1: ', i_reference(1)
write(*,*) ' Last atom of group 1: ', i_reference(n_reference)
write(*,*) ' Mass of group 1: ', mass1
write(*,*) ' Number of atoms of group 2 (compute_rmsd_of): ', n_rmsd 
write(*,*) ' First atom of group 2: ', i_rmsd(1)
write(*,*) ' Last atom of group 2: ', i_rmsd(n_rmsd)
write(*,*) ' Mass of group 2: ', mass2
write(*,*) ' Number of residues of group 2: ', number_of_residues
 
! Checking if dcd file contains periodic cell information

write(*,"( /,tr2,52('-') )")
write(*,*) ' Periodic cell data: Read carefully. '
call chkperiod(dcdfile(1),dcdaxis,.false.) 

! Stop if some group is empty

if(n_reference == 0 .or. n_rmsd == 0) then
  write(*,*) ' ERROR: Groups must contain at least one atom.'
  stop
end if

! Reading the total number of frames of all dcd files

nframes = 0
do idcd = 1, ndcd
  open(10,file=dcdfile(idcd),action='read',form='unformatted')
  call getnframes(10,nfdcd(idcd),dcdaxis,lastframe)
  close(10)
  write(*,*) ' Number of frames of DCD file ', idcd,':', nfdcd(idcd)
  nframes = nframes + nfdcd(idcd)
end do
if(lastframe == 0) lastframe = nframes

! Getting the structure relative to which the alignment will be made

! If the reference is the external pdbfile provided:

if(alignto == -2) then
  open(10,file=pdbfile,action='read')
  iatom = 0
  do while(iatom.lt.lastatom)
    read(10,"( a200 )",iostat=status) line
    if(status /= 0) then
      write(*,*) ' ERROR reading reference PDB file. '
      stop
    end if
    if(line(1:4) == "ATOM" .or. line(1:6) == "HETATM") then
      iatom = iatom + 1
      read(line(31:38),*,iostat=status) xdcd(iatom)
      read(line(39:46),*,iostat=status) ydcd(iatom)
      read(line(47:54),*,iostat=status) zdcd(iatom)
      if(status /= 0) then
        write(*,*) ' ERROR reading reference PDB file. '
        stop
      end if
    end if
  end do
  do i = 1, n_reference
    xref1(i) = xdcd(i_reference(i))
    yref1(i) = ydcd(i_reference(i))
    zref1(i) = zdcd(i_reference(i))
  end do
  do i = 1, n_rmsd
    xref2(i) = xdcd(i_rmsd(i))
    yref2(i) = ydcd(i_rmsd(i))
    zref2(i) = zdcd(i_rmsd(i))
    xrmsd_ini(i) = xdcd(i_rmsd(i))
    yrmsd_ini(i) = ydcd(i_rmsd(i))
    zrmsd_ini(i) = zdcd(i_rmsd(i))
  end do 
  close(10)

! If the reference is based on DCD coordinates

else
  write(*,"(tr2,52('-'),/ )")

  if(alignto >= 0) then

    if(alignto == 1) then
      write(*,*) ' Reading DCD reference structure for alignment: first frame'
    else if(alignto == 0) then
      alignto = nframes
      write(*,*) ' Reading DCD reference structure for alignment: last frame:',&
                 alignto
      write(*,*) ' This may take a while, since must read whole DCD file.'
    else if(alignto > 1) then
      write(*,*) ' Reading DCD reference structure for alignment: ', alignto
      write(*,*) ' This may take a while if the structure is at the end.'
    end if

    ifcount = 0
    do idcd = 1, ndcd

      ! Skip DCD files if the first frame is not in it

      if( alignto > ifcount + nfdcd(idcd) ) then
        ifcount = ifcount + nfdcd(idcd)
        cycle
      end if

      open(10,file=dcdfile(idcd),action='read',form='unformatted')
      read(10) dummyc, dummyi, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
      read(10) dummyi, dummyr
      read(10) ntotat
      if(ntotat /= natom) then
        write(*,"(a,/,a)") ' ERROR: Number of atoms in the dcd file does not',&
                          &'        match the number of atoms in the psf file'
        stop
      end if

      ! If the reference frame is some frame of the dcd file

      do i = 1, alignto - ifcount - 1
        if(dcdaxis) read(10) 
        read(10) 
        read(10) 
        read(10) 
      end do
      if(dcdaxis) read(10)
      read(10) (xdcd(k), k = 1, lastatom)
      read(10) (ydcd(k), k = 1, lastatom)            
      read(10) (zdcd(k), k = 1, lastatom)

      ! Save reference coordinates in the x-reference vectors

      do i = 1, n_reference
        xref1(i) = xdcd(i_reference(i))
        yref1(i) = ydcd(i_reference(i))
        zref1(i) = zdcd(i_reference(i))
      end do
      do i = 1, n_rmsd
        xref2(i) = xdcd(i_rmsd(i))
        yref2(i) = ydcd(i_rmsd(i))
        zref2(i) = zdcd(i_rmsd(i))
        xrmsd_ini(i) = xdcd(i_rmsd(i))
        yrmsd_ini(i) = ydcd(i_rmsd(i))
        zrmsd_ini(i) = zdcd(i_rmsd(i))
      end do
      exit

    end do
    close(10)
  end if

  ! If the alignment will be made relative to average structure

  if( alignto == -1 ) then
    write(*,"('  Computing average structure:   0.00%')",advance='no') 

    ! Allocate auxiliary arrays of these calculations

    allocate( xfirst1(n_reference), yfirst1(n_reference), zfirst1(n_reference),&
              xfirst2(n_rmsd), yfirst2(n_rmsd), zfirst2(n_rmsd) )

    ! The computation of the average structure involves, first, the alignment
    ! of each frame to some frame (here, the first frame). Then, after alignment,
    ! the aligned structure can be used for computing the average, thus it is
    ! summed up in the *ref arrays. 
 
    ifcount = 0
    dcdloop1: do idcd = 1, ndcd

      ! Skip frames is smaller than first frame 

      if ( firstframe > ifcount + nfdcd(idcd) ) then
        ifcount = ifcount + nfdcd(idcd)
        cycle
      end if 

      open(10,file=dcdfile(idcd),action='read',form='unformatted')
      read(10) dummyc, dummyi, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
      read(10) dummyi, dummyr
      read(10) ntotat

      frames: do iframe = 1, nfdcd(idcd)
        ifcount = ifcount + 1
        write(*,"( 7a, f6.2, '%' )",advance='no') (char(8),ic=1,7), 100.*(float(ifcount)/lastframe)

        ! If the last frame to be considered is reached, exit loop

        if ( ifcount > lastframe ) exit 

        ! Read coordinates from dcd file in this frame

        if(dcdaxis) read(10)
        read(10) (xdcd(k), k = 1, lastatom)
        read(10) (ydcd(k), k = 1, lastatom)            
        read(10) (zdcd(k), k = 1, lastatom)

        ! If the first frame was not reached yet or this frame is to be ignored, cycle

        if( mod(ifcount,stride) /= 0 .or. ifcount < firstframe) cycle

        ! On the first frame, just read the coordinates and initialize the
        ! average reference position vectors
        
        if ( ifcount/stride == 1 ) then
          do j = 1, n_reference
            xfirst1(j) = xdcd(i_reference(j))
            yfirst1(j) = ydcd(i_reference(j))
            zfirst1(j) = zdcd(i_reference(j))
            xref1(j) = 0.d0
            yref1(j) = 0.d0
            zref1(j) = 0.d0
          end do
          do j = 1, n_rmsd
            xfirst2(j) = xdcd(i_rmsd(j))
            yfirst2(j) = ydcd(i_rmsd(j))
            zfirst2(j) = zdcd(i_rmsd(j))
            xref2(j) = 0.d0
            yref2(j) = 0.d0
            zref2(j) = 0.d0
            xrmsd_ini(j) = xdcd(i_rmsd(j))
            yrmsd_ini(j) = ydcd(i_rmsd(j))
            zrmsd_ini(j) = zdcd(i_rmsd(j))
          end do
          call compute_cm_dcd(n_reference,i_reference,mass,mass1,&
                              xdcd,ydcd,zdcd,0,&
                              cm1rx,cm1ry,cm1rz)
          call compute_cm_dcd(n_rmsd,i_rmsd,mass,mass2,&
                              xdcd,ydcd,zdcd,0,&
                              cm2rx,cm2ry,cm2rz)
        end if

        ! Align current frame to the first frame: reference atoms

        call compute_cm_dcd(n_reference,i_reference,mass,mass1,&
                            xdcd,ydcd,zdcd,0,&
                            cm1x,cm1y,cm1z)
        call align(n_reference,i_reference,0,cm1x,cm1y,cm1z,xdcd,ydcd,zdcd,&
                   xfirst1,yfirst1,zfirst1,cm1rx,cm1ry,cm1rz,u,&
                   xm,ym,zm,xp,yp,zp) 
        do j = 1, n_reference
          call trans_rot(xdcd,ydcd,zdcd,i_reference(j),cm1x,cm1y,cm1z,u,&
                         cm1rx,cm1ry,cm1rz,xmove,ymove,zmove)
          xref1(j) = xref1(j) + xmove
          yref1(j) = yref1(j) + ymove
          zref1(j) = zref1(j) + zmove
        end do

        ! Align current frame to the first frame: compute_rmsd_of atoms

        call compute_cm_dcd(n_rmsd,i_rmsd,mass,mass2,&
                            xdcd,ydcd,zdcd,0,&
                            cm2x,cm2y,cm2z)
        call align(n_rmsd,i_rmsd,0,cm2x,cm2y,cm2z,xdcd,ydcd,zdcd,&
                   xfirst2,yfirst2,zfirst2,cm2rx,cm2ry,cm2rz,u,&
                   xm,ym,zm,xp,yp,zp) 
        do j = 1, n_rmsd
          call trans_rot(xdcd,ydcd,zdcd,i_rmsd(j),cm2x,cm2y,cm2z,u,&
                         cm2rx,cm2ry,cm2rz,xmove,ymove,zmove)
          xref2(j) = xref2(j) + xmove
          yref2(j) = yref2(j) + ymove
          zref2(j) = zref2(j) + zmove
        end do

      end do frames
      close(10)

    end do dcdloop1
    write(*,*)

    ! Computing average structure

    do j = 1, n_reference
      xref1(j) = float(stride) * xref1(j) / float(lastframe-firstframe+1)
      yref1(j) = float(stride) * yref1(j) / float(lastframe-firstframe+1)
      zref1(j) = float(stride) * zref1(j) / float(lastframe-firstframe+1)
    end do
    do j = 1, n_rmsd
      xref2(j) = float(stride) * xref2(j) / float(lastframe-firstframe+1)
      yref2(j) = float(stride) * yref2(j) / float(lastframe-firstframe+1)
      zref2(j) = float(stride) * zref2(j) / float(lastframe-firstframe+1)
    end do

    deallocate ( xfirst1, yfirst1, zfirst1,&
                 xfirst2, yfirst2, zfirst2 )

  end if

end if

! Computing the position of the center of mass of the references positions

call compute_cm(n_reference,i_reference,mass,mass1,&
                xref1,yref1,zref1,cm1rx,cm1ry,cm1rz)
write(*,"( t2,a,3(f8.3) )") ' CM of reference at reference coordinates: ',&
                            cm1rx, cm1ry, cm1rz

call compute_cm(n_rmsd,i_rmsd,mass,mass2,&
                xref2,yref2,zref2,cm2rx,cm2ry,cm2rz)
write(*,"( t2,a,3(f8.3) )") ' CM of compute_rmsd_of at reference coordinates: ',&
                            cm2rx, cm2ry, cm2rz

! Open optional output pdb files if required

if(writerel) open(30,file=pdboutrel)
if(writeinter) open(31,file=pdboutinter)

! Open output file and writes all information of this run

open(51,file=output(1:length(output)))
write(51,"( '#',/,&
           &'# Output of align.f90:',/,&
           &'# Input file: ',a )") inputfile(1:length(inputfile))
do idcd = 1, ndcd
  record2 = dcdfile(idcd)
  write(51,"( '# DCD file: ', i5,':',a )") idcd, record2(1:length(record2))
end do
write(51,"( '# Group file: ',a,/,&
           &'# PSF file: ',a,/,& 
           &'# First frame: ',i5,' Last frame: ',i5,' Stride: ',i5,/,&
           &'#',/,&
           &'# Number of atoms and mass of group 1: ',i6,f12.3,/,&
           &'# First and last atoms of group 1: ',i6,tr1,i6,/,&
           &'# Number of atoms and mass of group 2: ',i6,f12.3,/,&
           &'# First and last atoms of group 2: ',i6,tr1,i6,/,&
           &'#',/,&
           &'# Mobility of compute_rmsd_of group:',/,&
           &'#  TIME ',t24,'RMSD',t31,'INTERNAL-MOTION RMSD')")& 
           &groupfile(1:length(groupfile)),&
           &psffile(1:length(psffile)),&
           &firstframe, lastframe, stride,&
           &n_reference, mass1, i_reference(1), i_reference(n_reference),& 
           &n_rmsd, mass2, i_rmsd(1), i_rmsd(n_rmsd) 

! Reading the dcd files for doing the final alignment

write(*,"(tr2,52('-'),/ )")
write(*,"( /,'  Computing structural alignment ... ' )")

! Initializing rmsd per residue

do i = 1, number_of_residues
  rmsd_residue(i) = 0.
end do

ifcount = 0
dcdloop2:do idcd = 1, ndcd

  if ( firstframe > ifcount + nfdcd(idcd) ) then
    ifcount = ifcount + nfdcd(idcd)
    cycle
  end if

  record2 = dcdfile(idcd)
  open(10,file=dcdfile(idcd),action='read',form='unformatted')
  read(10) dummyc, dummyi, (dummyi,i=1,8), dummyr, (dummyi,i=1,9)
  read(10) dummyi, dummyr
  read(10) ntotat

  memframes = memory / ntotat
  ncycles = nfdcd(idcd) / memframes + 1
  memlast = nfdcd(idcd) - memframes * ( ncycles - 1 )

  ! Beautiful dcd reading output

  write(*,"( /, tr2,'DCD: ',i5,' Reading: ',f6.2,'% Computing: ',f6.2,'%' )",&
             advance='no') idcd, 0.00, 0.00

  ! Reading dcd file and aligning

  do icycle = 1, ncycles 
  
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
      if ( ifcount + j > lastframe ) exit
      write(*,"( 46a,'DCD: ',i5,' Reading: ',f6.2,'% Computing: ',f6.2,'%' )",&
                advance='no') (char(8),ic=1,46), idcd, 100.*(float(ifcount+j)/lastframe),&
                100.*(float(ifcount)/lastframe)
    end do

    ! Computing the alignment
  
    iatom = 0
    do k = 1, nfrcycle
      ifcount = ifcount + 1

      if ( ifcount > lastframe ) exit dcdloop2
  
      if(mod(ifcount,stride) == 0 .and. ifcount >= firstframe) then
        write(*,"( 7a, f6.2, '%' )",advance='no') (char(8),ic=1,7), 100.*(float(ifcount)/lastframe)

        ! Aligning compute_rmsd_of to itself in the reference and computing the RMSD of
        ! internal motions
  
        call compute_cm_dcd(n_rmsd,i_rmsd,mass,mass2,&
                            xdcd,ydcd,zdcd,iatom,&
                            cm2x,cm2y,cm2z)
        call align(n_rmsd,i_rmsd,iatom,cm2x,cm2y,cm2z,xdcd,ydcd,zdcd,&
                   xref2,yref2,zref2,cm2rx,cm2ry,cm2rz,u,&
                   xm,ym,zm,xp,yp,zp) 
        rmsdint = 0.d0
        if(writeinter) write(31,"( 'REMARK FRAME: ',i10 )") ifcount
        do i = 1, n_rmsd
          ii = iatom + i_rmsd(i)
          call trans_rot(xdcd,ydcd,zdcd,ii,cm2x,cm2y,cm2z,u,cm2rx,cm2ry,cm2rz,&
                         xmove,ymove,zmove)
          rmsdint = rmsdint + ( xmove - xref2(i) )**2 &
                            + ( ymove - yref2(i) )**2 &
                            + ( zmove - zref2(i) )**2
          if(writeinter)& 
          write(31,"( 'ATOM',t7,i5,t13,a4,t18,a3,t23,i4,t31,3(f8.3),t73,a3 )")& 
                i_rmsd(i),typeat(i_rmsd(i)),resat(i_rmsd(i)),resid(i_rmsd(i)),&
                xmove,ymove,zmove,segat(i_rmsd(i))
  
        end do 
        rmsdint = sqrt(rmsdint/float(n_rmsd))
        if(writeinter) write(31,"( 'END' )")
  
        ! Aligning 'reference' to its reference coordinates to obtain the 
        ! rotation matrix, and then move rmsd accordingly, to get the relative 
        ! movement of the groups
  
        call compute_cm_dcd(n_reference,i_reference,mass,mass1,&
                            xdcd,ydcd,zdcd,iatom,&
                            cm1x,cm1y,cm1z)
        call align(n_reference,i_reference,iatom,cm1x,cm1y,cm1z,xdcd,ydcd,zdcd,&
                   xref1,yref1,zref1,cm1rx,cm1ry,cm1rz,u,&
                   xm,ym,zm,xp,yp,zp)
        rmsd = 0.d0
        do i = 1, number_of_residues
          rmsd_residue_frame(i) = 0.
        end do
        if(writerel) write(30,"( 'REMARK FRAME: ',i10 )") ifcount
        do i = 1, n_rmsd
          ii = iatom + i_rmsd(i)
          call trans_rot(xdcd,ydcd,zdcd,ii,cm1x,cm1y,cm1z,u,cm1rx,cm1ry,cm1rz,&
                         xmove,ymove,zmove)
          rmsdnow =  ( xmove - xref2(i) )**2 &
                   + ( ymove - yref2(i) )**2 &
                   + ( zmove - zref2(i) )**2
          rmsd = rmsd + rmsdnow
          rmsd_residue_frame(i_residue(i)) = rmsd_residue_frame(i_residue(i)) + rmsdnow
  
          if(writerel)& 
          write(30,"( 'ATOM',t7,i5,t13,a4,t18,a3,t23,i4,t31,3(f8.3),t73,a3 )")& 
                i_rmsd(i),typeat(i_rmsd(i)),resat(i_rmsd(i)),resid(i_rmsd(i)),&
                xmove,ymove,zmove,segat(i_rmsd(i))
          
        end do 
        do i = 1, number_of_residues
          rmsd_residue(i) = rmsd_residue(i) + sqrt(rmsd_residue_frame(i)/natoms_per_residue(i))
        end do
        rmsd = sqrt(rmsd/float(n_rmsd))
  
        if(writerel) then
          do i = 1, n_reference
            ii = iatom + i_reference(i)
            call trans_rot(xdcd,ydcd,zdcd,ii,cm1x,cm1y,cm1z,u,cm1rx,cm1ry,cm1rz,&
                           xmove,ymove,zmove)
            write(30,"( 'ATOM',t7,i5,t13,a4,t18,a3,t23,i4,t31,3(f8.3),t73,a3 )")& 
                  i_reference(i),typeat(i_reference(i)),resat(i_reference(i)),&
                  resid(i_reference(i)),&
                  xmove,ymove,zmove,segat(i_reference(i))
          end do
          write(30,"( 'END' )")
        end if 
  
        ! Printing to output file in this frame 
   
        write(51,"( e12.6,tr2,f17.6,tr6,f17.6 )") scaletime*ifcount, rmsd, rmsdint
      end if

      iatom = iatom + ntotat
    end do
  end do
  close(10)
end do dcdloop2
close(51)

! Printint the rmsd per residue

if(write_per_residue) then
  open(32,file=write_per_residue_file)
  write(32,"( 'REMARK: Occupancy column contains average RMSD per residue.' )")
  write(32,"( 'REMARK: B-factor column contains relative average RMSD per residue.' )")
  max_rmsd = 0.d0
  min_rmsd = 1.d10
  do i = 1, number_of_residues
    rmsd_residue(i) = float(stride)*rmsd_residue(i)/float(lastframe-firstframe+1)
    max_rmsd = amax1(max_rmsd,rmsd_residue(i))
    min_rmsd = amin1(min_rmsd,rmsd_residue(i))
  end do
  do i = 1, n_rmsd
    if ( (max_rmsd-min_rmsd) > 1.d-7 ) then 
      relative_rmsd = (rmsd_residue(i_residue(i))-min_rmsd)/(max_rmsd-min_rmsd)
    else
      relative_rmsd = 1.
    end if
    write(32,"( 'ATOM',t7,i5,t13,a4,t18,a3,t23,i4,t31,3(f8.3),t55,f5.2,t61,f5.2,t73,a3)")&
          i_rmsd(i),typeat(i_rmsd(i)),resat(i_rmsd(i)),resid(i_rmsd(i)),&
          xrmsd_ini(i),yrmsd_ini(i),zrmsd_ini(i),&
          rmsd_residue(i_residue(i)),relative_rmsd,&
          segat(i_rmsd(i))
  end do
  close(32)
end if

if(writerel) close(30)
if(writeinter) close(31)

! Write final messages with names of output files and their content

time0 = etime(tarray) - time0
write(*,*)
write(*,"( tr2,52('-') )")
write(*,*)
write(*,*) ' OUTPUT FILES: ' 
write(*,*)
write(*,*) ' Wrote output file: ', output(1:length(output))
write(*,*)
write(*,*) ' Which contains the RMSD of the i_rmsd(compute_rmsd_of) group'
write(*,*) ' relative to i_reference(reference) and its internal mobility' 
write(*,*) ' as function of simulation frames. '
write(*,*) 
if(writerel) then
  write(*,*)  ' and file: ', pdboutrel(1:length(pdboutrel))
  write(*,*)  ' containing the groups aligned to the reference.'
  write(*,*)
end if
if(writeinter) then
  write(*,*)  ' and file: ', pdboutinter(1:length(pdboutinter))
  write(*,*)  ' containing internal motions of the reference.'
  write(*,*)
end if
if(write_per_residue) then
  write(*,*)  ' and file: ', write_per_residue_file(1:length(write_per_residue_file)) 
  write(*,*)  ' containing average RMSDs per residue.'
  write(*,*)
end if
write(*,*) ' Running time: ', time0
write(*,*) '####################################################'
write(*,*) 
write(*,*) '  END: Normal termination.  '
write(*,*) 
write(*,*) '####################################################'
write(*,*)        

end

!
! Subroutine trans_rot: Translates the coordinates of according to the 
!                       center of mass of the aligned group (as moving
!                       this group to the origin), then rotates according
!                       to the rotation matrix, then translates the atom
!                       according to the center of mass of the reference
!                       used for the alignment.
! 
! On input:
!   xdcd, ydcd, zdcd: Vector containing atoms' coordinates.
!   index: Index of the atom in the vectors of coordinates.
!   cmx, cmy, cmz: Center of mass of the aligned group.
!   u: rotation matrix (3x3)
!   cmrx, cmry, cmrz: Center of mass of the reference group.
!
! On output:
!   xmove, ymove, zmove: Coordinates rotated and translated.
!
!                        


subroutine trans_rot(xdcd,ydcd,zdcd,index,cmx,cmy,cmz,u,cmrx,cmry,cmrz,&
                     xmove,ymove,zmove)

implicit none
integer :: index
real :: xdcd(*), ydcd(*), zdcd(*)
real :: cmx, cmy, cmz, u(3,3), cmrx, cmry, cmrz,&
        xmove, ymove, zmove

xmove = ( u(1,1)*(xdcd(index)-cmx) +&
          u(1,2)*(ydcd(index)-cmy) +&
          u(1,3)*(zdcd(index)-cmz) ) + cmrx
ymove = ( u(2,1)*(xdcd(index)-cmx) +&
          u(2,2)*(ydcd(index)-cmy) +&
          u(2,3)*(zdcd(index)-cmz) ) + cmry
zmove = ( u(3,1)*(xdcd(index)-cmx) +& 
          u(3,2)*(ydcd(index)-cmy) +&
          u(3,3)*(zdcd(index)-cmz) ) + cmrz

return
end
             
!
! Subroutine align: given two sets of vectors, finds the best
!                   rotation that to align the two sets. Returns
!                   the second vector aligned to the first vector
!
!                 Method: S. K. Kearsley, 
!                         "On the orthogonal transformation used for
!                          structural comparisons"
!                         Acta Cryst. (1989) A45, 208-210 
!
!  Author: Leandro Martinez, IQ-UNICAMP, 26/10/2005
!
!  On input: n_align: the number of atoms of the group to be aligned
!            i_align: the indexes in vectors x-, y-, z- and mass of the atoms
!                     of the atoms to be aligned (1,...,n_align)
!            iatom: index of atom 0 (first atom-1) for xdcd, ydcd and zdcd
!            cmx,cmy,cmz: Center of mass of the atoms to be aligned
!            xdcd: The vector containing current x coordinates
!            ydcd: The vector containing current y coordinates
!            zdcd: The vector containing current z coordinates
!            xref: The vector containing referece x coordinates
!            yref: The vector containing referece y coordinates
!            zref: The vector containing referece z coordinates
!            cmxr, cmyr, cmzr: Center of mass of reference atoms.
!            Auxiliar arrays: xm, ym, zm, xp, yp, zp are vectors
!                             that have n_align positions. 
!
!  On return: d(3): movement of the center of mass to be applied 
!             u(3,3): rotation matrix to be applied
!
!  Obs: It is assumed that the reference atoms are already moved
!       in such a way that their center of mass is in the origin.
!

subroutine align(n_align,i_align,iatom,cmx,cmy,cmz,xdcd,ydcd,zdcd,&
                 xref,yref,zref,cmxr,cmyr,cmzr,u,xm,ym,zm,xp,yp,zp)

  implicit none
  integer :: i, j, iq, n_align, i_align(*), iatom, ii
  real :: xdcd(*), ydcd(*), zdcd(*), xref(*), yref(*), zref(*)
  real :: cmx, cmy, cmz, u(3,3), a(4,4), q(4,4), qmin,&
          xm(*), ym(*), zm(*), xp(*), yp(*), zp(*),&
          cmxr, cmyr, cmzr
  
  ! For ssyev
  real :: work(12)
  integer :: info
  
  ! If the number of atoms of this group is 1, return the identity matrix
  
  if( n_align == 1 ) then
    do i = 1, 3
      do j = 1, 3
        if(i == j) then
          u(i,j) = 1.d0
        else
          u(i,j) = 0.d0
        end if
      end do
    end do
    return
  end if
   
  ! Computing the quaternion matrix
  
  do i = 1, n_align
    ii = iatom + i_align(i)
    xm(i) = ( xref(i) - cmxr ) - ( xdcd(ii) - cmx )
    ym(i) = ( yref(i) - cmyr ) - ( ydcd(ii) - cmy )
    zm(i) = ( zref(i) - cmzr ) - ( zdcd(ii) - cmz )
    xp(i) = ( xref(i) - cmxr ) + ( xdcd(ii) - cmx )
    yp(i) = ( yref(i) - cmyr ) + ( ydcd(ii) - cmy )
    zp(i) = ( zref(i) - cmzr ) + ( zdcd(ii) - cmz )
  end do
  
  do i = 1, 4
    do j = 1, 4
      q(i,j) = 0.d0
    end do
  end do
   
  do i = 1, n_align
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)  
     
  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix

  call ssyev('V','U',4,q,4,a,work,12,info)
   
  ! Computing the rotation matrix
  
  iq = 1
  u(1,1) = q(1,iq)**2 + q(2,iq)**2 - q(3,iq)**2 - q(4,iq)**2
  u(1,2) = 2. * ( q(2,iq)*q(3,iq) + q(1,iq)*q(4,iq) )
  u(1,3) = 2. * ( q(2,iq)*q(4,iq) - q(1,iq)*q(3,iq) )  
  u(2,1) = 2. * ( q(2,iq)*q(3,iq) - q(1,iq)*q(4,iq) )  
  u(2,2) = q(1,iq)**2 + q(3,iq)**2 - q(2,iq)**2 - q(4,iq)**2 
  u(2,3) = 2. * ( q(3,iq)*q(4,iq) + q(1,iq)*q(2,iq) )  
  u(3,1) = 2. * ( q(2,iq)*q(4,iq) + q(1,iq)*q(3,iq) )  
  u(3,2) = 2. * ( q(3,iq)*q(4,iq) - q(1,iq)*q(2,iq) )  
  u(3,3) = q(1,iq)**2 + q(4,iq)**2 - q(2,iq)**2 - q(3,iq)**2 
  
return
end subroutine align



 
