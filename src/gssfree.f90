!
! This program reads the output from gss.f90 and computes the 
! free energy profile (-RT ln gss) and the integrated solvation
! free energy ( \int_0^\infty gss(r)*n(r) dr )
!
! The user must provide the output data file of gss.f90 and the 
! temperature. The output is provided in kcal/mol.
!
! Optionally, the user can provide a precision, which defines from
! which distance the gss is converged. The default is to consider
! the gss converged when it reached the final value with a precision
! of 1%. Then, the average of this range will be taken as the limiting
! value of the gss.
!
! L. Martinez
! Institute of Chemistry, University of Campinas
! Nov 10, 2015
! http://leandro.iqm.unicamp.br
!

program gssfree

  implicit none
  integer :: iargc, ioerr, i, ifinal, j, ndata
  double precision :: d, gss, dummyd, rtlngss, gsolv, &
                      temp, gsslimit, average, precision, &      
                      dfinal, cutoff
  double precision, allocatable :: dlist(:), gsslist(:), nmols(:)
  character(len=200) :: file, record

  ! Gas constant, kcal/mol

  double precision, parameter :: R = 1.9872035d-3

  ! Reading input paramters

  if ( iargc() /= 2 .and. iargc() /= 3 .and. iargc() /= 4 ) then
    write(*,*) ' Run with: gssfree input.dat 298.15 [cutoff] [precision] > output.dat '
    stop
  end if

  call version()

  call getarg(1,file)
  call getarg(2,record)
  read(record,*,iostat=ioerr) temp
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not read temperature. '
    write(*,*) ' Run with: gssfree input.dat 298.15 [cutoff] [precision] > output.dat '
    stop
  end if
  
  ! Define distance cutoff 
  cutoff = -1.d0
  if ( iargc() > 2 ) then
    call getarg(3,record)
    read(record,*,iostat=ioerr) cutoff
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read cutoff. '
      write(*,*) ' Run with: gssfree input.dat 298.15 [cutoff] [precision] > output.dat '
      stop
    end if
  end if
  if ( cutoff < 0.d0 ) then
    write(*,"(a)") "# GSS cutoff: end of data file. "
  else
    write(*,"(a,f12.5)") "# GSS cutoff: ", cutoff
  end if

  ! Define precision of gss convergence to unity (Default: 1%)
  precision = 0.01
  if ( iargc() == 4 ) then
    call getarg(4,record)
    read(record,*,iostat=ioerr) precision
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not read precision. '
      write(*,*) ' Run with: gssfree input.dat 298.15 [cutoff] [precision] > output.dat '
      stop
    end if
  end if
  write(*,"(a,f12.5)") "# Expected GSS convergence precision: ", precision

  open(10,file=file,action="read",iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open data file: ', trim(file)
    stop
  end if

  ! Printint output

  write(*,"(a)") "#"
  write(*,"(a)") "# Output of gssfree "
  write(*,"(a)") "#"
  write(*,"(a,a)") "# GSS output file: ", trim(adjustl(file))
  write(*,"(a,f8.3)") "# Temperature: ", temp
  write(*,"(a,f8.3)") "# Cutoff: ", cutoff
  write(*,"(a)") "#"

  ! First, checking the limit of the distribution function (should
  ! be as close to unity as possible)

  ndata = 0
  average = 0.d0
  do
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) d, (dummyd,i=1,6)
    if ( ioerr /= 0 ) cycle
    if ( cutoff < 0.d0 .or. d <= cutoff ) then 
      ndata = ndata + 1
      read(record,*) d, gss
    end if
  end do
  gsslimit = gss
  rewind(10)
  allocate( dlist(ndata), gsslist(ndata), nmols(ndata) )
  
  i = 0
  do
    read(10,"( a200 )",iostat=ioerr) record
    if ( ioerr /= 0 ) exit
    read(record,*,iostat=ioerr) d, (dummyd,j=1,6)
    if ( ioerr /= 0 ) cycle
    if ( cutoff < 0.d0 .or. d <= cutoff ) then
      i = i + 1
      read(record,*) dlist(i), gsslist(i), nmols(i)
      if ( abs(gsslist(i)-gsslimit)/gsslimit > precision*gsslimit ) then
        ifinal = i
        dfinal = dlist(i)
      end if
    end if
  end do
  write(*,"(a,f12.5)") "# Will consider the gss converged at distance: ", dfinal
  close(10)

  ! Computing the average of the converged region

  j = 0
  average = 0.d0
  do i = ifinal, ndata
    average = average + gsslist(i)
  end do
  average = average / ( ndata - ifinal + 1 )
  write(*,"(a,f12.5)") "# Average of the converged limit: ", average
  average = average - 1.d0

  ! Computing the free-energy parameters 

  write(*,"(a)") "#" 
  write(*,"(a)") "#    Distance    -RT ln(gss)     DG_solv"
  gsolv = 0.d0
  do i = 1, ndata
    if ( gsslist(i) > 1.d-10 ) then
      gsslist(i) = gsslist(i) - average
      if ( gsslist(i) > 0.d0 ) then
        rtlngss = -1.d0*R*temp*dlog(gsslist(i))
        gsolv = gsolv + rtlngss*nmols(i)
        write(*,"( 3(tr2,e12.5) )") dlist(i), rtlngss, gsolv
      end if
    end if
  end do

end program gssfree
