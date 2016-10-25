!
! perresidueavg: A program that reads the output of solvation when it prints
!                the solvation per residue and outputs the average solvation
!                per residue.
!
!
! L. Martinez, April 17, 2015. 
! Institute of Chemistry - University of Campinas
!

program perresidueavg

  implicit none
  integer :: i, ncol, ioerr, nres, ndata
  real :: t
  real, allocatable :: solvres(:), add(:)
  character(len=200) :: file
  character(len=50000) :: line

  call getarg(1,file)
 
  open(10,file=file,status='old',action='read')

  ! Get the number of residues

  do
    read(10,"( a50000 )",iostat=ioerr) line
    if ( ioerr /= 0 ) exit
    if ( line(1:1) == "#" ) cycle
    ncol = 1
    do
      read(line,*,iostat=ioerr) (t,i=1,ncol)
      if ( ioerr /= 0 ) exit
      ncol = ncol + 1
    end do
    exit
  end do
  nres = ncol - 2
  write(*,"( a, a )") '# File = ', trim(file)
  write(*,"( a, i10 )") '# Number of residues = ', nres

  allocate(solvres(nres),add(nres))

  do i = 1, nres
    solvres(i) = 0.
  end do

  rewind(10)
  ndata = 0
  lines: do
    read(10,"( a50000 )",iostat=ioerr) line
    if ( ioerr /= 0 ) exit
    if ( line(1:1) == "#" ) cycle
    backspace(10)
    do
      read(10,"( a50000 )",iostat=ioerr) line
      if ( ioerr /= 0 ) exit lines
      read(line,*,iostat=ioerr) t, (add(i),i=1,nres)
      if ( ioerr /= 0 ) exit lines
      ndata = ndata + 1
      do i = 1, nres
        solvres(i) = solvres(i) + add(i)
      end do
    end do
  end do lines
  close(10)

  write(*,"( a, i10 )") '# Number of data points = ', ndata

  do i = 1, nres
    solvres(i) = solvres(i) / ndata
  end do

  write(*,"( a )") '# RESIDUE_INDEX   AVERAGE_SOLVATION'
  do i = 1, nres
    write(*,*) i, solvres(i)
  end do

end program perresidueavg

