!
! Program data_average: Computes the average curve of several
! data sets, in interpolated points
!
! L. Martinez
! Institute of Chemistry, University of Campinas
! Dec 2, 2016
! http://leandro.iqm.unicamp.br
!

program data_average

  use file_operations
  implicit none

  type data_type
    character(len=200) :: file
    integer :: n
    double precision :: xmin, xmax
    double precision, allocatable :: x(:), y(:), yint(:)
  end type data_type

  integer :: narg, iargc, ioerr
  integer :: i, j, k, nomit
  integer :: ndata, xcol, ycol, nint, nj
  double precision :: xmin, xmax, step, xread, yread, xlast, xint, ytemp, median
  double precision :: interpolate, sd
  double precision, allocatable :: average(:), ysort(:)
  character(len=200) :: record, record2, lineformat
  type(data_type), allocatable :: data(:)

  ! Maximum number of data sets for which the output is extensive:
  nomit = 20

  narg = iargc()
  if ( narg < 4 ) then
    write(*,*) ' Run with: data_average xcol ycol data1.dat data2.dat .... ' 
    stop
  end if
 
  ndata = narg - 2
  allocate( data(ndata) )
  
  call getarg(1,record)
  read(record,*,iostat=ioerr) xcol 
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: First argument must be xcol '
    stop
  end if
  call getarg(2,record)
  read(record,*,iostat=ioerr) ycol 
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: First argument must be ycol '
    stop
  end if

  ! Read properties of data sets

  j = 2
  xmax = -1.d30
  xmin = 1.d30
  step = 1.d30
  do i = 1, ndata

    data(i)%n = 0
    data(i)%xmin = 1.d30
    data(i)%xmax = -1.d30

    j = j + 1
    call getarg(j,data(i)%file)
    open(10,file=data(i)%file,iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open file: ', trim(adjustl(data(i)%file))
      stop
    end if

    ! Read data from file

    xlast = -1.d30
    do

      read(10,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      if ( comment(record) ) cycle

      ! Read x value from this line

      read(record,*,iostat=ioerr) (record2, k = 1, xcol)
      if ( ioerr /= 0 ) cycle
      read(record2,*,iostat=ioerr) xread
      if ( ioerr /= 0 ) cycle

      ! Read y value
      
      read(record,*,iostat=ioerr) (record2, k = 1, ycol)
      if ( ioerr /= 0 ) cycle
      read(record2,*,iostat=ioerr) yread
      if ( ioerr /= 0 ) cycle

      data(i)%n = data(i)%n + 1
      data(i)%xmin = min(data(i)%xmin,xread)
      data(i)%xmax = max(data(i)%xmax,xread)

      xmin = min(xmin,xread)
      xmax = max(xmax,xread)
      if ( xlast < -1.d20 ) then
        xlast = xread
      else
        step = min(step,dabs(xread-xlast))
        xlast = xread
      end if
      if ( step == 0.d0 ) then
        write(*,*) ' ERROR: Two consecutive x coordinates are identical. '
        write(*,*) '  File: ', trim(adjustl(data(i)%file))
        stop
      end if

    end do
    close(10)
    if ( ndata <= nomit ) then
      lineformat = "(a,i8,tr2,a,a,i8,a,f17.5,a,f17.5,a)"
      if ( abs(data(i)%xmin) > 1d5 .or. abs(data(i)%xmax) > 1d5 ) then
        lineformat = "(a,i8,tr2,a,a,i8,a,e17.5,a,e17.5,a)"
      end if
      write(*,lineformat) &
              "# Data set ", i, trim(adjustl(data(i)%file)), &
              " N = ", data(i)%n, &
              ", Range = [", data(i)%xmin, ",", data(i)%xmax, "]"
    end if
    if ( ndata > nomit ) then
      if ( i <= 5 .or. i > ndata - 5 ) then
        write(*,"(a,a,a,i8)") "# Number of data points in file: ", &
                              trim(adjustl(data(i)%file)), ": ", data(i)%n
      end if
      if ( i == 5 ) then
        write(*,"('# ... (omitting details)')")
      end if
    end if

    allocate(data(i)%x(data(i)%n))
    allocate(data(i)%y(data(i)%n))

  end do

  write(*,"(a,i8)") "# Number of data sets: ", ndata
  write(*,"(a,f12.5)") "# Minimum x value: ", xmin
  write(*,"(a,f12.5)") "# Maximum x value: ", xmax
  write(*,"(a,f12.5)") "# Minimum step in x: ", step

  ! Number of interpolated points

  nint = int((xmax - xmin) / step) + 1
  write(*,"(a,i8)") "# Number of x points in interpolated grid: ", nint

  ! Now, actually reading the data sets

  do i = 1, ndata

    open(10,file=data(i)%file,iostat=ioerr)
    j = 0
    do

      read(10,"(a200)",iostat=ioerr) record
      if ( ioerr /= 0 ) exit
      if ( comment(record) ) cycle

      ! Read x value from this line

      read(record,*,iostat=ioerr) (record2, k = 1, xcol)
      if ( ioerr /= 0 ) cycle
      read(record2,*,iostat=ioerr) xread
      if ( ioerr /= 0 ) cycle

      ! Read y value
      
      read(record,*,iostat=ioerr) (record2, k = 1, ycol)
      if ( ioerr /= 0 ) cycle
      read(record2,*,iostat=ioerr) yread
      if ( ioerr /= 0 ) cycle

      j = j + 1
      data(i)%x(j) = xread
      data(i)%y(j) = yread

    end do
    close(10)

  end do

  ! Interpolating data points

  allocate(average(nint))
  do i = 1, ndata
    allocate(data(i)%yint(nint))
  end do
  do i = 1, nint
    xint = xmin + step*(i-1)
    average(i) = 0.d0
    nj = 0
    do j = 1, ndata
      if ( xint >= data(j)%xmin .and. xint <= data(j)%xmax ) then
        data(j)%yint(i) = interpolate(data(j)%n,data(j)%x,data(j)%y,xint) 
        average(i) = average(i) + data(j)%yint(i)
        nj = nj + 1
      end if
    end do
    average(i) = average(i) / nj
  end do
  
  ! Computing standard devation, and printing output

  allocate(ysort(nint))
  write(*,"('#            X     Average Y            SD    SD/sqrt(N)        Median')")
  do i = 1, nint
    xint = xmin + step*(i-1)
    sd = 0.d0
    nj = 0

    ! Compute median

    do j = 1, ndata
      ysort(j) = data(j)%yint(i)
    end do
    do k = 1, ndata-1
      j = k + 1
      do while( ysort(j-1) < ysort(j) )
        ytemp = ysort(j-1)
        ysort(j-1) = ysort(j) 
        ysort(j) = ytemp
        j = j - 1
        if ( j == 1 ) exit
      end do
    end do
    if ( mod(ndata,2) == 0 ) then
      median = 0.5d0*(ysort(ndata/2) + ysort(ndata/2+1))
    else
      median = ysort(ndata/2+1)
    end if
    
    do j = 1, ndata
      if ( xint >= data(j)%xmin .and. xint <= data(j)%xmax ) then
        sd = sd + (data(j)%yint(i)-average(i))**2
        nj = nj + 1
      end if
    end do
    sd = dsqrt(sd/(nj-1))
    write(*,"(5(tr2,f12.5))") xint, average(i), sd, sd/dsqrt(dble(nj)), median
  end do
  
end program data_average

double precision function interpolate(n,x,y,xint)

  integer :: i
  double precision :: x(n), y(n), xint

  i = 1
  do while( x(i) < xint ) 
    i = i + 1
    if ( i == n ) then
      interpolate = y(n)
      return
    end if
  end do
  if ( i == 1 ) then
    interpolate = y(1)
    return
  end if
  interpolate = y(i-1)+((y(i)-y(i-1))/(x(i)-x(i-1)))*(xint-x(i-1))

end function interpolate





















