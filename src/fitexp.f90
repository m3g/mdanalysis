!
! Version 16.323
!
!
! Module containing the problem data
!

module problem_data
  implicit none
  integer :: nterms, ndata, n
  double precision :: startat, stopat
  double precision, allocatable :: xdata(:), ydata(:)
end module problem_data

!
! fitexp: Fit multiple exponentials to data
!

program fitexp

use problem_data
implicit none
integer :: i, j, icol, status, ntrials, it, narg, nrepeat, maxtrial, xcol, ycol
double precision :: x_temp, y_temp, f, drand, seed, fbest, f_val,&
                    y0, c, c_precision, r_precision
character(len=1) :: acol
character(len=200) :: inputfile, record, datafile, outputfile, keyword,&
                      value
double precision, allocatable :: x(:), g(:), l(:), u(:), xbest(:)
logical set_bounds

call version

! Some default parameters

ntrials = 5
maxtrial = 1000
startat = 1
stopat = 100000
nterms = 1
outputfile = "fitexp.dat"
set_bounds = .false.
c_precision = 1.d-1
r_precision = 1.d-10
xcol = 1
ycol = 2

! Check arguments

narg = iargc()
if ( narg < 1 .or. narg > 3 ) then
  write(*,*) ' Run with: fitexp input.inp [data.dat] [output.dat] '
  stop
end if

! Reading input file

call getarg(1,inputfile)
open(10,file=inputfile,action="read")
do 
  read(10,"( a200 )",iostat=status) record
  if( status /= 0 ) exit
  if(record(1:1) == "#" .or. len_trim(record) < 1 ) cycle
  select case ( keyword(record) )
    case ("data")
      datafile = value(record)
    case ("output")
      outputfile = value(record)
    case ("startat")
      record = value(record)
      read(record,*) startat
    case ("stopat")
      record = value(record)
      read(record,*) stopat
    case ("n_exp_terms") 
      record = value(record)
      read(record,*) nterms
    case ("n_trials")
      record = value(record)
      read(record,*) ntrials
    case ("max_trials")
      record = value(record)
      read(record,*) maxtrial
    case ("xcol")
      record = value(record)
      read(record,*) xcol
    case ("ycol")
      record = value(record)
      read(record,*) ycol
    case ("linear_term_precision")
      record = value(record)
      read(record,*) c_precision
    case ("exponent_lower_bound",&
          "exponent_upper_bound",&
          "linear_lower_bound",&
          "linear_upper_bound")
      set_bounds = .true.
    case default
      write(*,*) ' ERROR: Unrecognized keyword: ', keyword(record)
      stop
  end select
end do

! Read data file from the command line

if ( narg > 1 ) then
  call getarg(2,record)
  datafile = record
end if

! Read output file from command line

if ( narg > 2 ) then
  call getarg(3,record)
  outputfile = record
end if

! Number of variables

n = nterms * 2

! Allocate x, g, l and u

allocate( x(n), g(n), l(n), u(n), xbest(n) )

! Read upper and lower bounds for exponents (first set default values)

do i = 1, nterms
  l(i) = 0.d0
  u(i) = 1.d0
  l(nterms+i) = 1.d-5
  u(nterms+i) = 1.d4
end do
if( set_bounds ) then
  rewind(10)
  do 
    read(10,"( a200 )",iostat=status) record
    if( status /= 0 ) exit
    select case ( keyword(record) )
      case ("linear_lower_bound")
        read(record(21:200),*) i, x_temp
        l(i) = x_temp 
      case ("linear_upper_bound")
        read(record(21:200),*) i, x_temp
        u(i) = x_temp 
      case ("exponent_lower_bound")
        read(record(21:200),*) i, x_temp
        l(nterms+i) = x_temp 
      case ("exponent_upper_bound")
        read(record(21:200),*) i, x_temp
        u(nterms+i) = x_temp 
    end select
  end do
end if
close(10)

write(*,"( '# fitexp output ',/,&
          &'# Data file: ', a,/,&
          &'# Start at: ',d14.6,/,& 
          &'# Stop at: ',d14.6,/,& 
          &'# Number of terms: ',i8 )")& 
          datafile(1:len_trim(datafile)),&
          startat,stopat
do i = 1, nterms
  write(*,"( '# Lower and upper bound for exponent ',i4,':',&
            &f14.6,tr2,f14.6 )") i, l(nterms+i), u(nterms+i)
end do

! Reading the number or valid data points

open(10,file=datafile,action="read")
ndata = 0
do
  read(10,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(record(1:1) == "#") cycle
  read(record,*,iostat=status) (acol,icol=1,xcol-1), x_temp
  if(status /= 0) cycle 
  read(record,*,iostat=status) (acol,icol=1,ycol-1), y_temp
  if(status /= 0) cycle 
  if( x_temp >= startat .and. &
      x_temp <= stopat ) then 
    ndata = ndata + 1
    if ( ndata == 1 ) y0 = y_temp
  end if
end do
write(*,"(a,i4)") "# Read X values from column: ", xcol
write(*,"(a,i4)") "# Read Y values from column: ", ycol
write(*,"(a,i10)") "# Number of data points: ", ndata
if ( ndata < 3 ) then
  write(*,*) ' ERROR: Number of data points read = ', ndata
  close(10)
  stop
end if
rewind(10)

! Allocate arrays and read data into them

allocate( xdata(ndata), ydata(ndata) )
ndata = 0
do
  read(10,"( a200 )",iostat=status) record
  if(status /= 0) exit
  if(record(1:1) == "#") cycle
  read(record,*,iostat=status) (acol,icol=1,xcol-1), x_temp
  if(status /= 0) cycle 
  read(record,*,iostat=status) (acol,icol=1,ycol-1), y_temp
  if(status /= 0) cycle 
  if( x_temp >= startat .and. &
      x_temp <= stopat ) then 
    ndata = ndata + 1
    xdata(ndata) = x_temp
    ydata(ndata) = y_temp / y0
  end if
end do
write(*,"( '# Number of data points: ',i8 )") ndata
close(10)
write(*,"( '#',/,'#',55('-'),/,35('-') )",advance="no")

! Running optimization

seed = 1.825242d0
fbest = 1.d20
nrepeat = 0
it = 0 
do while( nrepeat < ntrials .and. it <= maxtrial )
  it = it + 1

  ! Random initial point of this trial

  do i = 1, n
    x(i) = l(i) + drand(seed) * ( u(i) - l(i) )
  end do

  ! Call solver

  call call_algencan(n,x,l,u,f)

  ! Checking if this point is feasible

  c = -1.d0
  do i = 1, nterms
    c = c + x(i)
  end do
  if ( dabs(c) > c_precision ) then
    write(*,"( 35a )",advance="no") (char(8),i=1,35)
    write(*,"( '# TRIAL ',i6,11(' '),'UNFEASIBLE' )",advance="no") it
    cycle
  end if

  ! Function value of this trial

  call func(x,f)
  write(*,"( 35a )",advance="no") (char(8),i=1,35)
  write(*,"( '# TRIAL ', i6, ' ERROR = ', e12.6 )",advance="no") it, f

  ! If the solution is the same as before, increase nrepeat

  if ( abs(f-fbest)/n < r_precision ) then
    nrepeat = nrepeat + 1
  else
    if ( f < fbest ) then
      nrepeat = 0
    end if
  end if

  ! Save best solution

  if( f < fbest ) then

    fbest = f
    do i = 1, n
      xbest(i) = x(i)
    end do

    ! Write best point fit to the screen

    write(*,*) ' <== BEST UP TO NOW: '
    write(*,"( '#  y = ',f14.6,'*exp( -x / ',f14.6,')' )")&
                   xbest(1), 1.d0 / xbest(nterms+1)
    do i = 2, nterms
      write(*,"( '#',tr5'+',f14.6,'*exp( -x / ',f14.6,')' )")&
                xbest(i), 1.d0 / xbest(nterms+i)
    end do
    write(*,"( '#',/,'#',55('-'),/,'#',34(' ') )",advance="no")

    fbest = f
    do i = 1, n
      xbest(i) = x(i)
    end do

    ! Writting output file with best point up to now

    open(10,file=outputfile)
    write(10,"( '# fitexp output ',/,&
              &'# Data file: ', a,/,&
              &'# Start at: ',d14.6,/,& 
              &'# Stop at: ',d14.6,/,& 
              &'# Number of terms: ',i8 )")& 
              datafile(1:len_trim(datafile)),&
              startat,stopat
    do i = 1, nterms
      write(10,"( '# Lower and upper bound for ',i6,'th exponent: ',&
                &e14.6,tr2,e14.6 )") i, l(nterms+i), u(nterms+i)
    end do
    write(10,"( '# Best fit found: ' )")
    write(10,"( '# Trial = ',i6 )"), it
    write(10,"( '# Error = ',e14.6 )"), fbest
    write(10,"( '# ',tr5,' y = ' )")
    do i = 1, nterms
      write(10,&
             "( '# ',tr5' +',f14.6,'*exp( -x / ',f14.6,')' )")&
                 xbest(i), 1.d0 / xbest(nterms+i)
    end do
    write(10,"( '#',t11,'X',t35,'DATA',t62,'FIT' )")
    do i = 1, ndata 
      f_val = 0.d0
      do j = 1, nterms
        f_val = f_val + xbest(j) * exp( -xdata(i) * xbest(nterms+j) )   
      end do
      write(10,*) xdata(i), ydata(i), f_val 
    end do
    close(10)
  end if

end do
write(*,*)
write(*,"(a,i5,a)") "# The best solution was found ",ntrials," times."
write(*,"(a)") '#'
write(*,"(a)") '# Decay rates, ordered from higher to lower: '
do i = nterms + 1, 2*nterms-1
  j = i + 1
  do while( xbest(j-1) > xbest(j) )
    ! Rates (used for sorting)
    f = xbest(j-1)
    xbest(j-1) = xbest(j)
    xbest(j) = f
    ! fractions
    f = xbest(j-1-nterms)
    xbest(j-1-nterms) = xbest(j-nterms)
    xbest(j-nterms) = f
    j = j - 1
    if ( j == nterms+1 ) exit
  end do
end do 
write(*,"(a,3(tr2,f14.6))") '#   RATES: ', (1.d0/xbest(i),i=nterms+1,2*nterms)
write(*,"(a,3(tr2,f14.6))") '# WEIGHTS: ', (xbest(i),i=1,nterms)
write(*,"(a)") '#'
write(*,"(a,a)") '# Wrote output file: ', trim(adjustl(outputfile))
write(*,"(a)") '#'
write(*,"(a)") "# END. "

end

!
! Subroutine func: Computes the error of the fitting  
!
! On input: ndata: The number of experimental data points
!           x: the x value of the experimental data
!           y: the y value of the experimental data for the
!              corresponding x
!           nterms: Number of terms of the exponential to
!                   to be fitted.
!           x: vector containing the parameters of the
!              exponential
!                  y =   x(1)*exp(-xdata*x(4))
!                      + x(2)*exp(-xdata*x(5))
!                      + x(3)*exp(-xdata*x(6))
!                          |               |
!           positions:     j           nterms + j
!
! On return: f: The current quadratic error of the least-square 
!               fitting
!

subroutine func(x,f) 

use problem_data
implicit none
integer i, j
double precision :: x(*), f, f_val

f = 0.d0
do i = 1, ndata
  f_val = 0.d0
  do j = 1, nterms
    f_val = f_val + x(j) * exp( -xdata(i) * x(nterms+j) )   
  end do
  f = f + ( f_val - ydata(i) )**2
end do 

return
end

!
! Subroutine grad: Computes the gradient of the error
!
! On input: ndata: number of experimental data points
!           x: the x value of the data points
!           y: the y value of the data points
!           nterms: Number of terms of the exponential fit
!           x: vector containing the parameters of the 
!                  exponential fit (variables)
!
! On return: g: the gradient
! 
!

subroutine grad(x,g)

use problem_data
implicit none
integer i, j
double precision :: x(*), g(*), f_val, dpf, glin

! Reseting the gradient vector

do i = 1, n
  g(i) = 0.d0
end do

! Computing the gradient

do i = 1, ndata

  f_val = 0.d0
  do j = 1, nterms
    f_val = f_val + x(j) * exp( -xdata(i) * x(nterms+j) )   
  end do

! The linear term
   
  glin = 2.d0 * ( f_val - ydata(i) )

! Exponential terms

  do j = 1, nterms
    dpf = glin * exp(-xdata(i) * x(nterms+j))
    g(j) = g(j) + dpf
    g(nterms+j) = g(nterms+j) &
                  - dpf * x(j) * xdata(i)
  end do

end do

return
end

!
! Subroutine that calls algencan:
!

subroutine call_algencan(n,x,l,u,f)

implicit none
integer :: n, m, iprint, ncomp, inform
double precision :: x(n), l(n), u(n), f, epsfeas, epsopt, lambda(n),&
                    cnorm, snorm, nlpsupn
logical :: coded(10), checkder, equatn(n), linear(n)

call param(epsfeas,epsopt,iprint,ncomp)

coded( 1) = .true.  ! evalf    
coded( 2) = .true.  ! evalg   
coded( 3) = .false.  ! evalh   
coded( 4) = .true.  ! evalc   
coded( 5) = .true.  ! evaljac
coded( 6) = .false.  ! evalhc  
coded( 7) = .false. ! evalfc
coded( 8) = .false. ! evalgjac
coded( 9) = .false. ! evalhl
coded(10) = .false. ! evalhlp  

m = 0
linear(1) = .true.
equatn(1) = .true.
lambda(1) = 0.5d0

checkder = .false.   

iprint = 1
call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,&
              linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)  

return
end

! Subroutine for algencan to compute the function value

subroutine evalf(n,x,f,flag)

implicit none
integer :: n, flag
double precision :: x(n), f

flag = 0

call func(x,f)

return
end

! Subroutine for algencan to compute the gradient

subroutine evalg(n,x,g,flag)

implicit none
integer :: n, flag
double precision :: x(n), g(n) 

flag = 0
call grad(x,g)

return
end

! Subroutine for algencan that computes the constraint

subroutine evalc(n,x,ind,c,flag)

implicit none
integer :: i, ind, flag, n
double precision :: c
double precision :: x(*)

flag = 0

c = -1.d0
do i = 1, n / 2
  c = c + x(i)
end do

return
end

! Subroutine for algencan that computes the gradient of the constraint

subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

implicit none
integer :: i, flag, ind, jcnnz, n
integer :: jcvar(*)
double precision :: x(*), jcval(*)

flag = 0
jcnnz = n / 2

do i = 1, n / 2
  jcvar(i) = i
  jcval(i) = 1.d0
end do

return
end

