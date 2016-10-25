!
! program jarzynski: Computes the free energy from the work performed
!                    by the spring in a Steered Molecular Dynamics
!                    simulation, performed with NAMD (reads the NAMD
!                    output). Actually, reads several log files of
!                    different simulations, computes the work in each
!                    simulation, and outputs: 
!
!                    1. The free energy as a function of distance (PMF),
!                       computed from all simulations, according to
!                       to the exponential average and also using 
!                    2. the second order cumulant expansion.
!                    3. An estimative of the error of the free energy 
!                       at each point.
!                    4. The work performed by the force in each simulation
!
! Input information: A file containing a list of namd log files
!
! Output: A table containing in the first column the streching of
!         of spring, in the second column the free energy, in the
!         third column the estimated error of the free energy and
!         in subsequent columns the work as a function of spring
!         streching for each of the namd log files provided.
!
! L. Martinez, Institute of Chemistry, State University of Campinas
!              Mar 5, 2009. leandromartinez98@gmail.com
!
! Version 16.146
!

program jarzynski

implicit none
integer :: i, ilog, status, nlogs, maxdata, idata, stride, print,&
           smdfreq
integer, allocatable :: ndata(:)
real :: cumenergy, beta, w2av, waverage, &
        conv, velocity, v(3), vnorm, wexpav, expenergy
real, allocatable :: work(:,:), x(:), y(:), z(:),&
                     fx(:), fy(:), fz(:), ts(:)
character :: dummyc
character(len=200) :: inputfile, line
character(len=200), allocatable :: namdlog(:)

call version

! Getting the name of the input file from the command line

call getarg(1,inputfile)

! Trying to open the input file

open(10,file=inputfile,action='read',iostat=status)
if( status /= 0 ) then
  write(*,*) ' ERROR: Failed to open input file. '
  stop
end if

! Counting the number of log files provided

nlogs = 0
do 
  read(10,"( a1 )",iostat=status) line
  if ( status /= 0 ) exit
  if ( line(1:1) /= "#" ) then
    nlogs = nlogs + 1
  end if
end do
rewind(10)
write(*,"( '# Number of namd log files: ', i5 )") nlogs
allocate ( namdlog(nlogs), ndata(nlogs) )

! Reading the name of the log files, and checking the number
! of force values in each of them

nlogs = 0
maxdata = 0
do 
  read(10,"( a200 )",iostat=status) line
  if ( status /= 0 ) exit
  if ( line(1:1) /= "#" .and. len_trim(line) > 0 ) then
    nlogs = nlogs + 1
    namdlog(nlogs) = line


    open(20,file=namdlog(nlogs),action="read",iostat=status)
    if(status /= 0) then
      line = namdlog(nlogs)
      write(*,*) ' ERROR: Failed to open log file: ',&
                 line(1:len_trim(line))
      stop
    end if
    ndata(nlogs) = 0
    do
      read(20,"( a4 )",iostat=status) line
      if ( status /= 0 ) exit
      if(line(1:4) == "SMD ") then
        ndata(nlogs) = ndata(nlogs) + 1
      end if
    end do     
    close(20)
    maxdata = max(ndata(nlogs),maxdata)

    line = namdlog(nlogs)
    write(*,"( '# log file ',i5,': ',tr1,a,' ndata = ',i10 )") &
               nlogs, line(1:len_trim(line)), ndata(nlogs)
  
  end if
end do
close(10)
write(*,"( '# Maximum number of force values in one log: ', i10)") maxdata

! Allocating the arrays required for computing the work in each frame

allocate ( work(nlogs,maxdata),&
           ts(maxdata), x(maxdata), y(maxdata), z(maxdata),&
           fx(maxdata), fy(maxdata), fz(maxdata) )

! Conversion factor from pN x A to kcal/mol 

conv = 1.d0 / 69.479d0

! Use SMD output jumping some of them or not (1 uses all)

stride = 1

! Looping over the log files for computing the work in each of them

do ilog = 1, nlogs

  write(*,"( '# Reading log file: ', i6 )") ilog  
  open(10,file=namdlog(ilog),action="read")
  idata = 0
  do
    read(10,"( a200 )",iostat=status) line
    if( status /= 0 ) exit
    if(line(1:4) == "SMD ") then
      idata = idata + 1
      read(line,*) dummyc, ts(idata), x(idata), y(idata), z(idata),&
                           fx(idata), fy(idata), fz(idata)
    end if
    if(line(1:18) == "Info: SMD VELOCITY") then
      read(line(19:200),*) velocity
      write(*,"( '# SMD velocity found: ', f17.10 )") velocity
    end if
    if(line(1:19) == "Info: SMD DIRECTION") then
      read(line(20:200),*) v(1), v(2), v(3)
      write(*,"( '# SMD direction found: ', 3(tr2,f10.4) )") &
            v(1), v(2), v(3)
    end if 
    if(line(1:26) == "Info: SMD OUTPUT FREQUENCY") then
      read(line(27:200),*) smdfreq
      write(*,"( '# SMD Output frequency: ', i6 )") smdfreq
    end if
  end do
  close(10)

! Normalizing the velocity vector found for this log

  vnorm = sqrt( v(1)**2 + v(2)**2 + v(3)**2 )
  v(1) = v(1) / vnorm
  v(2) = v(2) / vnorm
  v(3) = v(3) / vnorm

! Computing the cummulative work

  do i = 1, ndata(ilog), stride
    work(ilog,i) = conv * &
                   ( fx(i)*v(1) + &
                     fy(i)*v(2) + &
                     fz(i)*v(3) ) * &
                   velocity*stride*smdfreq
    if( i > 1 ) work(ilog,i) = work(ilog,i) + work(ilog,i-stride)
  end do  
  do i = ndata(ilog) + 1, maxdata
    work(ilog,i) = work(ilog,i-1)
  end do

end do

! Computing the free energy as a function of the distance

beta = 1.688934106d0

! Output table title

write(*,"( '# Column 1: Displacement of the string ',/,&
           '# Column 2: Free energy computed using exponential average. ',/,&
           '# Column 3: Variance of works. ',/,&
           '# Column 4: Free energy computed with second order cumulant ',/,&
           '#           approximation of exponentials. ',/,&
           '# Column 5 to N: Work for each individual simulation. ')")


print = 100
do i = 1, maxdata, max( maxdata / print , 1 )

! Computing work averages for this distance

  waverage = 0.d0
  wexpav = 0.d0
  w2av = 0.d0
  do ilog = 1, nlogs
    waverage = waverage + work(ilog,i)
    wexpav = wexpav + exp( - beta * work(ilog,i) ) 
    w2av = w2av + work(ilog,i)**2
  end do
  waverage = waverage / dfloat(nlogs)
  w2av = w2av / dfloat(nlogs)
  wexpav = wexpav / dfloat(nlogs)

! Computing the exponential average

  if( wexpav > 1.d-10 ) then
    expenergy = -log( wexpav ) / beta
  else
    expenergy = 0.d0
  end if

! Computing the cumulant expansion - second order

  cumenergy = waverage - 0.5d0 * beta * ( w2av - waverage**2 )

! Writting data for this line

  write(*,"( 1000(tr2,e12.5) )") &
        velocity*ts(i), expenergy, sqrt(w2av - waverage**2), &
        cumenergy, (work(ilog,i),ilog=1,nlogs)  

end do

end
