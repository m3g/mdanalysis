C
C lm-Version 16.146
C
C     ******************************************************************
C     ******************************************************************

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalh(n,x,hlin,hcol,hval,hnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     This subroutine might compute the Hessian matrix of the objective 
C     function. For achieving this objective YOU MAY MODIFY it according 
C     to your problem. To modify this subroutine IS NOT MANDATORY. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     hnnz     integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hlin     integer hlin(hnnz),
C              see below,
C
C     hcol     integer hcol(hnnz),
C              see below,
C
C     hval     double precision hval(hnnz),
C              the non-null value of the (hlin(k),hcol(k)) position 
C              of the Hessian matrix of the objective function must 
C              be saved at hval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the objective funtion. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
C     HESSIAN MATRIX OF YOUR OBJECTIVE FUNCTION: 
C     ******************************************************************

      flag = 0

      hnnz = 0 

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALH. 
C     ******************************************************************
 
      end

c    C     ******************************************************************
c    C     ******************************************************************
c     
c          subroutine evalc(n,x,ind,c,flag)
c    
c          implicit none
c    
c    C     SCALAR ARGUMENTS
c          integer ind,flag,n
c          double precision c
c    
c    C     ARRAY ARGUMENTS
c          double precision x(n)
c    
c    C     This subroutine must compute the ind-th constraint of your 
c    C     problem. For achieving this objective YOU MUST MOFIFY it 
c    C     according to your problem. See below the places where your 
c    C     modifications must be inserted.
c    C
c    C     Parameters of the subroutine:
c    C
c    C     On Entry:
c    C
c    C     n        integer,
c    C              number of variables,
c    C
c    C     x        double precision x(n),
c    C              current point,
c    C
c    C     ind      integer,
c    C              index of the constraint to be computed,
c    C
c    C     On Return
c    C
c    C     c        double precision,
c    C              ind-th constraint evaluated at x,
c    C
c    C     flag     integer
c    C              You must set it to any number different of 0 (zero) if 
c    C              some error ocurred during the evaluation of the 
c    C              constraint. (For example, trying to compute the square 
c    C              root of a negative number, dividing by zero or a very 
c    C              small number, etc.) If everything was o.k. you must set 
c    C              it equal to zero.
c     
c    C     Constraints
c    
c    C     ******************************************************************
c    C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR CONSTRAINTS 
c    C     ******************************************************************
c    
c          flag = 0
c    
c          if ( ind .eq. 1 ) then
c              c = x(1) ** 2 + 1 - x(n)
c              return
c          end if
c    
c          if ( ind .eq. 2 ) then
c              c = 2.0d0 - x(1) - x(n)
c              return
c          end if
c    
c    C     ******************************************************************
c    C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. 
c    C     ******************************************************************
c     
c          end
c    
c    C     ******************************************************************
c    C     ******************************************************************
c     
c          subroutine evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)
c    
c          implicit none
c    
c    C     SCALAR ARGUMENTS
c          integer flag,ind,jcnnz,n
c    
c    C     ARRAY ARGUMENTS
c          integer jcvar(n)
c          double precision x(n),jcval(n)
c    
c    C     This subroutine must compute the gradient of the ind-th constraint.
c    C     For achieving these objective YOU MUST MODIFY it in the way 
c    C     specified below.
c    C
c    C     Parameters of the subroutine:
c    C
c    C     On Entry:
c    C
c    C     n        integer,
c    C              number of variables,
c    C
c    C     x        double precision x(n),
c    C              current point,
c    C
c    C     ind      integer,
c    C              index of the constraint whose gradient will be computed,
c    C
c    C     On Return
c    C
c    C     jcnnz    integer,
c    C              number of perhaps-non-null elements of the computed 
c    C              gradient,
c    C
c    C     jcvar    integer jcvar(jcnnz),
c    C              see below,
c    C
c    C     jcval    double precision jcval(jcnnz),
c    C              the non-null value of the partial derivative of the 
c    C              ind-th constraint with respect to the jcvar(k)-th 
c    C              variable must be saved at jcval(k).
c    C
c    C     flag     integer
c    C              You must set it to any number different of 0 (zero) if 
c    C              some error ocurred during the evaluation of the 
c    C              constraint. (For example, trying to compute the square 
c    C              root of a negative number, dividing by zero or a very 
c    C              small number, etc.) If everything was o.k. you must set 
c    C              it equal to zero.
c    
c    C     Sparse gradient vector of the ind-th constraint
c    
c    C     ******************************************************************
c    C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS  
c    C     OF YOUR CONSTRAINTS: 
c    C     ******************************************************************
c    
c          flag = 0
c    
c          if ( ind .eq. 1 ) then
c    
c              jcnnz = 2
c    
c              jcvar(1) = 1
c              jcval(1) = 2.0d0 * x(1)
c              jcvar(2) = 2
c              jcval(2) = - 1.0d0
c    
c              return
c    
c          end if
c    
c          if ( ind .eq. 2 ) then
c    
c              jcnnz = 2
c    
c              jcvar(1) = 1
c              jcval(1) = - 1.0d0
c              jcvar(2) = 2
c              jcval(2) = - 1.0d0
c    
c              return
c    
c          end if
c    
c    C     ******************************************************************
c    C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. 
c    C     ******************************************************************
c     
c          end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hcnnz,ind,n

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     This subroutine might compute the Hessian matrix of the ind-th
C     constraint. For achieving this objective YOU MAY MODIFY it 
C     according to your problem. To modify this subroutine IS NOT 
C     MANDATORY. See below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose Hessian will be computed,
C
C     On Return
C
C     hcnnz    integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hclin    integer hclin(hcnnz),
C              see below,
C
C     hccol    integer hccol(hcnnz),
C              see below,
C
C     hcval    double precision hcval(hcnnz),
C              the non-null value of the (hclin(k),hccol(k)) position 
C              of the Hessian matrix of the ind-th constraint must 
C              be saved at hcval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the ind-th constraint. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
C     HESSIANS OF YOUR CONSTRAINTS: 
C     ******************************************************************

      flag = 0

      if ( ind .eq. 1 ) then

          hcnnz = 1

          hclin(1) = 1
          hccol(1) = 1
          hcval(1) = 2.0d0

          return

      end if

      if ( ind .eq. 2 ) then

          hcnnz = 0

          return

      end if

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalfc(n,x,f,m,c,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhl(n,x,m,lambda,scalef,scalec,hllin,hlcol,hlval,
     +hlnnz,flag)

      implicit none

C     SCALAR ARGUMENTS
      integer flag,hlnnz,m,n
      double precision scalef

C     ARRAY ARGUMENTS
      integer hlcol(*),hllin(*)
      double precision hlval(*),lambda(m),scalec(m),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

C     This subroutine might compute the product of the Hessian of the
C     Lagrangian times vector p (just the Hessian of the objective 
C     function in the unconstrained or bound-constrained case). 
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     m        integer,
C              number of constraints,
C
C     lambda   double precision lambda(m),
C              vector of Lagrange multipliers,
C
C     p        double precision p(n),
C              vector of the matrix-vector product,
C
C     goth     logical,
C              can be used to indicate if the Hessian matrices were
C              computed at the current point. It is set to .false.
C              by the optimization method every time the current
C              point is modified. Sugestion: if its value is .false. 
C              then compute the Hessians, save them in a common 
C              structure and set goth to .true.. Otherwise, just use 
C              the Hessians saved in the common block structure,
C
C     On Return
C
C     hp       double precision hp(n),
C              Hessian-vector product,
C
C     goth     logical,
C              see above,
C              
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              Hessian-vector product. (For example, trying to compute 
C              the square root of a negative number, dividing by zero 
C              or a very small number, etc.) If everything was o.k. you 
C              must set it equal to zero.

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine endp(n,x,l,u,m,lambda,equatn,linear)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     This subroutine can be used to do some extra job after the solver
C     has found the solution,like some extra statistics, or to save the
C     solution in some special format or to draw some graphical
C     representation of the solution. If the information given by the
C     solver is enough for you then leave the body of this subroutine
C     empty.
C     
C     Parameters of the subroutine:
C
C     The paraemters of this subroutine are the same parameters of
C     subroutine inip. But in this subroutine there are not output
C     parameter. All the parameters are input parameters.

      end

cC     ******************************************************************
cC     ******************************************************************
c
c      program algencanma
c
c      implicit none
c
cC     PARAMETERS
c
c      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax
c
c      parameter ( mmax      =   500000 )
c      parameter ( nmax      =   500000 )
c      parameter ( nsmax     =     1000 )
c      parameter ( jcnnzmax  = 10000000 )
c      parameter ( hnnzmax   = 10000000 )
c      parameter ( fnnzmax   = 10000000 )
c      parameter ( wintmax   = 10000000 )
c      parameter ( nsysmax   =   100000 )
c
cC     LOCAL SCALARS
c      logical checkder
c      integer inform,iprint,m,n,ncomp
c      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm
c
cC     LOCAL ARRAYS
c      logical coded(10),equatn(mmax),linear(mmax)
c      double precision l(nmax),lambda(mmax),u(nmax),x(nmax)
c
cC     EXTERNAL SUBROUTINES
c      external algencan,param
c
cC     SET UP PROBLEM DATA
c
c      call param(epsfeas,epsopt,iprint,ncomp)
c
c      call inip(n,x,l,u,m,lambda,equatn,linear,coded,checkder)
c
c      call algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,equatn,
c     +linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)
c
c      call endp(n,x,l,u,m,lambda,equatn,linear)
c
c      stop
c
c      end

C     ******************************************************************
C     ******************************************************************

      subroutine param(epsfeas,epsopt,iprint,ncomp)

C     SCALAR ARGUMENTS
      integer iprint,ncomp
      double precision epsfeas,epsopt

      epsfeas  = 1.0d-08
      epsopt   = 1.0d-08

      iprint   = 10
      ncomp    = 6

      end

C     ******************************************************************
C     ******************************************************************

      subroutine algencan(epsfeas,epsopt,iprint,ncomp,n,x,l,u,m,lambda,
     +equatn,linear,coded,checkder,f,cnorm,snorm,nlpsupn,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,iprint,m,n,ncomp
      double precision cnorm,epsfeas,epsopt,f,nlpsupn,snorm

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer alinfo,geninfo,i,iter,j,maxit,nwcalls,nwtotit,outiter,
     +        solinfo,totiter,msqcalls,msqtotit
      double precision cnormb,cnormu,cnormub,fb,fu,fub,nlpsupnb,seed,
     +        drand
      real time
 
C     LOCAL ARRAYS
      double precision nl(nmax),c(mmax),rho(mmax)
      real dum(2)

C     DATA STATEMENTS
      data dum/0.0,0.0/

C     EXTERNAL FUNCTIONS AND SUBROUTINES
      external fparam,checkd,auglag,gencan,drand

C     ==================================================================
C     Start timing
C     ==================================================================

      time = dtime(dum)

C     ==================================================================
C     Initialization
C     ==================================================================

C     Set machine-dependent constants

      bignum    = 1.0d+99
      macheps   = 1.0d-16
      macheps12 = sqrt( macheps )
      macheps13 = macheps ** ( 1.0d0 / 3.0d0 )
      macheps23 = macheps ** ( 2.0d0 / 3.0d0 )

C     Set global counters

      fcnt    = 0
      efcnt   = 0
      efccnt  = 0
      egcnt   = 0
      egjccnt = 0
      ehcnt   = 0
      ehlcnt  = 0
      ehlpcnt = 0

      do j = 1,m
          eccnt(j)  = 0
          ejccnt(j) = 0
          ehccnt(j) = 0
      end do

C     ==================================================================
C     Set default values for algoritmic parameters
C     ==================================================================

C     Set user-provided subroutines indicators

      fcoded    = coded(1)
      gcoded    = coded(2)
      hcoded    = coded(3)
      ccoded    = coded(4)
      jaccoded  = coded(5)
      hccoded   = coded(6)
      fccoded   = coded(7)
      gjaccoded = coded(8)
      hlcoded   = coded(9)
      hlpcoded  = coded(10)

      innercall = .false.
      useustp   = .false.

C     Set indicator of whether the true Hessian of the Lagrangian can be
C     computed or not

      truehl = .false.
      if ( hlcoded .or. ( hcoded .and. ( hccoded .or. m .eq. 0 ) )) then
          truehl = .true.
      end if

C     Hessian-vector product strategy: HAPPRO, INCQUO or TRUEHL (TRUEHL
C     is the default option. If the proper subroutines were not coded by
C     the user, then HAPPRO is used instead.)

      if ( truehl .or. hlpcoded ) then
          hptype = 'TRUEHL'
      else
          hptype = 'HAPPRO'
      end if

C     Ignore objective function (to only find a feasible point by 
C     minimizing 1/2 of the squared infeasibility)
      ignoref = .false.

C     Inner-to-the-face minimization algorithm (CG is the default option)
      avoidds = .true.

C     Skip acceleration step
      skipacc = .true.

C     Scaling of linear systems
      sclsys = .false.

C     Slacks for inequality constraints
      slacks = .false.

C     Remove fixed variables (with identical lower and upper bounds)
      rmfixv = .true.

C     Scale objective function and constraints
      if ( m .gt. 0 ) then
          scale = .true.
      else
          scale = .false.
      end if

C     Main output control (silent-mode?)

      iprintctl(1) = .true. ! Banner
      iprintctl(2) = .true. ! Parameters and problem processing
      iprintctl(3) = .true. ! Warnings and errors messages
      iprintctl(4) = .true. ! Screen-mirror file algencan.out
      iprintctl(5) = .true. ! Solution file solution.txt
      iprintctl(6) = .true. ! Statistic files with table lines

      open(10,err=100,file='/dev/null')
      close(10)

      do i = 1,6
          iprintctl(i) = .false.
      end do

      iprint = 0

 100  continue

      if ( iprintctl(4) ) then
          open(unit=10,file='algencan.out',status='replace')
      else
          open(unit=10,                    status='scratch')
      end if

C     ==================================================================
C     Set solver arguments using the specification file
C     ==================================================================

      call fparam(epsfeas,epsopt,iprint,ncomp)

C     Outer and inner iterations output detail

      iprintout = iprint / 10
      iprintinn = mod( iprint, 10 )

C     Error tracker

      inform = 0

C     ==================================================================
C     Initialize problem data structures
C     ==================================================================

      call sinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

      nprint = min( n, ncomp )
      mprint = min( m, ncomp )

c      seed = 123456.0d0
c      do i = 1,n
c          x(i) = x(i) + max( abs( x(i) ), 1.0d-08 ) * 
c     +                  ( 2.0d-8 * drand(seed) - 1.0d-08 )
c          x(i) = x(i) + 1.0d-08 * 
c     +                  ( 2.0d-8 * drand(seed) - 1.0d-08 )
c          x(i) = max( l(i), min( x(i), u(i) ) )
c      end do

C     ==================================================================
C     Call the solver
C     ==================================================================

C     ALGENCAN for PNL problems

      if ( .not. ignoref .and. m .gt. 0 ) then
          call auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,f,c,
     +    cnorm,snorm,nl,nlpsupn,fu,cnormu,fub,cnormub,fb,cnormb,
     +    nlpsupnb,outiter,totiter,nwcalls,nwtotit,msqcalls,msqtotit,
     +    alinfo,inform)

          solinfo = alinfo

C     GENCAN for box-constrained problems

      else
          maxit = 1000

C         Used in feasibility problems (ignoref=true). With lambda=0 and
C         rho=1, to minimize 1/2 of the squared infeasibility coincides
C         with minimizing the augmented Lagrangian.
          do j = 1,m
              lambda(j) = 0.0d0
              rho(j)    = 1.0d0
          end do

          call gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas,epsopt,
     +    maxit,iter,f,nl,nlpsupn,cnorm,cnormu,geninfo,inform)

          solinfo  = geninfo

          outiter  = 0
          totiter  = iter
          nwcalls  = 0
          nwtotit  = 0
          msqcalls = 0
          msqtotit = 0

          if ( scale ) then
              fu = f * sf
          else
              fu = f
          end if

          fb       = f
          fub      = fu
          cnormb   = cnorm
          cnormub  = cnormu
          nlpsupnb = nlpsupn
      end if

      if ( inform .lt. 0 ) return

C     Close output file

      close(10)

C     ==================================================================
C     End problem data structures
C     ==================================================================

      call sendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     Stop timing
C     ==================================================================

      time = dtime(dum)
      time = dum(1)

      if ( iprintout .eq. 1 ) then
          write(*,9000) time
      end if

C     ==================================================================
C     Write statistics
C     ==================================================================

      if ( iprintctl(6) ) then
          open(20,file='algencan-tabline.out')
          write(20,9040) fu,cnormu,f,cnorm,nlpsupn,fub,cnormub,fb,
     +                   cnormb,nlpsupnb,inform,solinfo,n,m,outiter,
     +                   totiter,fcnt,nwcalls,nwtotit,msqcalls,msqtotit,
     +                   time
          close(20)
      end if

C     ==================================================================
C     NON-EXECUTABLE STATEMENTS
C     ==================================================================

 9000 format(/,1X,'Total CPU time in seconds: ',F8.2)
 9040 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,I3,1X,I1,1X,I6,1X,I6,1X,I2,1X,I7,1X,I7,1X,I2,1X,I7,
     +       1X,I7,1X,I7,0P,F8.2)


      end
C     ******************************************************************
C     ******************************************************************

      subroutine auglag(n,x,l,u,m,lambda,equatn,linear,epsfeas,epsopt,
     +f,c,cnorm,snorm,nl,nlpsupn,fu,cnormu,fub,cnormub,fb,cnormb,
     +nlpsupnb,outiter,totiter,nwcalls,nwtotit,msqcalls,msqtotit,alinfo,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer alinfo,inform,m,msqcalls,msqtotit,n,nwcalls,nwtotit,
     +        outiter,totiter
      double precision cnorm,cnormb,cnormu,cnormub,epsfeas,epsopt,f,fb,
     +        fu,fub,nlpsupn,nlpsupnb,snorm

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision c(m),l(n),lambda(m),nl(n),rho(m),u(n),x(n)

C     Solves the nonlinear programming problem
C
C     min f(x)
C
C     subject to
C
C             c_j(x)  = 0, j in E,
C             c_j(x) <= 0, j in I,
C             l <= x <= u,
C
C     where E is the set of indices of the equality constraints, I is
C     the set of indices of the inequality constraints, and there are
C     n variables and m constraints, using the method of multipliers
C     described in
C
C     R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
C     "On Augmented Lagrangian methods with general lower-level
C     constraints", SIAM Journal on Optimization 18, pp. 1286-1309, 2007.

C     alinfo:
C
C     0: Feasibility, optimality and complementarity satisfied
C     1: Maximum number of algencan iterations reached
C     2: Infeasible problem?

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      character innstp
      logical scaletmp
      integer geninfo,i,iter,maxit,msqiter,nwinfo,nwiter,outitnp,
     +        rhoincr,xind
      double precision al,cnorma,epsfeas12,epsfeas14,epsopt12,epsopt14,
     +        epsopk,fa,nlpi,nlpsupna,rsupn,snormb,snormprev,xsupn

C     LOCAL ARRAYS
      double precision lambar(mmax),lambdaa(mmax),lambdab(mmax),
     +        sigma(mmax),sigpre(mmax),xa(nmax),xb(nmax)

C     DATA BLOCKS
      character genstp(0:9)
      data genstp(0) /'C'/
      data genstp(1) /'M'/
      data genstp(2) /' '/
      data genstp(3) /'F'/
      data genstp(4) /'G'/
      data genstp(5) /'P'/
      data genstp(6) /'U'/
      data genstp(7) /'S'/
      data genstp(8) /' '/
      data genstp(9) /' '/

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprintout .ge. 1 ) then
          write(* ,1000) n,m
          write(10,1000) n,m

          if ( iprintout .ge. 4 .and. nprint .ne. 0 ) then
              write(* ,1010) nprint,(l(i),i=1,nprint)
              write(* ,1020) nprint,(u(i),i=1,nprint)

              write(10,1010) nprint,(l(i),i=1,nprint)
              write(10,1020) nprint,(u(i),i=1,nprint)
          end if

      end if

C     ==================================================================
C     Initialization
C     ==================================================================

C     Counters

      outiter = 0
      outitnp = 0

      totiter = 0

      nwcalls = 0
      nwtotit = 0

      msqcalls = 0
      msqtotit = 0

C     Constants

      epsopt12  = sqrt( epsopt    )
      epsopt14  = sqrt( epsopt12  )
      epsfeas12 = sqrt( epsfeas   )
      epsfeas14 = sqrt( epsfeas12 )

C     Project initial point

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute current point norm

      xsupn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind  = i
          end if
      end do

C     Compute objective function and constraints

      call sevalobjc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

C     Set penalty parameters and compute maximum rho

      if ( rhoauto ) then
          call comprhoini(f,m,c,equatn,rho)
      end if

      rsupn = 0.0d0
      do i = 1,m
          rsupn = max( rsupn, rho(i) )
      end do

C     Compute safeguarded Lagrange multipliers

      do i = 1,m
          lambar(i) = max( lammin, min( lambda(i), lammax ) )
      end do

C     Compute complementarity and feasibility violations

      snormprev = bignum

      cnorm  = 0.0d0
      snorm  = 0.0d0
      do i = 1,m
          sigpre(i) = bignum

          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( c(i) ) )
              sigma(i) = c(i)

          else
              cnorm = max( cnorm, c(i) )
              sigma(i) = max( c(i), - lambar(i) / rho(i) )
          end if

          snorm = max( snorm, abs( sigma(i) ) )
      end do

C     Compute unscaled objective function and constraints norm

      if ( scale ) then
          fu = f * sf

          cnormu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  cnormu = max( cnormu, abs( sc(i) * c(i) ) )
              else
                  cnormu = max( cnormu, sc(i) * c(i) )
              end if
          end do

      else
          fu     = f
          cnormu = cnorm
      end if

C     Compute continuous projected Lagrangian gradient norm

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

      nlpsupn = 0.0d0
      do i = 1,n
          nlpi = x(i) - nl(i)
          if ( l(i) .le. nlpi .and. nlpi .le. u(i) ) then
              nlpi = - nl(i)
          else
              nlpi = max( l(i), min( nlpi, u(i) ) ) - x(i)
          end if
          nlpsupn = max( nlpsupn, abs( nlpi ) )
      end do

C     Save best solution

      fub      = fu
      cnormub  = cnormu

      fb       = f
      cnormb   = cnorm
      snormb   = snorm
      nlpsupnb = nlpsupn

      do i = 1,n
          xb(i) = x(i)
      end do

      do i = 1,m
          lambdab(i) = lambda(i)
      end do

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Print information of this iteration
C     ==================================================================


      if ( iprintout .eq. 1 ) then

          if ( outiter .gt. 0 ) then
              innstp = genstp(geninfo)
          else
              innstp = ' '
          end if

          if ( .not. scale ) then
              if ( mod(outiter,10) .eq. 0 ) then
                  write(* ,1030)
                  write(10,1030)
              end if
              write(* ,1040) outiter,f,cnorm,snorm,nlpsupn,xsupn,rsupn,
     +                       totiter,innstp,nwcalls,nwtotit
              write(10,1040) outiter,f,cnorm,snorm,nlpsupn,xsupn,rsupn,
     +                       totiter,innstp,nwcalls,nwtotit
          else
              if ( mod(outiter,10) .eq. 0 ) then
                  write(* ,1050)
                  write(10,1050)
              end if
              write(* ,1060) outiter,fu,cnormu,f,cnorm,snorm,nlpsupn,
     +                       rsupn,totiter,innstp,nwcalls,nwtotit
              write(10,1060) outiter,fu,cnormu,f,cnorm,snorm,nlpsupn,
     +                       rsupn,totiter,innstp,nwcalls,nwtotit
          end if

      else if ( iprintout .ge. 2 ) then
          write(* ,1070) outiter
          write(10,1070) outiter

          if ( iprintout .ge. 3 ) then
              write(* ,1080) totiter,nwcalls,nwtotit
              write(10,1080) totiter,nwcalls,nwtotit
          end if

          if ( .not. scale ) then
              write(* ,1090) f,cnorm
              write(10,1090) f,cnorm
          else
              write(* ,1100) f,fu,cnorm,cnormu
              write(10,1100) f,fu,cnorm,cnormu
          end if

          write(* ,1110) snorm,nlpsupn,xind,xsupn,rsupn
          write(10,1110) snorm,nlpsupn,xind,xsupn,rsupn

          if ( iprintout .ge. 4 ) then
              if ( nprint .ne. 0 ) then
                  write(* ,1120) nprint,(x(i),i=1,nprint)
                  write(10,1120) nprint,(x(i),i=1,nprint)
              end if

              if ( mprint .ne. 0 ) then
                  write(* ,1130) mprint,(lambda(i),i=1,mprint)
                  write(* ,1140) mprint,(rho(i),i=1,mprint)

                  write(10,1130) mprint,(lambda(i),i=1,mprint)
                  write(10,1140) mprint,(rho(i),i=1,mprint)
              end if
          end if

      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          open(20,file='algencan-tabline.out')
          write(20,1400) fu,cnormu,f,cnorm,nlpsupn,fub,cnormub,fb,
     +                   cnormb,nlpsupnb,inform,9,n,m,outiter,totiter,0,
     +                   nwcalls,nwtotit,msqcalls,msqtotit,999.9d0
          close(20)
      end if

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     Test feasibility, optimality and complementarity

      if ( max( snorm, cnormu ) .le. epsfeas .and.
     +     nlpsupn .le. epsopt ) then
          alinfo = 0

          if ( iprintout .ge. 1 ) then
              write(*, 1300)
              write(10,1300)
          end if

          return
      end if

C     Test whether the number of iterations is exhausted

      if ( outiter .ge. maxoutit ) then
          alinfo = 1

          if ( iprintout .ge. 1 ) then
              write(*, 1310)
              write(10,1310)
          end if

          return
      end if

C     Test whether the problem seems to be infeasible

      if ( snorm .gt. epsfeas .and. snorm .ge. snormprev ) then
          outitnp = outitnp + 1
      else
          outitnp = 0
      end if

      if ( outitnp .ge. maxoutitnp .or. rsupn .gt. rhomax ) then

          alinfo = 2

          if ( iprintout .ge. 1 ) then
              write(*, 1320)
              write(10,1320)
          end if

          return
      end if

C     ==================================================================
C     Near the solution, try to solve the KKT system by Newton's method
C     ==================================================================

      if ( .not. skipacc .and. truehl .and. ( nwcalls .gt. 0 .or.
     +     ( snorm .le. epsfeas12 .and. nlpsupn .le. epsopt12 ) .or.
     +     ( snorm .le. epsfeas14 .and. nlpsupn .le. epsopt14 .and.
     +       outiter .gt. 0 .and. geninfo .ne. 0 ) ) ) then

          nwcalls = nwcalls + 1

          do i = 1,n
              xa(i) = x(i)
          end do

          do i = 1,m
              lambdaa(i) = lambar(i)
          end do

          scaletmp = scale

          if ( scale ) then
              do i = 1,m
                  lambdaa(i) = lambdaa(i) * sf / sc(i)
              end do

              scale = .false.
          end if

          call newtonkkt(n,xa,l,u,m,lambdaa,equatn,linear,epsfeas,
     +    epsopt,fa,cnorma,nlpsupna,nwiter,msqiter,nwinfo,inform)

          nwtotit = nwtotit  + nwiter

          if ( nwinfo .ge. 6 ) then
              skipacc = .true.
          end if

          if ( msqiter .gt. 0 ) then
              msqcalls = msqcalls + 1
              msqtotit = msqtotit + msqiter
          end if

          scale = scaletmp

          if ( inform .lt. 0 ) return

C         Save best solution
 
          if ( ( cnormub .gt. epsfeas .and. cnorma .le. cnormub ) .or. 
     +         ( cnormub .le. epsfeas .and. cnorma .le. epsfeas .and. 
     +           fa .le. fub ) ) then

              fub      = fa
              cnormub  = cnorma

              fb       = fa
              cnormb   = cnorma
              snormb   = 0.0d0
              nlpsupnb = nlpsupna

              do i = 1,n
                  xb(i) = xa(i)
              end do

              do i = 1,m
                  lambdab(i) = lambdaa(i)
              end do

          end if

C         Test feasibility, optimality and complementarity

          if ( cnorma .le. epsfeas .and. nlpsupna .le. epsopt ) then

              fu      = fa
              cnormu  = cnorma

              f       = fa
              cnorm   = cnorma
              snorm   = 0.0d0
              nlpsupn = nlpsupna

              do i = 1,n
                  x(i) = xa(i)
              end do

              do i = 1,m
                  lambda(i) = lambdaa(i)
              end do

              alinfo = 0

              if ( iprintout .ge. 1 ) then
                  write(*, 1300)
                  write(10,1300)
              end if

              return

          end if

      end if

C     ==================================================================
C     Iteration
C     ==================================================================

      outiter = outiter + 1

C     ==================================================================
C     Solve the augmented Lagrangian subproblem
C     ==================================================================

C     Set optimality requeriment for the subproblem

      if ( outiter .gt. 1 .and. 
     +     snorm .le. epsfeas12 .and. nlpsupn .le. epsopt12 ) then
          epsopk = max( epsopt, 0.1d0 * epsopk )
      else
          epsopk = sqrt( epsopt )
      end if

      if ( outiter .eq. 1 ) then
          maxit =   10
      else
          maxit = 1000
      end if

C     Call the inner-solver

      call gencan(n,x,l,u,m,lambar,equatn,linear,rho,epsfeas,epsopk,
     +maxit,iter,al,nl,nlpsupn,cnorm,cnormu,geninfo,inform)

      totiter = totiter + iter

      if ( inform .lt. 0 ) return

C     ==================================================================
C     Prepare for the next iteration
C     ==================================================================

C     Compute current point norm

      xsupn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind  = i
          end if
      end do

C     Compute objective function and constraints

      call sevalobjc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

C     Compute complementarity and feasibility violations

      snormprev = snorm

      cnorm  = 0.0d0
      snorm  = 0.0d0
      do i = 1,m
          sigpre(i) = sigma(i)

          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( c(i) ) )
              sigma(i) = c(i)

          else
              cnorm = max( cnorm, c(i) )
              sigma(i) = max( c(i), - lambar(i) / rho(i) )
          end if

          snorm = max( snorm, abs( sigma(i) ) )
      end do

C     Compute unscaled objective function and constraints norm

      if ( scale ) then
          fu = f * sf

          cnormu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  cnormu = max( cnormu, abs( sc(i) * c(i) ) )
              else
                  cnormu = max( cnormu, sc(i) * c(i) )
              end if
          end do

      else
          fu     = f
          cnormu = cnorm
      end if

C     Update Lagrange multipliers approximation

      do i = 1,m
          call evaldpdy(c(i),rho(i),lambar(i),equatn(i),lambda(i))
          lambar(i) = max( lammin, min( lambda(i), lammax ) )
      end do

C     Compute continuous projected Lagrangian gradient norm

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

      nlpsupn = 0.0d0
      do i = 1,n
          nlpi = x(i) - nl(i)
          if ( l(i) .le. nlpi .and. nlpi .le. u(i) ) then
              nlpi = - nl(i)
          else
              nlpi = max( l(i), min( nlpi, u(i) ) ) - x(i)
          end if
          nlpsupn = max( nlpsupn, abs( nlpi ) )
      end do

C     Save best solution

      if ( ( cnormub .gt. epsfeas .and. cnormu .le. cnormub ) .or. 
     +     ( cnormub .le. epsfeas .and. cnormu .le. epsfeas .and. 
     +       f .lt. fb ) ) then

          fub      = fu
          cnormub  = cnormu

          fb       = f
          cnormb   = cnorm
          snormb   = snorm
          nlpsupnb = nlpsupn

          do i = 1,n
              xb(i) = x(i)
          end do

          do i = 1,m
              lambdab(i) = lambda(i)
          end do

      end if

C     Update penalty parameters

      if ( outiter .eq. 1 ) then

          call comprhoini(f,m,c,equatn,rho)

          if ( iprintout .ge. 5 ) then
              write(*, 1200)
              write(10,1200)
          end if

      else

C         Note that |sigma(i)| <= eps implies c(i) <= eps.

          if ( rhoiden ) then

              if ( max( snorm, cnormu ) .gt. epsfeas .and.
     +             snorm .gt. rhofrac * snormprev ) then

                  do i = 1,m
                      rho(i) = rhomult * rho(i)
                  end do

                  if ( iprintout .ge. 5 ) then
                      write(*, 1210) rhomult
                      write(10,1210) rhomult
                  end if

              else

                  if ( iprintout .ge. 5 ) then
                      write(*, 1230)
                      write(10,1230)
                  end if

              end if

          else

              rhoincr = 0
              do i = 1,m
                  if ( max(abs(sigma(i)), sc(i)*c(i)) .gt. epsfeas .and.
     +                 abs(sigma(i)) .gt. rhofrac*abs(sigpre(i)) ) then
                      rho(i)  = rhomult * rho(i)
                      rhoincr = rhoincr + 1
                  end if
              end do

              if ( iprintout .ge. 5 ) then

                  if ( rhoincr .ne. 0 ) then
                      write(*, 1220) rhoincr,m,rhomult
                      write(10,1220) rhoincr,m,rhomult
                  else
                      write(*, 1230)
                      write(10,1230)
                  end if

              end if

          end if

      end if

C     Compute maximum rho

      rsupn = 0.0d0
      do i = 1,m
          rsupn = max( rsupn, rho(i) )
      end do

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,1X,'Entry to ALGENCAN.',
     +       /,1X,'Number of variables  : ',I7,
     +       /,1X,'Number of constraints: ',I7)
 1010 format(/,1X,'Lower bounds  (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1020 format(/,1X,'Upper bounds  (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))

 1030 format(/,1X,'It',15X,'f',4X,'cnrm',4X,'snrm',2X,'nlpnrm',4X,
     +         'xnrm',1X,'penalty',2X,'innit',1X,'nwcalls',2X,'nwit')
 1040 format(I3,1P,D16.8,5(1X,1P,D7.1),1X,I5,A1,2(1X,I6))

 1050 format(/,1X,'It',4X,'unscaled f & cnrm',6X,'scaled f & cnrm',
     +         3X,'snrm',1X,'nlpnrm',1X,'penalty',1X,'innit',2X,'nwkkt')
 1060 format(I3,2(1P,D14.6,1X,1P,D6.0),3(1X,1P,D6.0),1X,I5,A1,1X,I2,1X,
     +         I3)

 1070 format(/,1X,'ALGENCAN OUTER ITERATION                      = ',I7)

 1080 format(/,1X,'Up-to-now total number of iterations          = ',I7,
     +       /,1X,'Up-to-now acceleration trials                 = ',I7,
     +       /,1X,'Up-to-now total number of Newton iterations   = ',I7)

 1090 format(/,1X,'Functional value                              = ',
     +             1P,D24.16,
     +       /,1X,'Sup-norm of constraints                       = ',
     +             17X,1P,D7.1)

 1100 format(/,1X,'Functional value          (scaled = ',1P,D8.1,
     +            ') = ',1P,D24.16,
     +       /,1X,'Sup-norm of constraints   (scaled = ',1P,D8.1,
     +            ') = ',17X,1P,D7.1)

 1110 format(  1X,'Sup-norm of complementarity-feasibility       = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of the Lagrangian projected gradient = ',
     +             17X,1P,D7.1,
     +       /,1X,'Sup-norm of x (attained at x_{', I7,'})       = ',
     +             17X,1P,D7.1,
     +       /,1X,'Largest penalty parameter                     = ',
     +             17X,1P,D7.1)

 1120 format(/,1X,'Current point (first ',I7,' components):',
     +       /,6(1X,1P,D11.4))
 1130 format(/,1X,'Updated Lagrange multipliers (first ',I7,
     +         1X,'components):',
     +       /,6(1X,1P,D11.4))
 1140 format(/,1X,'Updated penalty parameters (first ',I7,
     +         1X,'components):',
     +       /,6(1X,1P,D11.4))

 1200 format(/,1X,'Penalty parameters are re-initiated after the ',
     +            'resolution of the first',/,1X,'subproblem.')
 1210 format(/,1X,'The desired infeasibility was not achieved. ',
     +            'The penalty',/,1X,'parameter will be ',
     +            'increased multiplying by rhofrac = ',1PD11.4,'.')
 1220 format(/,1X,'The desired infeasibility was not achieved in ',I7,
     +         1X,'over ',I7,/,1X,'constraints.',/,1X,' Penalty ',
     +            'parameters will be increased multiplying by ',
     +            'rhofrac = ',1PD11.4,'.')
 1230 format(/,1X,'Desired feasibility improvement was achieved, ',
     +            'so, penalty parameters will not',/,1X,'be ',
     +            'modified.')

 1300 format(/,1X,'Flag of ALGENCAN: Solution was found.')
 1310 format(/,1X,'Flag of ALGENCAN: Maximum of iterations reached.')
 1320 format(/,1X,'Flag of ALGENCAN: The problem seems to be ',
     +            'infeasible.')

 1400 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,I3,1X,I1,1X,I6,1X,I6,1X,I2,1X,I7,1X,I7,1X,I2,1X,I7,
     +       1X,I7,1X,I7,0P,F8.2)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine comprhoini(f,m,c,equatn,rho)

C     SCALAR ARGUMENTS
      integer m
      double precision f

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision c(m),rho(m)

C     Consider the Augmented Lagrangian function
C
C         al = f(x) + \sum_i P(lambda_i,c_i(x),rho),
C
C     where
C
C         P(lambda,y,rho) = y ( lambda + 0.5 rho y ),
C
C     If c_i(x) is an equality constraint or lambda + rho y > 0, and
C
C         P(lambda,y,rho) = - 0.5 lambda^2 / rho,
C
C     otherwise.
C
C     Assuming that lambda_i = 0 for all i, it is clear that
C
C         P(lambda_i,c_i(x),rho) = 0.5 rho c_i(x)^2
C
C     and that the value of  rho that balances f(x) and
C
C     \sum_i P(lambda_i,c_i(x),rho) is given by
C
C     rho = f(x) / ( 0.5 \sum_i c_i(x)^2 ).

C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer i
      double precision rhoini,sumc

      sumc = 0.0d0
      do i = 1,m
          if ( equatn(i) .or. c(i) .gt. 0.0d0 ) then
              sumc = sumc + 0.5d0 * c(i) ** 2
          end if
      end do

      rhoini = 10.0d0 * max( 1.0d0, abs( f ) ) / max( 1.0d0, sumc )

      rhoini = max( macheps12, min( rhoini, 1.0d0 / macheps12 ) )

      do i = 1,m
          rho(i) = rhoini
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision cnorm,cnormu
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

C     Subroutine sevalnal is always called by GENCAN with the same point
C     before calling this subroutine. Then, constraints are computed and
C     saved.

C     Compute infeasibility

      cnorm = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( c(i) ) )
          else
              cnorm = max( cnorm, c(i) )
          end if
      end do

      if ( scale ) then
          cnormu = 0.0d0
          do i = 1,m
              if ( equatn(i) ) then
                  cnormu = max( cnormu, abs( sc(i) * c(i) ) )
              else
                  cnormu = max( cnormu, sc(i) * c(i) )
              end if
          end do

      else
          cnormu = cnorm
      end if

      end
C     ******************************************************************
C     ******************************************************************

      subroutine gencan(n,x,l,u,m,lambda,equatn,linear,rho,epsfeas,
     +epsopt,maxit,iter,f,g,gpsupn,cnorm,cnormu,geninfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer geninfo,inform,iter,m,maxit,n
      double precision cnorm,cnormu,epsfeas,epsopt,f,gpsupn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision g(n),l(n),lambda(m),rho(m),u(n),x(n)

C     Solves the box-constrained minimization problem
C
C            Minimize f(x) subject to l <= x <= u
C     
C     using the method described in 
C
C     E. G. Birgin and J. M. Martinez, ''Large-scale active-set box-
C     constrained optimization method with spectral projected 
C     gradients'', Computational Optimization and Applications 23, pp. 
C     101-125, 2002.  

C     geninfo:
C
C     0: Small continuous-projected-gradient norm
C     1: Maximum number of gencan iterations reached
C     3: Lack of progress in the objective function value
C     4: Lack of progress in the continuous-projected-gradient norm
C     5: Lack of progress in the current point
C     6: Unbounded objective function?
C     7: Too small backtracking step. Wrong gradient?

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/
C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/
C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/

C     LOCAL SCALARS
      character * 6 precond
      character * 2 innittype
      logical forceoi,memfail,vustop
      integer cgcnt,cginfo,cgiter,cgmaxit,chcnt,cinit,i,innit,iscnt,
     +        itnfp,itngp,itnxp,lsinfo,mfit,nind,nindprev,nwcnt,lvfit,
     +        rbdnnz,trcnt,trinfo,xind,cisit
      double precision acgeps,amax,bcgeps,cgeps,cgdel,dsupn,fdiff,fplus,
     +        gieucn,gisupn,gpeucn,gpi,lamspg,mslamb,newtrdel,ssupn,tmp,
     +        trdel,tsmall,ysupn,xeucn,xsupn

C     LOCAL ARRAYS
      integer rbdind(nmax)
      character rbdtype(nmax)
      double precision d(nmax),gplus(nmax),xplus(nmax)

C     DATA BLOCKS
      character * 19 ittext(0:5)
      integer BDSREACHED,FSTORDPNT,SECORDPNT,SMALLSTEP,SMALLTRRAD,
     +        UNBOUNDED,UNDEFSTEP

      data ittext(0) /'(Initial point)   '/
      data ittext(1) /'(SPG step)        '/
      data ittext(2) /'(Truncated Newton)'/
      data ittext(3) /'(MS system)       '/
      data ittext(4) /'(Trust region)    '/
      data ittext(5) /'(Inner SPG)       '/

      data BDSREACHED /2/
      data FSTORDPNT  /6/
      data SECORDPNT  /7/
      data SMALLSTEP  /3/
      data SMALLTRRAD /4/
      data UNBOUNDED  /2/
      data UNDEFSTEP  /5/

C     EXTERNAL
      logical sstop

C     EXTERNAL SUBROUTINES
      external sevalal

C     ==================================================================
C     Initialization
C     ==================================================================

C     Set some initial values:

C     to record a memory failure in the direct solver
      memfail = .false.

C     for Conjugate Gradients
      precond = 'QNCGNA'

C     for the first inner-to-the-face minimization algorithm
      innittype = 'CG'

C     for testing lack of progress checking f, g and x
      itnfp = 0
      itngp = 0
      itnxp = 0

C     to force a leaving-face iteration when lack of progress is 
C     detected within a face
      forceoi = .false.

C     for calculating More-Sorensen's direction
      mslamb =  0.0d0

C     for counting number of iterations as well as inner-to-the-face
C     and leaving-face iterations
      iter  = 0
      cgcnt = 0
      nwcnt = 0
      trcnt = 0
      iscnt = 0
      chcnt = 0
      lvfit = 0
      innit = 0
      cinit = 0
      cisit = 0

C     just to print "Initial point" in the first ouput
      ittype = 0

C     Print problem information

      if ( iprintinn .ge. 1 ) then
          write(* ,1000) n
          write(10,1000) n
      end if

      if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
          write(* ,1010) nprint,(l(i),i=1,nprint)
          write(* ,1020) nprint,(u(i),i=1,nprint)
          write(* ,1030) nprint,(x(i),i=1,nprint)

          write(10,1010) nprint,(l(i),i=1,nprint)
          write(10,1020) nprint,(u(i),i=1,nprint)
          write(10,1030) nprint,(x(i),i=1,nprint)
      end if

C     Project initial guess. If the initial guess is infeasible, 
C     projection puts it into the box.

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute function and gradient at the initial point

      call sevalal(n,x,m,lambda,rho,equatn,linear,f,inform)
      if ( inform .lt. 0 ) return

      call sevalnal(n,x,m,lambda,rho,equatn,linear,g,inform)
      if ( inform .lt. 0 ) return

C     Compute x Euclidian and sup norms, continuous-project-gradient 
C     Euclidian and sup norms, internal gradient Euclidian norm. Set
C     nind as the number of free variables and save in array ind their 
C     identifiers. Also set variable sameface (same face) indicating 
C     that the "previous iterate" does not belong to the current face.

      sameface = .false.

      nind   = 0
      xsupn  = 0.0d0
      xeucn  = 0.0d0
      gpsupn = 0.0d0
      gpeucn = 0.0d0
      gisupn = 0.0d0
      gieucn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind = i
          end if
          xeucn = xeucn + x(i) ** 2
          gpi = x(i) - g(i)
          if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
              gpi = - g(i)
          else
              gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
          end if
          gpsupn = max( gpsupn, abs( gpi ) )
          gpeucn = gpeucn + gpi ** 2
          if ( x(i) - l(i) .gt. macheps23 * max( 1.0d0, l(i) ) .and. 
     +         u(i) - x(i) .gt. macheps23 * max( 1.0d0, u(i) ) ) then
              gisupn    = max( gisupn, abs( gpi ) )
              gieucn    = gieucn + gpi ** 2
              nind      = nind + 1
              ind(nind) = i
          end if
      end do
      xeucn  = sqrt( xeucn  )
      gpeucn = sqrt( gpeucn )
      gieucn = sqrt( gieucn )

C     Compute trut-region radius, in case betra is used as inner solver

      trdel    = max( trdelmin, trdelini * max( 1.0d0, xeucn ) )
      newtrdel = 0.0d0
	
C     To be used in feasibility problems, compute constraints norm

      call sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)
      if ( inform .lt. 0 ) return

C     Initial spectral steplength

C     Compute a small step and set the point at which the auxiliary
C     gradient will be computed

      if ( gpsupn .ne. 0.0d0 ) then 
          tsmall = macheps12 * max( 1.0d0, xsupn / gpsupn )
      else
          tsmall = 0.0d0
      end if

      do i = 1,n
          gpi = x(i) - g(i)
          if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
              gpi = - g(i)
          else
              gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
          end if
          s(i) = x(i) + tsmall * gpi
      end do

C     Compute the gradient at the auxiliary point

      call ievalnalu(n,s,m,lambda,rho,equatn,linear,.false.,y,inform)
      if ( inform .lt. 0 ) return

C     Compute s = x_{1/2} - x_0 and y = g_{1/2} - g_0

      sts = 0.0d0
      sty = 0.0d0
      ssupn = 0.0d0
      ysupn = 0.0d0
      yeucn = 0.0d0
      do i = 1,n
          s(i) = s(i) - x(i)
          y(i) = y(i) - g(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
          ssupn = max( ssupn, abs( s(i) ) )
          ysupn = max( ysupn, abs( y(i) ) )
          yeucn = yeucn + y(i) ** 2
      end do
      seucn = sqrt( sts )
      yeucn = sqrt( yeucn )

C     Compute a linear relation between gpsupn and cgeps, i.e.,
C     scalars a and b such that 
c
C         a * log10( ||g_P(x_ini)|| ) + b = log10(cgeps_ini) and
c
C         a * log10( ||g_P(x_fin)|| ) + b = log10(cgeps_fin),
c
C     where cgeps_ini and cgeps_fin are provided. Note that if 
C     cgeps_ini is equal to cgeps_fin then cgeps will be always 
C     equal to cgeps_ini and cgeps_fin.

      acgeps = log10( cgepsf / cgepsi ) / log10( cggpnf / gpsupn )
      bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )

C     Print initial information

      if ( iprintinn .eq. 1 ) then
          if ( mod(iter,10) .eq. 0 ) then
              write(* ,1040)
              write(10,1040)
          end if
          write(* ,1050) iter,f,gpsupn,xsupn
          write(10,1050) iter,f,gpsupn,xsupn

      else if ( iprintinn .ge. 2 ) then
          write(* ,1070) iter,ittext(ittype)
          write(10,1070) iter,ittext(ittype)

          if ( iprintinn .ge. 3 ) then
              write(* ,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
              write(10,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
          end if

          write(* ,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn
          write(10,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn

          if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
              write(* ,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(* ,1110) nprint,(x(i),i=1,nprint)
              write(* ,1130) nprint,(g(i),i=1,nprint)
              write(* ,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
              write(10,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(10,1110) nprint,(x(i),i=1,nprint)
              write(10,1130) nprint,(g(i),i=1,nprint)
              write(10,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
          end if
      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          open(20,file='gencan-tabline.out')
          write(20,1400) f,0.0d0,f,0.0d0,gpsupn,f,0.0d0,f,0.0d0,gpsupn,
     +                   inform,9,n,m,0,iter,fcnt,0,0,0,0,999.9d0
          close(20)
      end if

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     Test user-provided stopping criterion

      if ( useustp ) then

          vustop = sstop(n,x,m,lambda,rho,equatn,linear,inform)
          if ( inform .lt. 0 ) return

          if ( vustop ) then
              geninfo = 8

              if ( iprintinn .ge. 1 ) then
                  write(*, 1280)
                  write(10,1280)
              end if

              return
          end if

      end if

C     Test whether the continuous-projected-gradient sup-norm
C     is small enough to declare convergence

      if ( .not. ignoref ) then
          if ( gpsupn .le. epsopt ) then
              geninfo = 0

              if ( iprintinn .ge. 1 ) then
                  write(*, 1200)
                  write(10,1200)
              end if

              return
          end if

      else
          if ( cnormu .le. epsfeas ) then
              geninfo = 0

              if ( iprintinn .ge. 1 ) then
                  write(*, 1200)
                  write(10,1200)
              end if

              return
          end if
      end if

C     Test whether the number of iterations is exhausted

      if ( iter .ge. maxit ) then
          geninfo = 1

          if ( iprintinn .ge. 1 ) then
              write(*, 1210)
              write(10,1210)
          end if

          return
      end if

C     ==================================================================
C     Test stopping criteria related to lack of progress
C     ==================================================================

      if ( iter .gt. 0 ) then

C     Test whether we have performed many iterations without 
C     moving from the current point, by checking the functional value

      if ( abs( fdiff ) .le. macheps * max( 1.0d0, abs( f ) ) ) then

          itnfp = itnfp + 1

          if ( itnfp .ge. maxinnitnp ) then

              if ( itnfp .gt. cinit ) then
                  geninfo = 3

                  if ( iprintinn .ge. 1 ) then
                      write(*, 1230)
                      write(10,1230)
                  end if

                  return

              else
                  forceoi = .true.
              end if
          end if

      else
          itnfp = 0
      end if

C     Test whether we have performed many iterations without 
C     moving from the current point, by checking their gradients

      if ( ysupn .le. macheps * max( 1.0d0, gpsupn ) .or.
     +     yeucn .le. macheps * max( 1.0d0, gpeucn ) ) then

          itngp = itngp + 1

          if ( itngp .ge. maxinnitnp ) then

              if ( itngp .gt. cinit ) then
                  geninfo = 4

                  if ( iprintinn .ge. 1 ) then
                      write(*, 1240)
                      write(10,1240)
                  end if

                  return

              else
                  forceoi = .true.
              end if
          end if

      else
          itngp = 0
      end if

C     Test whether we have performed many iterations without 
C     moving from the current point, by checking the step norm

      if ( ssupn .le. macheps * max( 1.0d0, xsupn ) .or.
     +     seucn .le. macheps * max( 1.0d0, xeucn ) ) then

          itnxp = itnxp + 1

          if ( itnxp .ge. maxinnitnp ) then

              if ( itnxp .gt. cinit ) then

                  geninfo = 5

                  if ( iprintinn .ge. 1 ) then
                      write(*, 1250)
                      write(10,1250)
                  end if

                  return

              else
                  forceoi = .true.
              end if
          end if

      else
          itnxp = 0
      end if

      end if

C     ==================================================================
C     Iteration
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Compute new iterate
C     ==================================================================

C     We abandon the current face if the norm of the internal gradient
C     (here, internal components of the continuous projected gradient)
C     is smaller than eta times the norm of the continuous 
C     projected gradient. Using eta = 0.1 is a rather conservative 
C     strategy in the sense that internal iterations are preferred over 
C     SPG iterations. Replace eta = 0.1 by other tolerance in (0,1) if 
C     you find it convenient. 

      if ( gieucn .le. eta * gpeucn .or. 
     +     gisupn .le. eta * gpsupn .or. forceoi ) then

C         ==============================================================
C         Some constraints should be abandoned. Compute the new iterate 
C         using an SPG iteration
C         ==============================================================

          forceoi = .false.
          lvfit = lvfit + 1
          cinit = 0

C         Compute spectral steplength

          if ( sty .le. 0.0d0 ) then
              lamspg = max( 1.0d0, xsupn / gpsupn )
          else
              lamspg = sts / sty
          end if
          lamspg = min( lspgma, max( lspgmi, lamspg ) )

C         Perform safeguarded quadratic interpolation along the 
C         spectral continuous projected gradient

          call spgls(n,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg,
     +    xplus,fplus,tmp,d,sevalal,lsinfo,inform) 
          if ( inform .lt. 0 ) return

          call sevalnal(n,xplus,m,lambda,rho,equatn,linear,gplus,inform)
          if ( inform .lt. 0 ) return

C         Set iteration type

          ittype = 1

      else

C         ==============================================================
C         The new iterate will belong to the closure of the current face
C         ==============================================================

          innit = innit + 1
          cinit = cinit + 1

C         Shrink the point, its gradient and the bounds

          call shrink(nind,x)
          call shrink(nind,g)
          call shrink(nind,l)
          call shrink(nind,u)

C         Save values of fixed variables for further evaluations

          do i = 1,n - nind
              xcomplement(i) = x(nind + i)
          end do

          nt = n

C         Compute the new point using a trust-region approach

          if ( innittype .eq. 'TR' ) then

              call betra(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,
     +        trdel,newtrdel,mslamb,epsopt,xeucn,xsupn,gpsupn,d,rbdind,
     +        rbdtype,xplus,fplus,gplus,trcnt,iscnt,chcnt,memfail,
     +        trinfo,inform)

              if ( inform .lt. 0 ) return

              if ( memfail ) then
                  mfit = iter
              end if
          
          end if

C         Compute the descent direction solving the newtonian system

          if ( innittype .eq. 'DS' ) then

C             Solve the Newtonian system
C
C                 ( H + rho A^T A ) x = b
C
C             by solving the Martinez-Santos system
C
C                 H x + A^t y = b        
C                 A x - y/rho = 0
C
C             with a direct solver.

              call newtd(nind,x,l,u,g,m,rho,equatn,d,memfail,inform)
              if ( inform .lt. 0 ) return

              nwcnt = nwcnt + 1

              if ( memfail ) then
                  mfit = iter
              end if

C             Set iteration type

              ittype = 3

          end if

          if ( innittype .eq. 'CG' .or. 
     +       ( innittype .eq. 'DS' .and. memfail ) .or.
     +       ( innittype .eq. 'TR' .and. memfail ) .or.
     +       ( innittype .eq. 'TR' .and. trinfo .eq. UNDEFSTEP ) ) then

C             Compute "trust-region radius"

              if ( iter .eq. 1 ) then
                  cgdel = max( 100.0d0, 100.0d0 * xsupn ) 
              else
                  cgdel = max( delmin, 10.0d0 * ssupn )
              end if

C             Set conjugate gradient stopping criteria. 

C             cgmaxit = min( 2 * nind, 10000 )
              cgmaxit = min( nind, 10000 )
              cgeps   = 10.0d0 ** ( acgeps * log10( gpsupn ) + bcgeps )
              cgeps   = max( cgepsf, min( cgepsi, cgeps ) )

C             Call Conjugate Gradients to solve the Newtonian system

              call cgm(nind,x,m,lambda,rho,equatn,linear,l,u,g,cgdel,
     +        cgeps,cgmaxit,precond,d,rbdnnz,rbdind,rbdtype,cgiter,
     +        cginfo,inform)

              cgcnt = cgcnt + cgiter

              if ( inform .lt. 0 ) return

C             Set iteration type

              ittype = 2

          end if

          if ( ittype .le. 3 ) then

C             Compute maximum step

              if ( ittype .eq. 2 .and. cginfo .eq. BDSREACHED ) then
                  amax = 1.0d0
                  
              else
                  call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)
              end if

C             Perform line search

              call tnls(nind,x,l,u,m,lambda,rho,equatn,linear,f,g,amax,
     +        d,rbdnnz,rbdind,rbdtype,xplus,fplus,gplus,lsinfo,inform)
              if ( inform .lt. 0 ) return

          end if
              
C         Set iteration type based on the lenght of the maximum step

          dsupn = 0.0d0 
          do i = 1,nind
              dsupn = max( dsupn, abs( d(i) ) )
          end do
          
          if ( avoidds .or. 
     +         ( memfail .and. iter - mfit .le. 3 ) ) then
              innittype = 'CG'
              
          else
              if ( innittype .eq. 'DS' .and. 
     +             amax .le. macheps12 * dsupn ) then
                  innittype = 'CG'
              elseif ( cisit .ge. 3 ) then
                  innittype = 'CG'
              else
                  innittype = 'DS'
              end if
          end if
          
C         Expand the point, its gradient and the bounds. Also expand 
C         xplus and gplus.

          call expand(nind,x)
          call expand(nind,g)
          call expand(nind,l)
          call expand(nind,u)

          call expand(nind,xplus)
          call expand(nind,gplus)
      end if

C     Compute s = xplus - x and y = gplus - g

      fdiff = fplus - f

      sts = 0.0d0
      sty = 0.0d0
      ssupn = 0.0d0
      ysupn = 0.0d0
      yeucn = 0.0d0
      do i = 1,n
          s(i) = xplus(i) - x(i)
          y(i) = gplus(i) - g(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
          ssupn = max( ssupn, abs( s(i) ) )
          ysupn = max( ysupn, abs( y(i) ) )
          yeucn = yeucn + y(i) ** 2
      end do
      seucn = sqrt( sts )
      yeucn = sqrt( yeucn )

C     Set new point

      f = fplus

      do i = 1,n
          x(i) = xplus(i)
          g(i) = gplus(i)
      end do

C     ==================================================================
C     Prepare for the next iteration 
C     ==================================================================

C     This adjustment/projection is ''por lo que las putas pudiera''

c      do i = 1,n
c          if ( x(i) - macheps23 * max( 1.0d0, abs( x(i) ) ) .le. 
c     +         l(i) ) then
c              x(i) = l(i)
c          else if ( x(i) + macheps23 * max( 1.0d0, abs( x(i) ) ) .ge. 
c     +              u(i) ) then
c              x(i) = u(i)
c          end if
c      end do

C     Compute continuous-project-gradient Euclidian and Sup norms,
C     internal gradient Euclidian norm, and store in nind the number of
C     free variables and in array ind their identifiers. Verify whether 
C     x and xprev belong to the same face.

      sameface = .true.
      nindprev = nind

      nind   = 0
      xsupn  = 0.0d0
      xeucn  = 0.0d0
      gpsupn = 0.0d0
      gpeucn = 0.0d0
      gisupn = 0.0d0
      gieucn = 0.0d0
      do i = 1,n
          if ( abs( x(i) ) .ge. xsupn ) then
              xsupn = abs( x(i) )
              xind = i
          end if
          xeucn = xeucn + x(i) ** 2
          gpi = x(i) - g(i)
          if ( l(i) .le. gpi .and. gpi .le. u(i) ) then
              gpi = - g(i)
          else
              gpi = max( l(i), min( gpi, u(i) ) ) - x(i)
          end if
          gpsupn = max( gpsupn, abs( gpi ) )
          gpeucn = gpeucn + gpi ** 2
          if ( x(i) - l(i) .gt. macheps23 * max( 1.0d0, l(i) ) .and. 
     +         u(i) - x(i) .gt. macheps23 * max( 1.0d0, u(i) ) ) then
              gisupn    = max( gisupn, abs( gpi ) )
              gieucn    = gieucn + gpi ** 2
              nind      = nind + 1
              if ( nind .gt. nindprev .or. ind(nind) .ne. i ) then
                  sameface = .false.
              end if
              ind(nind) = i
          end if
      end do
      xeucn  = sqrt( xeucn  )
      gpeucn = sqrt( gpeucn )
      gieucn = sqrt( gieucn )

      if ( ittype .eq. 5 .and. sameface ) then
          cisit = cisit + 1
      else
          cisit = 0
      end if

C     To be used in feasibility problems, compute constraints norm

      call sevalfeas(n,x,m,equatn,cnorm,cnormu,inform)
      if ( inform .lt. 0 ) return

C     Print information of this iteration

      if ( iprintinn .eq. 1 ) then
          if ( mod(iter,10) .eq. 0 ) then
              write(* ,1040)
              write(10,1040)
          end if
          write(* ,1060) iter,f,gpsupn,xsupn,ssupn,nind,lvfit,innit,fcnt
          write(10,1060) iter,f,gpsupn,xsupn,ssupn,nind,lvfit,innit,fcnt

      else if ( iprintinn .ge. 2 ) then
          write(* ,1070) iter,ittext(ittype)
          write(10,1070) iter,ittext(ittype)

          if ( iprintinn .ge. 3 ) then
              write(* ,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
              write(10,1080) lvfit,innit,fcnt,cgcnt,nwcnt,trcnt,iscnt,
     +                       chcnt
          end if

          write(* ,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn
          write(10,1090) f,gpsupn,gpeucn,gisupn,gieucn,xind,xsupn,xeucn,
     +                   ssupn,seucn

          if ( iprintinn .ge. 4 .and. nprint .ne. 0 ) then
              write(* ,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(* ,1110) nprint,(x(i),i=1,nprint)
              write(* ,1130) nprint,(g(i),i=1,nprint)
              write(* ,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
              write(10,1100) min(nprint,nind),nind,
     +                       (ind(i),i=1,min(nprint,nind))
              write(10,1110) nprint,(x(i),i=1,nprint)
              write(10,1130) nprint,(g(i),i=1,nprint)
              write(10,1140) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),
     +                       i=1,nprint)
          end if
      end if

C     Save intermediate data for crash report

      if ( iprintctl(6) ) then
          open(20,file='gencan-tabline.out')
          write(20,1400) f,0.0d0,f,0.0d0,gpsupn,f,0.0d0,f,0.0d0,gpsupn,
     +                   inform,9,n,m,0,iter,fcnt,0,0,0,0,999.9d0
          close(20)
      end if

C     ==================================================================
C     Test line-search and trust-region related stopping criteria
C     ==================================================================

      if ( ittype .ge. 4 ) then
          
C         Test whether the functional value is unbounded

          if ( trinfo .eq. UNBOUNDED ) then
              geninfo = 6
              
              if ( iprintinn .ge. 1 ) then
                  write(*, 1260)
                  write(10,1260)
              end if
              
              return
          end if

C         Test whether the line search did a very small step
          
          if ( trinfo .eq. SMALLSTEP ) then
              forceoi = .true.
          end if

C         Test whether the trust-region radius became too small
          
          if ( trinfo .eq. SMALLTRRAD ) then
              forceoi = .true.
          end if

C         Test whether the point is stationary in the face
          
          if ( trinfo .eq. FSTORDPNT ) then
              forceoi = .true.
          end if

          if ( trinfo .eq. SECORDPNT ) then
              forceoi = .true.
          end if

      else

C         Test whether the functional value is unbounded

          if ( lsinfo .eq. UNBOUNDED ) then
              geninfo = 6
              
              if ( iprintinn .ge. 1 ) then
                  write(*, 1260)
                  write(10,1260)
              end if
              
              return
          end if
          
C         Test whether the line search did a very small step
          
          if ( lsinfo .eq. SMALLSTEP ) then
              if ( ittype .eq. 1 ) then
                  geninfo = 7
                  
                  if ( iprintinn .ge. 1 ) then
                      write(*, 1270)
                      write(10,1270)
                  end if
                  
                  return
              else
                  forceoi = .true.
              end if
          end if

      end if

C     ==================================================================
C     Iterate 
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,5X,'Entry to GENCAN.',
     +       /,5X,'Number of variables: ',I7)

 1010 format(/,5X,'Lower bounds  (first ',I7,' components):',
     +       /,(5X, 6(1X,1P,D11.4)))
 1020 format(/,5X,'Upper bounds  (first ',I7,' components):',
     +       /,(5X, 6(1X,1P,D11.4)))
 1030 format(/,5X,'Initial point (first ',I7,' components):',
     +       /,(5X, 6(1X,1P,D11.4)))

 1040 format(/,5X,'It',16X,'f',2X,'gpsupn',3X,'xsupn',3X,'ssupn',4X,
     +            'nind',3X,'lvfit',3X,'innit',4X,'fcnt')
 1050 format(     I7,1X,1P,D16.8,2(1X,1P,D7.1))
 1060 format(     I7,1X,1P,D16.8,3(1X,1P,D7.1),4(1X,I7))

 1070 format(/,5X,'GENCAN ITERATION                  ',12X,'= ',I7,1X,
     +         A19)
 1080 format(/,5X,'Leaving-face iterations           ',12X,'= ',I7,
     +       /,5X,'Inner-to-the-face iterations      ',12X,'= ',I7,
     +       /,5X,'Functional evaluations            ',12X,'= ',I7,
     +       /,5X,'Conjugate gradient iterations     ',12X,'= ',I7,
     +       /,5X,'Newtonian system factorizations   ',12X,'= ',I7,
     +       /,5X,'Trust-region iterations           ',12X,'= ',I7,
     +       /,5X,'Inner SPG iterations              ',12X,'= ',I7,
     +       /,5X,'Trust-region matrix factorizations',12X,'= ',I7)

 1090 format(/,5X,'Functional value                              = ',
     +             1P,D24.16,
     +       /,5X,'Sup-norm of the continuous projected gradient = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1,
     +       /,5X,'Sup-norm of the internal projection of gp     = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1,
     +       /,5X,'Sup-norm of x (attained at x_{', I7,'})       = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1,
     +       /,5X,'Sup-norm of x - x_{prev}                      = ',
     +             1P,D7.1,' 2-norm = ',1P,D7.1)
 1100 format(/,5X,'Current free variables (first ',I7,', total ',
     +            'number ',I7,'): ',
     +       /,(5X, 6(1X,I7)))
 1110 format(/,5X,'Current point (first ',I7, ' components): ',
     +       /,(5X, 6(1X,1P,D11.4)))
 1130 format(/,5X,'Current gradient (first ',I7,' components): ',
     +       /,(5X, 6(1X,1P,D11.4)))
 1140 format(/,5X,'Current continuous projected gradient (first ',I7, 
     +         5X,'components): ',
     +       /,(5X, 6(1X,1P,D11.4)))

 1200 format(/,5X,'Flag of GENCAN: Solution was found.',/)
 1210 format(/,5X,'Flag of GENCAN: Maximum of iterations reached.',/)
 1230 format(/,5X,'Flag of GENCAN: Lack of progress in the ',
     +            'functional value.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1240 format(/,5X,'Flag of GENCAN: Lack of progress in the gradient ',
     +            'norm.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1250 format(/,5X,'Flag of GENCAN: Lack of progress in the current ',
     +            'point.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1260 format(/,5X,'Flag of GENCAN: Objective function seems to be ',
     +            'unbounded.',/)
 1270 format(/,5X,'Flag of GENCAN: Too small step in the line search.',
     +       /,5X,'Probably, an exaggerated small norm of the ',
     +            'continuous projected gradient',
     +       /,5X,'is being required for declaring convergence.',/)
 1280 format(/,5X,'Flag of GENCAN: User-provided stopping criterion ',
     +            'satisfied.',/)

 1400 format(1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,1P,D24.16,1X,1P,D7.1,1X,1P,D24.16,1X,1P,D7.1,1X,1P,D7.1,
     +       1X,I3,1X,I1,1X,I6,1X,I6,1X,I2,1X,I7,1X,I7,1X,I2,1X,I7,
     +       1X,I7,1X,I7,0P,F8.2)

      end


C     ******************************************************************
C     ******************************************************************

      subroutine compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,rbdnnz
      double precision amax

C     ARRAY ARGUMENTS
      integer rbdind(nind)
      character rbdtype(nind)
      double precision d(nind),l(nind),u(nind),x(nind)

C     Compute maximum step amax > 0 such that l <= x + amax d <= u.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer i
      double precision amaxi

      rbdnnz = 0

      amax = bignum

      do i = 1,nind
         if ( d(i) .gt. 0.0d0 ) then

             amaxi = ( u(i) - x(i) ) / d(i)

             if ( amaxi .lt. amax ) then
                 amax       = amaxi
                 rbdnnz     = 1
                 rbdind(1)  = i
                 rbdtype(1) = 'U'
             else if ( amaxi .eq. amax ) then
                 rbdnnz          = rbdnnz + 1
                 rbdind(rbdnnz)  = i
                 rbdtype(rbdnnz) = 'U'
             end if

         else if ( d(i) .lt. 0.0d0 ) then

             amaxi = ( l(i) - x(i) ) / d(i)

             if ( amaxi .lt. amax ) then
                 amax      = amaxi
                 rbdnnz     = 1
                 rbdind(1)  = i
                 rbdtype(1) = 'L'
             else if ( amaxi .eq. amax ) then
                 rbdnnz          = rbdnnz + 1
                 rbdind(rbdnnz)  = i
                 rbdtype(rbdnnz) = 'L'
             end if

         end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine newtd(nind,x,l,u,g,m,rho,equatn,d,memfail,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical memfail
      integer inform,nind,m

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision d(nind),g(*),l(*),rho(m),u(*),x(*)

C     This subroutine solves the Newtonian system
C
C             ( H + rho A^T A ) x = b
C
C     by solving the Martinez-Santos system
C
C             H x + A^t y = b
C             A x - y/rho = 0

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/

C     PARAMETERS
      integer dimmax
      parameter ( dimmax = nmax + mmax )

C     COMMON ARRAYS
      double precision mindiag(nmax)

C     LOCAL SCALARS
      integer dim,hnnz,i,iter,lssinfo,nneigv,pind
      double precision adsupn,diff,pval

C     LOCAL ARRAYS
      integer hdiag(dimmax),hlin(hnnzmax),hcol(hnnzmax)
      double precision adddiag(dimmax),hval(hnnzmax),sol(dimmax)

C     COMMON BLOCKS
      common /diadat/ mindiag
      save   /diadat/

C     ------------------------------------------------------------------
C     Presentation
C     ------------------------------------------------------------------

      if ( iprintinn .ge. 5 ) then
          write(* ,1000)
          write(10,1000)
      end if

C     ------------------------------------------------------------------
C     Initialization
C     ------------------------------------------------------------------

      iter    = 0
      memfail = .false.

      call lssini(sclsys,.true.,.false.)

C     ------------------------------------------------------------------
C     Compute ML matrix
C     ------------------------------------------------------------------

      call mlsyst(nind,x,g,m,rho,equatn,hlin,hcol,hval,hnnz,hdiag,sol,
     +dim,inform)

      if ( inform .lt. 0 ) return

C     ------------------------------------------------------------------
C     Analyse sparsity pattern
C     ------------------------------------------------------------------

      call lssana(dim,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      end if

C     ------------------------------------------------------------------
C     Main loop
C     ------------------------------------------------------------------

 100  continue

      iter = iter + 1

C     ------------------------------------------------------------------
C     Compute regularization
C     ------------------------------------------------------------------

      if ( iter .eq. 1 ) then

          if ( sameface .and. ittype .eq. 3 ) then

              do i = 1,nind
                  mindiag(i) = 0.1d0 * mindiag(i)
              end do

          else

              do i = 1,nind
                  if ( g(i) .eq. 0.0d0 ) then
                      mindiag(i) = 0.0d0
                  else
                      if ( g(i) .gt. 0.0d0 ) then
                          diff = x(i) - l(i)
                      else if ( g(i) .lt. 0.0d0 ) then
                          diff = u(i) - x(i)
                      end if
                      mindiag(i) = abs( g(i) / diff )
                  end if
              end do

          end if

      else

          do i = 1,nind
              if ( mindiag(i) .eq. 0.0d0 ) then
                  mindiag(i) = macheps23
              else
                  mindiag(i) = 10.0d0 * mindiag(i)
              end if
          end do

      end if

      do i = 1,nind
          adddiag(i) = max( macheps23, mindiag(i) - hval(hdiag(i)) )
      end do

      do i = nind + 1,dim
          adddiag(i) = 0.0d0
      end do

      adsupn = 0.0d0
      do i = 1,dim
          adsupn = max( adsupn, adddiag(i) )
      end do

      if ( iprintinn .ge. 5 ) then
          write(* ,1010) adsupn
          write(10,1010) adsupn
      end if

C     ------------------------------------------------------------------
C     Factorize matrix
C     ------------------------------------------------------------------

      call lssfac(dim,hnnz,hlin,hcol,hval,hdiag,adddiag,pind,pval,
     +nneigv,lssinfo)

      if ( lssinfo .eq. 0 .or. lssinfo .eq. 1 ) then

          if ( nneigv .ne. dim - nind ) then
C           ! WRONG INERTIA (SEE NOCEDAL AND WRIGHT)

C             Lemma 16.3 [pg. 447]: Assume that the Jacobian of the
C             constraints has full rank and that the reduced Hessian
C             Z^T H Z is positive definite. Then the Jacobian of the
C             KKT system has n positive eigenvalues, m negative
C             eigenvalues, and no zero eigenvalues.

C             Note that at this point we know that the matrix has no
C             zero eigenvalues. nneigv gives the number of negative
C             eigenvalues.

              if ( iprintinn .ge. 5 ) then
                  write(* ,1020) nneigv,dim - nind
                  write(10,1020) nneigv,dim - nind
              end if

              go to 100

          else

              if ( iprintinn .ge. 5 ) then
                  write(* ,1030)
                  write(10,1030)
              end if

          end if

      else if ( lssinfo .eq. 2 ) then
        ! SINGULAR JACOBIAN

          if ( iprintinn .ge. 5 ) then
              write(* ,1040)
              write(10,1040)
          end if

          go to 100

      else if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else ! if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if

C     ------------------------------------------------------------------
C     Solve
C     ------------------------------------------------------------------

      call lsssol(dim,sol)

      do i = 1,nind
          d(i) = sol(i)
      end do

      if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
          write(*, 1050) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
          write(10,1050) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
      end if

C     NON-EXECUTABLE STATEMENTS

 1000 format(/,5X,'Sparse factorization of the ML system.')
 1010 format(  5X,'Maximum value added to the diagonal: ',1P,D24.16)
 1020 format(  5X,'ML-matrix with wrong inertia.',
     +       /,5X,'Actual number of negative eigenvalues  = ',I16,'.',
     +       /,5X,'Desired number of negative eigenvalues = ',I16,'.')
 1030 format(  5X,'Direct solver finished successfully.')
 1040 format(  5X,'ML-matrix numerically singular.')
 1050 format(/,5X,'Newton direction (first ',I7,' components): ',
     +       /,1(5X,6(1X,1P,D11.4)))

      end

C     ******************************************************************
C     ******************************************************************

      subroutine mlsyst(nind,x,nal,m,rho,equatn,ulin,ucol,uval,unnz,
     +udiag,b,dim,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer dim,inform,m,nind,unnz

C     ARRAY ARGUMENTS
      logical equatn(m)
      integer ucol(*),udiag(*),ulin(*)
      double precision b(*),nal(nind),rho(m),uval(*),x(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer col,i,j,k,lin,var

C     LOCAL ARRAYS
      integer wi(nmax)

C     This subrotuine is called from the reduced space.

C     MATRIX

C     Compute Hessian of the Lagrangian

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

      call expand(nind,x)

      call sevalhl(nt,x,m,dpdc,ulin,ucol,uval,unnz,inform)
      if ( inform .lt. 0 ) return

      call shrink(nind,x)

C     Preparation for shrink (wi indicates, for each free variable x_i,
C     its rank within the set of free variables. wi(i)=0 if x_i is not
C     a free variable)

      do i = 1,nt
         wi(i) = 0
      end do

      do i = 1,nind
         wi(ind(i)) = i
      end do

C     Shrink Hessian of the Lagrangian and set diagonal-elements indices

      k = 0

      do i = 1,nind
          udiag(i) = 0
      end do

      do i = 1,unnz
          lin = wi(ulin(i))
          col = wi(ucol(i))

          if ( lin .ne. 0 .and. col .ne. 0 ) then
              k = k + 1
              ulin(k) = lin
              ucol(k) = col
              uval(k) = uval(i)

              if ( lin .eq. col ) then
                  udiag(lin) = k
              end if
          end if
      end do

      do i = 1,nind
          if ( udiag(i) .eq. 0 ) then
              k = k + 1
              ulin(k) = i
              ucol(k) = i
              uval(k) = 0.0d0

              udiag(i) = k
          end if
      end do

C     Shrink Jacobian and add diagonal matrix - 1.0 / rho

      dim = nind

      do j = 1,m
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

              dim = dim + 1

              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  var = wi(jcvar(i))

                  if ( var .ne. 0 ) then
                      k = k + 1
                      ulin(k) = dim
                      ucol(k) = var
                      uval(k) = jcval(i)
                  end if
              end do

              k = k + 1
              ulin(k) = dim
              ucol(k) = dim
              uval(k) = - 1.0d0 / rho(j)

              udiag(dim) = k
          end if
      end do

      unnz = k

C     RHS

      do i = 1,nind
          b(i) = - nal(i)
      end do

      do i = nind + 1,dim
          b(i) = 0.0d0
      end do

      end
C     ******************************************************************
C     ******************************************************************

      subroutine cgm(nind,x,m,lambda,rho,equatn,linear,l,u,g,delta,
     +eps,maxit,precond,d,rbdnnz,rbdind,rbdtype,iter,cginfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      character * 6 precond
      integer cginfo,inform,iter,m,maxit,nind,rbdnnz
      double precision delta,eps

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer rbdind(nind)
      character rbdtype(nind)
      double precision d(nind),g(nind),l(nind),lambda(m),rho(m),u(nind),
     +        x(*)

C     This subroutine implements the Conjugate Gradients method for 
C     minimizing the quadratic approximation q(d) of L(x,lambda,rho) 
C     at x
C
C     q(d) = 1/2 d^T H d + g^T d,
C
C     where H is an approximation of the Hessian matrix of the 
C     Augmented Lagrangian and g is its gradient vector,
C
C     subject to || d || <= delta and l <= x + d <= u.
C
C     In the constraint ''|| d || <= delta'', the norm will be the
C     Euclidian-norm if the input parameter trtype is equal to 0, and
C     it will be the sup-norm if trtype is equal to 1.
C
C     The method returns an approximation d of the solution such that
C 
C     (a) ||H d + g||_2 <= eps * ||g||_2, 
C
C     (b) ||d|| = delta or x + d is in the boundary of the box, or
C
C     (c) ( p such that p^t H p = 0 ) and ( d = - amax g if such p was 
C         found during the first CG iteration or the current point d 
C         of CG if such p was found in any other iteration ).
C
C     inform   integer
C              termination parameter:
C
C              0 = convergence with ||H d + g||_2 <= eps * ||g||_2;
C
C              1 = convergence to the boundary of ||d|| <= delta;
C
C              2 = convergence to the boundary of l <= x + d <= u;
C
C              3 = stopping with d = dk  such that <gk,dk> <= - theta 
C                  ||gk||_2 ||dk||_2 and <gk,d_{k+1}> > - theta 
C                  ||gk||_2 ||d_{k+1}||_2;
C
C              4 = not enough progress of the quadratic model during
C                  maxitnqmp iterations, i.e., during maxitnqmp 
C                  iterations | q - qprev | <= macheps * max( | q |, 1 )
C
C              6 = very similar consecutive iterates, for two 
C                  consecutive iterates x1 and x2 we have that
C
C                  | x2(i) - x1(i) | <= macheps * max ( | x1(i) |, 1 )
C
C                  for all i.
C
C              7 = stopping with p such that p^T H p = 0 and g^T p = 0;
C
C              8 = too many iterations;

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      character * 6 prectmp
      logical goth,gotp,negcur,restarted,samep
      integer i,itertmp,itnqmp
      double precision alpha,abox,amax,amaxbds,amaxbdsi,amaxdel,
     +        bestprog,bbeta,currprog,dnorm2,gnorm2,gtd,gtp,pnorm2,ptd,
     +        pthp,ptr,q,qprev,rnorm2,ztrprev,ztr,znorm2

C     LOCAL ARRAYS
      double precision p(nmax),hp(nmax),r(nmax),z(nmax)

C     ==================================================================
C     Initialization
C     ==================================================================

      restarted = .false.

 001  continue

      goth = .false.
      gotp = .false.

C     gnorm2 = norm2s(nind,g)
      gnorm2 = 0.0d0
      do i = 1,nind
          gnorm2 = gnorm2 + g(i) ** 2
      end do

      iter     =      0
      itnqmp   =      0
      qprev    = bignum
      bestprog =  0.0d0

      do i = 1,nind
          d(i) = 0.0d0
          r(i) =  g(i)
      end do

      q        =  0.0d0
      gtd      =  0.0d0
      dnorm2   =  0.0d0
      rnorm2   = gnorm2

      ztr      =  0.0d0

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprintinn .ge. 5 ) then
          write(*, 980) maxit,eps,delta,precond,hptype
          write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q

          write(10,980) maxit,eps,delta,precond,hptype
          write(10,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
      end if

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     if ||r||_2 = ||H d + g||_2 <= eps * ||g||_2 then stop

      if ( iter .ne. 0 .and. 
     +   ( ( rnorm2 .le. eps ** 2 * gnorm2 .and. iter .ge. 4 ) .or. 
     +     ( rnorm2 .le. macheps ) ) ) then 

          cginfo = 0

          if ( iprintinn .ge. 5 ) then
              write(*, 990)
              write(10,990)
          end if
  
          go to 500

      end if

C     if the maximum number of iterations was achieved then stop

      if ( iter .ge. max(4, maxit) ) then

          cginfo = 8

          if ( iprintinn .ge. 5 ) then
              write(*, 998)
              write(10,998)
          end if
  
          go to 500

      end if

C     ==================================================================
C     Preconditioner
C     ==================================================================

      if ( precond .eq. 'NONE' ) then

          do i = 1,nind
              z(i) = r(i)
          end do

          ztrprev = ztr
          ztr     = rnorm2
          znorm2  = rnorm2

      else if ( precond .eq. 'QNCGNA' ) then

          call capplyhpre(nind,m,rho,equatn,gotp,r,z)

          ztrprev = ztr

          ztr = 0.0d0
          do i = 1,nind
              ztr = ztr + z(i) * r(i)
          end do

C         znorm2 = norm2s(nind,z)
          znorm2 = 0.0d0
          do i = 1,nind
              znorm2 = znorm2 + z(i) ** 2
          end do

      end if

C     ==================================================================
C     Compute direction
C     ==================================================================

      if ( iter .eq. 0 ) then

          do i = 1,nind
              p(i) = - z(i)
          end do

          ptr    = - ztr
          pnorm2 =   znorm2

      else

          bbeta = ztr / ztrprev

          do i = 1,nind
              p(i) = - z(i) + bbeta * p(i)
          end do

          if ( precond .eq. 'NONE' ) then

              pnorm2 = rnorm2 - 2.0d0 * bbeta * ( ptr + alpha * pthp ) 
     +               + bbeta ** 2 * pnorm2
              ptr = - rnorm2 + bbeta * ( ptr + alpha * pthp )

          else if ( precond .eq. 'QNCGNA' ) then

              ptr = 0.0d0
              pnorm2 = 0.0d0
              do i = 1,nind
                  ptr = ptr + p(i) * r(i)
                  pnorm2 = pnorm2 + p(i) ** 2
              end do

          end if

      end if

C     Force p to be a descent direction of q(d), i.e.,
C     <\nabla q(d), p> = <H d + g, p> = <r, p> \le 0.

      if ( ptr .gt. 0.0d0 ) then

          do i = 1,nind
              p(i) = - p(i)
          end do

          ptr = - ptr

      end if

C     ==================================================================
C     Compute p^T H p
C     ==================================================================

C     hp = H p

      call calchalp(nind,x,m,lambda,rho,equatn,linear,p,hp,goth,inform)
      if ( inform .lt. 0 ) return

C     Compute p^T hp

      pthp = 0.0d0
      do i = 1,nind
          pthp = pthp + p(i) * hp(i)
      end do 

C     ==================================================================
C     Compute maximum steps
C     ==================================================================

      amaxdel = bignum

      do i = 1,nind
          if ( p(i) .gt. 0.0d0 ) then
              amaxdel = min( amaxdel,  (   delta - d(i) ) / p(i) )
          else if ( p(i) .lt. 0.0d0 ) then
              amaxdel = min( amaxdel,  ( - delta - d(i) ) / p(i) )
          end if
      end do

      amaxbds = bignum

C     do i = 1,nind
C        if ( p(i) .gt. 0.0d0 ) then

C            amaxbdsi = ( u(i) - d(i) - x(i) ) / p(i)

C            if ( amaxbdsi .lt. amaxbds ) then
C                amaxbds    = amaxbdsi
C                rbdnnz     = 1
C                rbdind(1)  = i
C                rbdtype(1) = 'U'
C            else if ( amaxbdsi .eq. amaxbds ) then
C                rbdnnz          = rbdnnz + 1
C                rbdind(rbdnnz)  = i
C                rbdtype(rbdnnz) = 'U'
C            end if

C        else if ( p(i) .lt. 0.0d0 ) then

C            amaxbdsi = ( l(i) - d(i) - x(i) ) / p(i)

C            if ( amaxbdsi .lt. amaxbds ) then
C                amaxbds    = amaxbdsi
C                rbdnnz     = 1
C                rbdind(1)  = i
C                rbdtype(1) = 'L'
C            else if ( amaxbdsi .eq. amaxbds ) then
C                rbdnnz          = rbdnnz + 1
C                rbdind(rbdnnz)  = i
C                rbdtype(rbdnnz) = 'L'
C            end if

C        end if
C     end do

      amax = min( amaxdel, amaxbds )

C     ==================================================================
C     Compute the step
C     ==================================================================

      negcur = .false.

C     If p^T H p > 0 then take the conjugate gradients step

      if ( pthp .gt. 0.0d0 ) then

          alpha = min( amax, ztr / pthp )

C     Else, if we are at iteration zero then take the maximum 
C     positive step in the minus gradient direction

      else if ( iter .eq. 0 ) then

          alpha = amax

          negcur = .true.

C     Otherwise, stop at the current iterate

      else

          cginfo = 7

          if ( iprintinn .ge. 5 ) then
              write(*, 997)
              write(10,997)
          end if
  
          go to 500

      end if

C     ==================================================================
C     Test the angle condition
C     ==================================================================

      ptd = 0.0d0
      gtp = 0.0d0
      do i = 1,nind
          ptd = ptd + p(i) * d(i)
          gtp = gtp + g(i) * p(i)
      end do

C     These are gtd and dnorm2 for the new direction d which was not 
C     computed yet.
      gtd = gtd + alpha * gtp
      dnorm2 = dnorm2 + alpha ** 2 * pnorm2 + 2.0d0 * alpha * ptd

      if ( gtd .gt. 0.0d0 .or. 
     +     gtd ** 2 .lt. theta ** 2 * gnorm2 * dnorm2 ) then

          if ( precond .ne. 'NONE' .and. iter .eq. 0 ) then

              if ( iprintinn .ge. 5 ) then
                  write(*, 986)
                  write(10,986)
              end if

              restarted = .true.
              itertmp   = iter
              prectmp   = precond
              precond   = 'NONE'
              go to 001

          end if

          cginfo = 3

          if ( iprintinn .ge. 5 ) then
              write(*, 993)
              write(10,993)
          end if

          go to 500

      end if

C     ==================================================================
C     Compute the quadratic model functional value at the new point
C     ==================================================================

      qprev = q

      q = q + 0.5d0 * alpha ** 2 * pthp + alpha * ptr

C     ==================================================================
C     Compute new d
C     ==================================================================

      do i = 1,nind
          d(i) = d(i) + alpha * p(i)
      end do

C     ==================================================================
C     Compute the residual r = H d + g
C     ==================================================================

      do i = 1,nind
          r(i) = r(i) + alpha * hp(i)
      end do

C     rnorm2 = norm2s(nind,r)
      rnorm2 = 0.0d0
      do i = 1,nind
          rnorm2 = rnorm2 + r(i) ** 2
      end do

C     ==================================================================
C     Increment number of iterations
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Print information of this iteration
C     ==================================================================

      if ( iprintinn .ge. 5 ) then
          write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q
          write(10,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
      end if

C     ==================================================================
C     Test other stopping criteria
C     ==================================================================

C     Boundary of the "trust region"

      if ( alpha .eq. amaxdel ) then

          cginfo = 1

          if ( iprintinn .ge. 5 ) then
              if ( negcur ) then
                  write(*, 987)
                  write(10,987)
              end if

              write(*, 991)
              write(10,991)
          end if
  
          go to 500

      end if

C     Boundary of the box constraints

C     if ( alpha .eq. amaxbds ) then

C         cginfo = 2

C         if ( iprintinn .ge. 5 ) then
C             if ( negcur ) then
C                 write(*, 987)
C                 write(10,987)
C             end if

C             write(*, 992)
C             write(10,992)
C         end if
  
C         go to 500

C     end if

C     Small useful proportion

      abox = bignum

      do i = 1,nind
         if ( d(i) .gt. 0.0d0 ) then
             abox = min( abox, ( u(i) - x(i) ) / d(i) )
         else if ( d(i) .lt. 0.0d0 ) then
             abox = min( abox, ( l(i) - x(i) ) / d(i) )
         end if
      end do

      if ( abox .le. 0.1d0 ) then

          cginfo = 5

          if ( iprintinn .ge. 5 ) then
              write(* ,995)
              write(10,995)
          end if
  
          go to 500

      end if

C     Two consecutive iterates are too much close

      samep = .true.
      do i = 1,nind
         if ( abs( alpha * p(i) ) .gt. 
     +        macheps * max( 1.0d0, abs( d(i) ) ) ) then 
             samep = .false.
          end if
      end do

      if ( samep ) then

          cginfo = 6

          if ( iprintinn .ge. 5 ) then
              write(*, 996)
              write(10,996)
          end if
  
          go to 500

      end if

C     Many iterations without good progress of the quadratic model 

      currprog = qprev - q
      bestprog = max( currprog, bestprog )

      if ( currprog .le. epsnqmp * bestprog ) then

          itnqmp = itnqmp + 1

          if ( itnqmp .ge. maxcgitnp ) then
              cginfo = 4

              if ( iprintinn .ge. 5 ) then
                  write(*, 994)
                  write(10,994)
              end if

              go to 500
          endif

      else
          itnqmp = 0
      endif

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     ==================================================================
C     Return
C     ==================================================================

 500  continue

C     Print final information

      if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
          write(*, 985) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
          write(10,985) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
      end if

      if ( restarted ) then
          iter = iter + itertmp
          precond = prectmp
      end if

      return

C     Non-executable statements

 980  format(/,5X,'Conjugate Gradients (maxit = ',I7,',',1X,'eps = ',
     +            1P,D7.1,',',1X,'delta = ',1P,D7.1,')',
     +       /,5X,'(Preconditioner: ',A6,',',1X,
     +            'Hessian-vector product type: ',A6,')')
 984  format(  5X,'CG iter = ',I7,' rnorm = ',1P,D10.4,' dnorm = ',
     +            1P,D10.4,' q = ',1P,D11.4)
 985  format(/,5X,'Truncated Newton direction (first ',I7,
     +            ' components): ',/,1(5X,6(1X,1P,D11.4)))
 986  format(  5X,'The first CG-PREC iterate did not satisfy the ',
     +            'angle condition. CG will be restarted without ',
     +            'preconditioner)')
 987  format(  5X,'p such that p^T H p = 0 was found. ',
     +            'Maximum step was taken.')

 990  format(  5X,'Flag of CG: Convergence with small residual.')
 991  format(  5X,'Flag of CG: Convergence to the trust region ',
     +            'boundary.')
C992  format(  5X,'Flag of CG: Convergence to the box boundary.')
 993  format(  5X,'Flag of CG: The next CG iterate will not satisfy ',
     +            'the angle condition.')
 994  format(  5X,'Flag of CG: Not enough progress in the quadratic ',
     +            'model.')
 995  format(  5X,'Flag of CG: The maximum step to remain within the ',
     +            'box is smaller than 0.1)')
 996  format(  5X,'Flag of CG: Very near consecutive iterates.')
 997  format(  5X,'Flag of CG: p such that p^T H p = 0 was found.')
 998  format(  5X,'Flag of CG: Maximum number of CG iterations ',
     +            'reached.')

      end
C     *****************************************************************
C     *****************************************************************

      subroutine betra(n,nind,x,l,u,m,lambda,rho,equatn,linear,f,g,
     +trdelta,newdelta,mslamb,epsg,xeucn,xsupn,gpsupn,d,rbdind,rbdtype,
     +xtrial,ftrial,gtrial,triter,ispgiter,chcnt,memfail,betinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical memfail
      integer betinfo,chcnt,inform,ispgiter,m,n,nind,triter
      double precision epsg,f,ftrial,gpsupn,mslamb,newdelta,trdelta,
     +        xeucn,xsupn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      character rbdtype(nind)
      integer rbdind(nind)
      double precision d(nind),g(n),gtrial(n),l(n),lambda(m),rho(m),
     +        u(n),x(n),xtrial(n)

C     Solves the "unconstrained" inner problem that arises from the 
C     active-set method to solve box-constrained minimization problem
C
C            Minimize f(x) subject to l <= x <= u
C     
C     described in 
C
C     M. Andretta, E. G. Birgin e J. M. Martínez. "Practical active-set 
C     Euclidian trust-region method with spectral projected gradients 
C     for bound-constrained minimization". Optimization 54, pp. 
C     305-325, 2005.

C     betinfo:
C
C     0: Sufficient decrease of the function
C     1: x hits the boundary
C     2: Unbounded objective function?
C     3: Either the search direction or the step lenght is too small
C     4: Trust-region radius too small
C     5: Undefined direction
C     6: x is a first-order stationary point close to the boundary
C     7: x is a second-order stationary point

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/
C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/

C     LOCAL SCALARS
      logical frstit,pd,samep
      integer dinfo,extinfo,i,col,lin,index,lsinfo,rbdnnz
      double precision amax,ared,dbound,deucn,gsupn,lamspg,phi,pred,
     +        stpl,tmp

C     LOCAL ARRAYS
      integer hdiag(nmax)

C     EXTERNAL SUBROUTINES
      external calcal

C     Print presentation information

      if ( iprintinn .ge. 5 ) then
          write(*, 1000)
          write(10,1000)
      end if
      
C     Initialization

      if ( .not. sameface ) then
          trdelta = max( trdelta, newdelta )
      end if
      
      newdelta = 0.0d0
      
      memfail  = .false.

      frstit   = .true.

      ittype   = 4

C     ==================================================================
C     Compute distance to the boundary
C     ==================================================================
      
C     step 1: calculate the distance between x and the boundary.
C             dbound is the largest positive delta such that the ball 
C             centered at x with radius delta is still inside the box 
C             (set of constraints).

      dbound = bignum

      do i = 1,nind
          dbound = min( dbound, x(i) - l(i) )
      end do
      
      do i = 1,nind
          dbound = min( dbound, u(i) - x(i) )
      end do

C     Calculate infinite-norm of gradient.

      gsupn = 0.0d0
      do i = 1,nind
          gsupn = max( gsupn, abs( g(i) ) )
      end do

C     ==================================================================
C     Close to the boundary: perform inner SPG iteration
C     ==================================================================
      
C     step 2: close to the boundary: if the gradient is null then 
C             algorithm stops. Else, perform an inner SPG iteration.

      if ( dbound .lt. 2.0d0 * trdelmin ) then

          if ( gsupn .le. epsg ) then

              do i = 1,n
                  xtrial(i) = x(i)
                  gtrial(i) = g(i)
              end do
              
              ftrial  = f

              betinfo = 6
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9010) 
                  write(10,9010) 
              end if
              
              return
          end if

C         Inner SPG iteration

          if ( iprintinn .ge. 5 ) then
              write(*, 2070)
              write(10,2070)
          end if

          ittype = 5

C         Compute spectral steplength
          
          if ( sty .le. 0.0d0 ) then
              lamspg = max( 1.0d0, xsupn / gpsupn )
          else
              lamspg = sts / sty
          end if
          lamspg = min( lspgma, max( lspgmi, lamspg ) )
          
          ispgiter = ispgiter + 1

          call spgls(nind,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg,
     +    xtrial,ftrial,stpl,d,calcal,lsinfo,inform)
          
          if ( inform .lt. 0 ) return

          if ( lsinfo .eq. 0 .or. lsinfo .eq. 1 ) then
              betinfo = 0
          elseif ( lsinfo .eq. 2 ) then
              betinfo = 2
          elseif ( lsinfo .eq. 3 ) then
              betinfo = 3
          end if

          if ( betinfo .eq. 0 ) then

              if ( stpl .ge. 1.0d0 .and. extrp2 ) then

C                 Calculate amax, which is the biggest positive scalar 
C                 such that x + amax*d is still in the box bounded by l 
C                 and u.

                  call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)

                  call extrapolation(nind,x,l,u,m,lambda,rho,equatn,
     +            linear,g,xtrial,ftrial,gtrial,d,stpl,amax,rbdnnz,
     +            rbdind,rbdtype,fmin,beta,etaext,maxextrap,extinfo,
     +            inform)

                  if ( inform .lt. 0 ) return

                  if ( extinfo .eq. 2 ) then
                      betinfo = 2
                  else
                      betinfo = 0
                  end if

                  if ( betinfo .eq. 2 ) then
                      if ( iprintinn .ge. 5 ) then
                          write(*, 9050) 
                          write(10,9050) 
                      end if
                      
                      return
                  else
                      if ( stpl .ge. amax ) then
                          
                          betinfo = 1
                          
                          if ( iprintinn .ge. 5 ) then
                              write(*, 9030) 
                              write(10,9030) 
                          end if
                          
                      else
                          
                          betinfo = 0
                          
                          if ( iprintinn .ge. 5 ) then
                              write(*, 9000) 
                              write(10,9000) 
                          end if
                      end if

                      return
                  end if
              end if
          end if
          
C         Compute the gradient at the new iterate
          
          call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,gtrial,
     +    inform)
          if ( inform .lt. 0 ) return

          if ( betinfo .eq. 2 ) then

              if ( iprintinn .ge. 5 ) then
                  write(*, 9050) 
                  write(10,9050) 
              end if
              
              return
              
          elseif ( betinfo .eq. 3 ) then 

              if ( iprintinn .ge. 5 ) then
                  write(*, 9060) 
                  write(10,9060) 
              end if
              
              return
              
          else
              
C             If new point x is in the boundary, the inner algorithm 
C             stops and returns x as solution.
              
              rbdnnz = 0
      
              do i = 1,nind
                  if ( xtrial(i) .le. l(i) + 
     +                 macheps23 * max( abs( l(i) ), 1.0d0 ) ) then
                      xtrial(i) = l(i)
                      rbdnnz = rbdnnz + 1
                  else if ( xtrial(i) .ge. u(i) - macheps23 * 
     +                     max( abs( u(i) ), 1.0d0 ) ) then
                      xtrial(i) = u(i)
                      rbdnnz = rbdnnz + 1
                  end if
              end do
              
              if ( rbdnnz .gt. 0 ) then
                  
                  betinfo = 1
                  
                  if ( iprintinn .ge. 5 ) then
                      write(*, 9030) 
                      write(10,9030) 
                  end if

              else
                  
                  betinfo = 0
                  
                  if ( iprintinn .ge. 5 ) then
                      write(*, 9000) 
                      write(10,9000) 
                  end if

              end if

              return

          end if
          
      end if

C     ==================================================================
C     Far from the boundary: perform trust-region iteration
C     ==================================================================
      
      triter = triter + 1

C     step 3: far from the boundary, solve trust-region subproblem. 
C             Evaluate function Hessian at x.

      call calchal(nind,x,m,lambda,rho,equatn,linear,hlin,hcol,hval,
     +hnnz,inform)
      if ( inform .lt. 0 ) return
      
      do i = 1,nind
          hdiag(i) = 0
      end do

      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)

          if ( lin .eq. col ) then
              if ( hdiag(lin) .eq. 0 ) then
                  hdiag(lin) = i
              else
                  hval(hdiag(lin)) = hval(hdiag(lin)) + hval(i)
                  hval(i) = 0.0d0
              end if 
          end if
      end do

      do i = 1,nind
         if ( hdiag(i) .eq. 0 ) then
            hnnz       = hnnz + 1            
            hlin(hnnz) = i
            hcol(hnnz) = i
            hval(hnnz) = 0.0d0
            hdiag(i)   = hnnz
         end if
      end do

C     step 4: solve the trust-region subproblem using More-Sorensen's 
C             algorithm to minimize "exactly" quadratics subjected to 
C             balls. 

C     If trust-region radius is too small, the inner algorithm stops.

 100  continue

      if ( trdelta .lt. macheps * max( 1.0d0, xeucn ) ) then

          do i = 1,n
              xtrial(i) = x(i)
              gtrial(i) = g(i)
          end do
          
          ftrial  = f

          betinfo = 4
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9040) 
              write(10,9040) 
          end if
          
          return
      end if
      
      call moresor(nind,g,hnnz,hlin,hcol,hval,hdiag,trdelta,mssig,
     +0.0d0,mseps,msmaxit,mslamb,pd,d,chcnt,memfail,dinfo)

      if ( memfail ) then
          return
      end if

C     If maximum allowed number of MEQB iterations is achieved, another
C     direction d is calculated.

      if ( dinfo .eq. 5 ) then 
          
          if ( iprintinn .ge. 5 ) then
              write(*, 2000)
              write(10,2000)
          end if

          call dogleg(nind,g,hnnz,hlin,hcol,hval,pd,trdelta,d,dinfo)

      end if
      
C     If both internal gradient and Hessian matrix are null, subroutines 
C     MEQB and dogleg stop with dinfo = 0 and then the inner algorithm 
C     stops declaring "second-order stationary point".

      if ( dinfo .eq. 0 ) then

          do i = 1,n
              xtrial(i) = x(i)
              gtrial(i) = g(i)
          end do
          
          ftrial  = f

          betinfo = 7
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9020) 
              write(10,9020) 
          end if
          
          return
      end if
      
C     Print direction

      if ( iprintinn .ge. 5 ) then
          write(*, 1020) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
          write(10,1020) min0(nind,nprint),(d(i),i=1,min0(nind,nprint))
      end if

C     Evaluate the quadratic model of the objective function at d.

      call squad(nind,d,g,hnnz,hlin,hcol,hval,phi)

C     If the value of the quadratic model at d is 0 it means that x is a 
C     second-order stationary point. In this case, inner algorithm stops
C     declaring this.
      
      if ( ( abs( phi ) .le. phieps ) .and. ( gsupn .le. epsg ) ) then

          do i = 1,n
              xtrial(i) = x(i)
              gtrial(i) = g(i)
          end do
          
          ftrial  = f

          betinfo = 7
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9020) 
              write(10,9020) 
          end if
          
          return
      end if
      
C     Calculate predicted decrease of objective function
      
      pred = abs( phi )
      
C     Calculate d Euclidian-norm
      
      deucn = 0.0d0
      do i = 1,nind
          deucn = deucn + d(i)**2
      end do
      deucn = sqrt( deucn )
      
C     To avoid NaN and Inf directions

      if ( .not. ( deucn .le. bignum ) ) then
      
          trdelta = 2.5d-1 * trdelta

          if ( iprintinn .ge. 5 ) then
              write(*, 2060) 
              write(10,2060) 
          end if

          if ( trdelta .lt. macheps * max( 1.0d0, xeucn ) ) then

              do i = 1,n
                  xtrial(i) = x(i)
                  gtrial(i) = g(i)
              end do
              
              ftrial  = f

              betinfo = 5
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9070) 
                  write(10,9070) 
              end if
              
              return
          end if

          triter  = triter + 1
          frstit  = .false.

          go to 100

      end if

C     Calculate point xtrial = x + d.
      
      do i = 1,nind
          xtrial(i) = x(i) + d(i)
      end do

      stpl = 1.0d0

C     Verify if xtrial is inside de box. If not, interior is set to
C     false and amax is set to the biggest positive scalar such that
C     x + amax*d is inside the box.

      call compamax(nind,x,l,u,d,amax,rbdnnz,rbdind,rbdtype)
      
C     ==================================================================
C     Point on the boundary
C     ==================================================================      
      
C     If xtrial is not interior to the box, xtrial = x + d is replaced
C     by xtrial = x + amax*d. Now xtrial is definitely interior. Actually,
C     it is in the boundary. If the objective function decreases in
C     xtrial, the inner algorithm stops with xtrial as a solution.
C     Otherwise, a new trust-region radius trdelta is chosen (smaller than
C     dbound) and a new quadratic model minimizer is calculated (which
C     is necessarily interior because of the choice of trdelta).

      if ( amax .le. 1.0d0 ) then 
          
          if ( iprintinn .ge. 5 ) then
              write(*, 2010) 
              write(10,2010) 
          end if
          
          do i = 1,nind
              xtrial(i) = x(i) + amax * d(i)
          end do
          
          stpl = amax

C         Set x(i) to l(i) or u(i) for the indices i that got to the 
C         boundary (to prevent errors).
          
          do i = 1,rbdnnz
              index = rbdind(i)
              if ( rbdtype(i) .eq. 'L' ) then
                  xtrial(index) = l(index)
              elseif ( rbdtype(i) .eq. 'U' ) then
                  xtrial(index) = u(index)
              end if
          end do
          
          call calcal(nind,xtrial,m,lambda,rho,equatn,linear,ftrial,
     +    inform)
          if ( inform .lt. 0 ) return
          
C         Print functional value

          if ( iprintinn .ge. 5 ) then
              write(*, 1030) trdelta,ftrial,fcnt
              write(10,1030) trdelta,ftrial,fcnt
          end if

C         Test whether f is very small
      
          if ( ftrial .le. fmin ) then
              
              call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,
     +        gtrial,inform)
              if ( inform .lt. 0 ) return

              betinfo = 2
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9050) 
                  write(10,9050) 
              end if
              
              return
          end if      
          
C         If the new point x + d is too close to the previous point x, 
C         inner algorithm stops

          samep = .true.
          do i = 1,nind
              if ( abs( x(i) - xtrial(i) ) .gt. 
     +             macheps * max( abs( xtrial(i) ), 1.0d0 ) ) then
                  samep = .false.
              end if
          end do
          
          if ( samep .and. ftrial - f .le. macheps23 * abs(f) ) then
              
              betinfo = 3
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 9060) 
                  write(10,9060) 
              end if
              
              return
          end if

C         Test if function value decreases at xtrial
        
          if ( ( ftrial .le. f ) .or. 
     +         ( ( deucn .le. macheps23 * xeucn ) .and.
     +         ( ftrial .le. f + macheps23 * abs( f ) ) ) ) then
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 2030) 
                  write(10,2030) 
              end if
              
              if ( extrp4 ) then
                  
                  call extrapolation(nind,x,l,u,m,lambda,rho,equatn,
     +            linear,g,xtrial,ftrial,gtrial,d,stpl,amax,rbdnnz,
     +            rbdind,rbdtype,fmin,beta,etaext,maxextrap,extinfo,
     +            inform)

                  if ( inform .lt. 0 ) return

                  if ( extinfo .eq. 2 ) then
                      betinfo = 2
                  else
                      betinfo = 0
                  end if

C                 Update the trust-region radius (which may or may not
C                 be used)

                  if ( frstit ) then
                      trdelta = 2.0d0 * max( trdelta, stpl * deucn )
                  end if
                  
                  if ( betinfo .eq. 2 ) then
                      if ( iprintinn .ge. 5 ) then
                          write(*, 9050) 
                          write(10,9050) 
                      end if
                      
                      return
                  end if
                  
                  betinfo = 1
                  
C                 Calculate actual reduction of objective function

                  ared = f - ftrial
                  
                  go to 200
                  
              else
                  
C                 Update the trust-region radius (which may or may not
C                 be used)
               
                  trdelta = 2.0d0 * trdelta
                  
C                 Calculate actual reduction of objective function

                  ared = f - ftrial
                  
C                 Compute gtrial.
                  
                  call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,
     +            gtrial,inform)
                  if ( inform .lt. 0 ) return
                  
                  betinfo = 1
                  
                  go to 200
              end if
              
          else
              tmp      = trdelmin + msrho * 
     +                   ( ( dbound / (1.0d0 + mssig ) ) - trdelmin )
              newdelta = trdelta
              trdelta  = max( trdelmin, tmp )

              triter   = triter + 1
              frstit   = .false.

              if ( iprintinn .ge. 5 ) then
                  write(*, 2040) 
                  write(10,2040) 
              end if

              go to 100
          end if
      end if

C     ==================================================================
C     Point interior to the box
C     ==================================================================
      
C     step 5: in this case xtrial is inside the box. Acceptance or
C             rejection of the trust-region subproblem solution.

      call calcal(nind,xtrial,m,lambda,rho,equatn,linear,ftrial,inform)
      if ( inform .lt. 0 ) return

C     Print functional value

      if ( iprintinn .ge. 5 ) then
          write(*, 1030) trdelta,ftrial,fcnt
          write(10,1030) trdelta,ftrial,fcnt
      end if

C     Test whether f is very small
      
      if ( ftrial .le. fmin ) then
          
          call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,gtrial,
     +    inform)
          if ( inform .lt. 0 ) return

          betinfo = 2
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9050) 
              write(10,9050) 
          end if
          
          return
      end if
      
C     If the new point x + d is too close to the previous point x, inner
C     algorithm stops

      samep = .true.
      do i = 1,nind
          if ( abs( x(i) - xtrial(i) ) .gt. 
     +         macheps * max( abs( xtrial(i) ), 1.0d0 ) ) then
              samep = .false.
          end if
      end do
      
      if ( samep .and. ftrial - f .le. macheps23 * abs(f) ) then
          
          betinfo = 3
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9060) 
              write(10,9060) 
          end if
          
          return
      end if
      
C     Calculate actual reduction of objective function
      
      ared = f - ftrial
      
C     If there is not sufficient decrease of the function, the 
C     trust-region radius is decreased and the new quadratic model 
C     minimizer will be calculated. 

      if ( iprintinn .ge. 5 ) then
          write(*, 1010) deucn,pred,ared
          write(10,1010) deucn,pred,ared
      end if

      if ( deucn .le. macheps23 * xeucn .and.
     +     ared  .le. macheps23 * abs( f ) ) then
         
          call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,gtrial,
     +    inform)
          if ( inform .lt. 0 ) return
          
          go to 200

      end if
      
      if ( ( pred .le. phieps .or. ared .ge. tralpha*pred ) .and.
     +     ( pred .gt. phieps .or. f .ge. ftrial ) ) then

C         If extrapolation at step 5 is not to be performed, point 
C         xtrial is accepted.

          if ( extrp5 ) then
          
              call extrapolation(nind,x,l,u,m,lambda,rho,equatn,linear,
     +        g,xtrial,ftrial,gtrial,d,stpl,amax,rbdnnz,rbdind,rbdtype,
     +        fmin,beta,etaext,maxextrap,extinfo,inform)
              
              if ( inform .lt. 0 ) return
              
              if ( extinfo .eq. 2 ) then
                  betinfo = 2
              else
                  betinfo = 0
              end if
              
              if ( frstit ) then
                  trdelta = max( trdelta, stpl * deucn )
              end if
              
              if ( betinfo .eq. 2 ) then
                  if ( iprintinn .ge. 5 ) then
                      write(*, 9050) 
                      write(10,9050) 
                  end if
                  
                  return
              end if
              
C             Update actual reduction of objective function
              
              ared = f - ftrial
              
          else
              
              call calcnal(nind,xtrial,m,lambda,rho,equatn,linear,
     +        gtrial,inform)
              if ( inform .lt. 0 ) return
              
          end if

      else

          trdelta = 2.5d-1 * deucn
          triter  = triter + 1
          frstit  = .false.

          if ( iprintinn .ge. 5 ) then
              write(*, 2050) 
              write(10,2050) 
          end if
          
          go to 100

      end if
      
C     ==================================================================
C     Prepare for next call to this routine
C     ==================================================================
      
C     Update the trust-region radius (which may or may not be used). 
C     This update can only be done when the current iteration was a 
C     trust-region iteration (and not an inner SPG one).

 200  continue

      if ( ared .lt. 2.5d-1 * pred ) then
          trdelta = max( 2.5d-1 * deucn, trdelmin )
      else
          if ( ( ared .ge. 0.5d0 * pred ) .and.
     +         ( deucn .ge. trdelta - macheps23 * 
     +         max( trdelta, 1.0d0 ) ) ) then
              trdelta = max( 2.0d0 * trdelta, trdelmin )
          else
              trdelta = max( trdelta, trdelmin )
          end if
      end if
      
C     If new point x is in the boundary, the inner algorithm stops
C     and returns x as solution.

      if ( stpl .ge. amax ) then
          
          betinfo = 1
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9030) 
              write(10,9030) 
          end if
          
      else
          
          betinfo = 0
          
          if ( iprintinn .ge. 5 ) then
              write(*, 9000) 
              write(10,9000) 
          end if
      end if
      
C     Non-executable statements
      
 1000 format(/,5X,'Trust-region iteration')
 1010 format(  5X,'deucn = ',1P,D7.1,' pred = ',1P,D11.4,
     +            ' ared = ',1P,D11.4)
 1020 format(/,5X,'Trust-region direction (first ',I7,
     +            ' components): ',/,1(5X,6(1X,1P,D11.4)))
 1030 format(/,5X,'delta = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)

 2000 format(/,5X,'Since the direction More-Sorensen calculation ',
     +       /,5X,'failed, dogleg direction will be computed.')
 2010 format(/,5X,'x+d is not interior to the box.')
 2030 format(  5X,'f(x+d) < f(x), x+d is accepted.')
 2040 format(  5X,'f(x+d) >= f(x), a new direction d will be computed.')
 2050 format(  5X,'x+d did not obtain suficcient functional reduction. '
     +       /,5X,'A new direction d will be computed.')
 2060 format(/,5X,'Direction is undefined. ',
     +            'A new direction d will be computed.')
 2070 format(/,5X,'Point close to the boundary. ',
     +            'An inner SPG iteration will be used.')

 9000 format(  5X,'Flag of TR: Sufficient decrease of function.')
 9010 format(  5X,'Flag of TR: ',
     +            'First-order stationary point close to boundary.')
 9020 format(/,5X,'Flag of TR: Second-order stationary point.')
 9030 format(  5X,'Flag of TR: Point on the boundary.')      
 9040 format(  5X,'Flag of TR: Trust-region radius too small.')
 9050 format(  5X,'Flag of TR: Unbounded objective function?')
 9060 format(  5X,'Flag of TR: Very similar consecutive points.')
 9070 format(  5X,'Flag of TR: Undefined direction.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine squad(nred,x,g,hnnz,hlin,hcol,hval,phi)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,nred
      double precision phi
      
C     ARRAY ARGUMENTS
      integer hcol(hnnz),hlin(hnnz)
      double precision g(nred),hval(hnnz),x(nred)

C     Evaluates the quadratic model phi(x) = 1/2 x^T H x + g^T x.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )

C     LOCAL SCALARS
      integer col,i,lin
      
C     LOCAL ARRAYS
      double precision wd(nmax)

      do i = 1,nred
          wd(i) = 0.0d0
      end do
      
      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          
          wd(lin) = wd(lin) + hval(i) * x(col)
          if ( lin .ne. col ) then
              wd(col) = wd(col) + hval(i) * x(lin)
          end if
      end do
      
      phi = 0.0d0
      do i = 1,nred
          phi = phi + wd(i) * x(i)
      end do
      
      phi = phi * 0.5d0
      
      do i = 1,nred
          phi = phi + g(i) * x(i)
      end do
      
      end      
      
C     *****************************************************************
C     *****************************************************************

C     Algorithm that finds a unconstrained minimizer of objective 
C     function inside the box of constraints, hits the boundary 
C     (obtaining function decrease), or finds an interior point where 
C     the objective function has sufficient decrease (compared to its 
C     value at x). Extrapolation may be done.
C
C     When the current point x is "close to" the boundary, a Spectral 
C     Projected Gradient (SPG) iteration is used to calculate the new 
C     point. If this new point is at the boundary, the algorithm stops. 
C     Otherwise, a new iteration begins.
C
C     When x is "far from" the boundary, trust-region radius is 
C     determined and d is calculated using More-Sorensen algorithm to 
C     solve the trust-region subproblem (which is to find a minimizer a 
C     to a function quadratic model provided that the minimizer's 
C     Euclidian-norm is smaller than a given delta). The new point is
C     xtrial = x + d. 
C
C     If xtrial lies outside the box of constraints, it is truncated on 
C     the boundary. This new y on the boundary will be candidate to be a 
C     solution. If function value at new xtrial is smaller than function 
C     value at x, inner algorithm stops with xtrial. Otherwise, the 
C     trust-region radius is decreased so that the new solution d' to 
C     the trust-region subproblem makes x + d' be interior to the box.
C     More-Sorensen algorithm is used to calculate d' too.
C
C     If xtrial lies inside the box, sufficient decrease of objective 
C     function is tested. If it is true, xtrial is accepted as a solution
C     candidate. If xtrial in on the boundary, inner algorithm stops and 
C     if it is interior, a new iteration begins. If sufficient decrease is 
C     not obtained, trust-region radius is decreased and a new quadratic 
C     model minimizer is calculated (as in a classical trust-region 
C     algorithm for unconstrained minimization).
C
C     If the user wants, after calculating the candidate solution xtrial, 
C     extrapolation may be performed. For this, set extrpi to true, 
C     where i is the step of the inner algorithm that can call 
C     extrapolation procedure. If extrp4 is true, extrapolation will be
C     tried after xtrial hit the boundary. And if extrp5 is true,
C     extrapolation will be tried after xtrial is calculated by
C     trust-region algorithm, when xtrial is interior and provides
C     sufficient decrease of objective function.
C
C     If gradient at current point is null, inner algorithm stops 
C     declaring "first-order stationary point". If quadratic model 
C     minimum is 0, inner algorithm stops declaring "second-order 
C     stationary point".
C
C     M. Andretta, E. G. Birgin and J. M. Martinez, ''Practical active-set
C     Euclidian trust-region method with spectral projected gradients for
C     bound-constrained minimization'', Optimization 54, pp. 305-325, 2005.
C
C     On Entry
C      
C     n        integer
C              dimension of full space
C
C     nind     integer
C              dimension of reduced space
C
C     x        double precision x(n)
C              initial point, interior to the current face
C
C     l        double precision l(n)
C              lower bounds on x
C
C     u        double precision u(n)
C              upper bounds on x
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              BETRA and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if BETRA is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.). Finally,
C              linear is logical array that, for each constraint, 
C              indicates whether the constraint is a linear constraint
C              (.true.) or a nonlinear constraint (.false.)
C
C     f        double precision
C              objective function value at x
C
C     g        double precision g(n)
C              gradient at x
C
C     trdelta  double precision 
C              trust-region radius
C
C     newdelta double precision
C              trust-region radius is set to the maximum between newdelta 
C              and trdelta
C     
C     mslamb   double precision
C              value that More-Sorensen algorithm calculates to find 
C              the trust-region subproblem solution (MEQB)
C
C     epsg     double precision
C              allowed error for projected gradient norm
C
C     xeucn    double precision
C              x Euclidian norm
C
C     xsupn    double precision
C              x sup-norm
C
C     gpeucn   double precision
C              projected-gradient Euclidian norm
C
C     On Return
C      
C     trdelta  double precision
C              updated trut-region radius
C
C     newdelta double precision
C              when the trust-region radius trdelta is decreased so that 
C              the point x + d fit the current face, newdelta is set to 
C              the previous value of trdelta. Otherwise, it is set to 0
C
C     mslamb   double precision
C              updated value for next iteration (see entry parameter)
C
C     d        double precision d(nind)
C              direction computed such that xtrial = x + d
C
C     rbdind   integer rbdind(n)
C              indices of variables that reached their bounds
C     
C     rbdtype  character rbdtype(n)
C              if variable rbdind(i) reached its lower bound, 
C              rbdtype(i) = 'L'. If variable rbdind(i) reached its upper
C              bound, rbdtype(i) = 'U'. 
C
C     xtrial   double precision xtrial(n)
C              solution candidate, with inform as described bellow
C
C     ftrial   double precision
C              objective function value at xtrial
C
C     gtrial   double precision gtrial(n)
C              gradient at xtrial
C
C     triter   integer
C              number of trust-region iterations
C
C     ispgiter integer
C              number of inner SPG iterations
C
C     chcnt    integer
C              number of Cholesky decompositions
C
C     memfail  logical
C              true iff linear solver failed because of lack of memory
C
C     betinfo  integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C      
C              0 = Sufficient decrease of the function;
C
C              1 = x hits the boundary;
C
C              2 = Unbounded objective function?
C
C              3 = Either the search direction or the step lenght is too
C                  small;
C
C              4 = Trust-region radius too small;
C
C              5 = Undefined direction;
C
C              6 = x is a first-order stationary point close to the
C                  boundary;
C
C              7 = x is a second-order stationary point;
C
C     inform   0 = no error occurred;
C             <0 = error in some function, gradient or Hessian routine.
C     ******************************************************************
C     ******************************************************************

      subroutine dogleg(n,g,hnnz,hlin,hcol,hval,pd,delta,p,dlinfo)
      
      implicit none

C     SCALAR ARGUMENTS
      logical pd
      integer dlinfo,hnnz,n
      double precision delta
      
C     ARRAY ARGUMENTS
      integer hcol(hnnz),hlin(hnnz)
      double precision g(n),hval(hnnz),p(n)

C     Compute approximate minimizer of the quadratic model using Dogleg
C     method when the Hessian is positive definite and Cauchy point
C     otherwise.

C     dlinfo:
C
C     0: successfull exit. Both H and g are null;
C     1: successfull exit.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer col,i,lin
      double precision a,b,c,coef,d,delta2,geucn,geucn2,gthg,pbeucn2,
     +        pueucn,pueucn2,putpb,r

C     LOCAL ARRAYS
      double precision pu(nmax)
C     Presentation

      if ( iprintinn .ge. 5 ) then
          write(* ,1000)
          write(10,1000)
      end if

C     Initialization

      dlinfo = 1
      delta2 = delta**2

C     If H is not positive definite, compute Cauchy point

      if ( .not. pd ) then 
          go to 100
      end if
      
      pbeucn2 = 0.0d0
      do i = 1,n
          pbeucn2 = pbeucn2 + p(i)**2
      end do
      
C     If Newton step is inside the trust region, this step is taken
      
      if ( pbeucn2 .le. delta2 ) then
          go to 500
      end if
      
C     If Newton step is outside the trust region, compute the
C     unconstrained minimizer pu of the quadratic function
      
C     Compute g^T H g
      
      do i = 1,n
          pu(i) = 0.0d0
      end do
      
      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          
          pu(lin) = pu(lin) + hval(i) * g(col)
          if ( lin .ne. col ) then
              pu(col) = pu(col) + hval(i) * g(lin)
          end if
      end do
      
      gthg = 0.0d0
      do i = 1,n
          gthg = gthg + g(i) * pu(i)
      end do
      
C     Compute g^T g.
      
      geucn2 = 0.0d0
      do i = 1,n
          geucn2 = geucn2 + g(i)**2
      end do
      
C     Compute pu = coef * g.
      
      coef = - geucn2 / gthg 
      do i = 1,n
          pu(i) = coef * g(i)
      end do
      
C     If uncontrained minimizer is outside the trust region, it is
C     truncated at the border
      
      pueucn2 = 0.0d0
      do i = 1,n
          pueucn2 = pueucn2 + pu(i)**2
      end do
      pueucn = sqrt( pueucn2 )
      
      if ( pueucn2 .ge. delta2 ) then
          r = delta / pueucn
          do i = 1,n
              p(i) = r * pu(i)
          end do
          go to 500
      end if
      
C     Compute step length in directions pu and pb. Direction p is a
C     linear combination of pu and pb. To compute the step length in
C     each direction, we have to solve (in r):
C     \| pu + (r - 1) (pb - pu) \|^2 = delta^2.
      
      putpb = 0.0d0
      do i = 1,n
          putpb = putpb + pu(i) * p(i)
      end do
      
      a = pbeucn2 + pueucn2 - 2.0d0 * putpb
      b = - 2.0d0 * pbeucn2 - 4.0d0 * pueucn2 + 6.0d0 * putpb
      c = pbeucn2 + 4.0d0 * pueucn2 - 4.0d0 * putpb - (delta**2)
      d = b**2 - 4.0d0 * a * c
      
      if ( ( 2.0d0 * abs( a ) .lt. macheps23 ) .or. 
     +     ( d .lt. 0.0d0 ) ) then
          go to 200
      end if
      
      r = (- b + sqrt( d )) / (2.0d0 * a)
      
      if ( ( r .lt. 1.0d0 ) .or. ( r .gt. 2.0d0 ) ) then
          r = (- b - sqrt( d )) / (2.0d0 * a)
      end if
      r = r - 1.0d0
      
      do i = 1,n
          p(i) = r * p(i) + (1.0d0 - r) * pu(i)
      end do
      
      go to 500
      
C     Compute Cauchy point, when H is not positive definite

 100  continue
      
C     Compute g^T H g

      do i = 1,n
          pu(i) = 0.0d0
      end do
      
      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          
          pu(lin) = pu(lin) + hval(i) * g(col)
          if ( lin .ne. col ) then
              pu(col) = pu(col) + hval(i) * g(lin)
          end if
      end do
      
      gthg = 0.0d0
      do i = 1,n
          gthg = gthg + g(i) * pu(i)
      end do
      
C     Compute g^T g
      
      geucn2 = 0.0d0
      do i = 1,n
          geucn2 = geucn2 + g(i)**2
      end do
      geucn = sqrt( geucn2 )
      
C     Compute step length coef
      
 200  if ( ( abs( gthg ) .le. macheps23 ) .and. 
     +     ( geucn .le. macheps23 ) ) then

          dlinfo = 0
          
          if ( iprintinn .ge. 5 ) then
              write(* ,1020)
              write(10,1020)
          end if

          return
      end if
      
      if ( ( gthg .le. 0.0d0 ) .or. 
     +     ( geucn2*geucn .lt. delta*gthg ) ) then
          coef = - delta / geucn
      else
          coef = - geucn2 / gthg 
      end if

C     Compute p = coef * g
      
      do i = 1,n
          p(i) = coef * g(i)
      end do

C     Termination

 500  continue

      if ( iprintinn .ge. 5 ) then
          write(* ,1010)
          write(10,1010)
      end if

      if ( iprintinn .ge. 5 .and. nprint .ne. 0 ) then
          write(*, 1030) min0(n,nprint),(p(i),i=1,min0(n,nprint))
          write(10,1030) min0(n,nprint),(p(i),i=1,min0(n,nprint))
      end if

 1000 format(/,5X,'Computation of Dogleg direction.')
 1010 format(  5X,'Dogleg computed successfully.')
 1020 format(  5X,'Null direction was computed.')
 1030 format(/,5X,'Dogleg direction (first ',I7,' components): ',
     +       /,1(5X,6(1X,1P,D11.4)))

      end
      
C     ******************************************************************
C     ******************************************************************

C     Compute approximate minimizer of the quadratic model using Dogleg
C     method when the Hessian is positive definite and Cauchy point
C     otherwise.

C     On Entry
C
C     n        integer
C              dimension
C
C     g        double precision g(n)
C              vector used to define the quadratic function
C
C     hnnz     integer
C              number of nonzero elements of H
C
C     hlin     integer hlin(hnnz)
C              row indices of nonzero elements of H
C
C     hcol     integer hcol(hnnz)
C              column indices of nonzero elements of H
C
C     hval     double precision hval(hnnz)
C              nonzero elements of H, that is,
C              H(hlin(i),hcol(i)) = hval(i). Since H is symmetric, just
C              one element H(i,j) or H(j,i) must appear in its sparse
C              representation. If more than one pair corresponding to
C              the same position of H appears in the sparse
C              representation, the multiple entries will be summed.
C              If any hlin(i) or hcol(i) is out of range, the entry will
C              be ignored
C     
C     pd       logical
C              indicates if the last Cholesky decomposition of moresor 
C              was successfull. That is, if the last matrix used by 
C              moresor was positive definite
C     
C     delta    double precision
C              trust-region radius
C     
C     On Return
C     
C     p        double precision p(n)      
C              solution to problem
C              minimize     psi(w)
C              subjected to ||w|| <= delta
C      
C     dlinfo   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C
C              0 = successfull exit. Both H and g are null;
C
C              1 = successfull exit.
C     ******************************************************************
C     ******************************************************************

      subroutine moresor(n,g,bnnz,blin,bcol,bval,bdiag,delta,sigma1,
     +sigma2,eps,maxit,l,pd,p,chcnt,memfail,msinfo)
      
      implicit none

C     SCALAR ARGUMENTS
      logical memfail,pd
      integer bnnz,chcnt,maxit,msinfo,n
      double precision delta,eps,l,sigma1,sigma2

C     ARRAY ARGUMENTS
      integer bcol(bnnz),bdiag(n),blin(bnnz)
      double precision bval(bnnz),g(n),p(n)

C     Solves the problem
C
C     minimize    psi(w) = 1/2 w^TBw + g^Tw 
C     subject to  ||w|| <= delta
C
C     Using the method described in "Computing a trust region step", 
C     by More and Sorensen.

C     msinfo:
C
C     0: both g and B are null
C     1: first convergence criterion is satisfied
C     2: second convergence criterion is satisfied
C     3: third convergence criterion is satisfied
C     5: maximum allowed number of iterations is achieved

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer col,dum,i,idx,iter,lin,lssinfo
      double precision b1n,d,delta2,geucn,ll,llant,ls,lsant,lu,luant,
     +        peucn,peucn2,ptz,rpeucn2,rzeucn2,sgn,tau,teucn2,tmp,ueucn2

C     LOCAL ARRAYS
      integer wi(nmax)
      double precision t(nmax),wd(nmax),z(nmax)

      delta2 = delta**2

      msinfo  = 0

      memfail = .false.

      call lssini(.false.,.false.,.true.)

C     step 1: initialize ls (lower bound on l) with max{-bii}, where bii 
C     are the elements of the diagonal of B
      
      ls = -bval(bdiag(1))
      
      do i = 2,n
          lin = bdiag(i)
          ls  = max( ls, -bval(lin) )
      end do
      
C     Calculate ||B||1, B sparse
      
      do i = 1,n
          wd(i) = 0.0d0
      end do     
      
      do i = 1,bnnz
          lin = blin(i)
          col = bcol(i)
          
          if ( ( lin .le. n ) .and. ( col .le. n ) ) then
              wd(col) = wd(col) + abs( bval(i) )
              if ( lin .ne. col ) then
                  wd(lin) = wd(lin) + abs( bval(i) )
              end if
          end if
      end do
      
      b1n = wd(1)
      do i = 2,n
          b1n = max( b1n, wd(i) )
      end do
      
C     step 2: initialize ll (lower bound on l) with 
C             max{0, ls, ||g||/delta - ||B||1}, where ||B||1 is the 
C             1-norm of the matrix B

      geucn = 0.0d0
      do i = 1,n
          geucn = geucn + g(i)**2
      end do
      geucn = sqrt( geucn )
      
      ll = (geucn / delta) - b1n 
      ll = max( 0.0d0, ll )
      ll = max( ls, ll )

C     step 3: initialize lu (upper bound on l) with ||g||/delta + ||B||1

      lu = (geucn / delta) + b1n

C     If the matrix is null, there is nothing to be done

      if ( ( abs( ll ) .le. macheps23 ) .and.
     +     ( abs( lu ) .le. macheps23 ) .and. 
     +     ( abs( ls ) .le. macheps23 ) ) then
          msinfo = 0
          go to 21
      end if
      
C     step 4: initialize iteration counter

      iter = 1

      call lssana(n,bnnz,blin,bcol,bval,bdiag,lssinfo)
      
      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if

C     step 5: safeguard of l (ensures that l is bigger than ll)
      
 5    continue

      l = max( l, ll )

C     step 6: safeguard of l (ensures that l is smaller than lu)

      l = min( l, lu )

C     step 7: safeguard of l

      if ( ( l .le. ls + macheps23 * max( abs( ls ), 1.0d0 ) )
     +     .and. ( iter .ne. 1 ) ) then
          l = max( 1.0d-3 * lu, sqrt( ll*lu ) )
      end if

C     step 8: try to use the Cholesky decomposition: (B +lI) = R^TR.
C             If the decomposition is successfull, R is stored in the 
C             upper triangular portion of B (including the diagonal) and 
C             pd is set to true.
C             If the decomposition fails, d and idx are set as explained 
C             before, pd is set to false, and the Euclidian-norm of u is 
C             calculated (see explanation of variable ueucn2)

      do i = 1,n
          wd(i) = l
      end do

      call lssfac(n,bnnz,blin,bcol,bval,bdiag,wd,idx,d,dum,lssinfo)
      chcnt = chcnt + 1
      
      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if

      if ( ( lssinfo .eq. 1 ) .or. ( lssinfo .eq. 2 ) ) then
          pd = .false.
      else ! lssinfo .eq. 0
          pd = .true.
      end if

C     In this case (B + lI) is not positive definite, and d and idx are 
C     calculated. Because p cannot be calculated (it is not possible to 
C     solve the system using Cholesky factorization), the values of l, 
C     ll and ls are updated, the iteration counter is increased and a
C     new iteration is started

      if ( .not. pd ) then

C         Print information (current iteration)

          if ( iprintinn .ge. 5 ) then
              write(*, 1000) iter
              write(*, 1010) ls,ll,lu,l
              write(*, 1070)
              
              write(10,1000) iter
              write(10,1010) ls,ll,lu,l
              write(10,1070) 
          end if
          
          llant = ll
          luant = lu
          lsant = ls
          
          if ( lu-l .le. macheps23 * max( abs( lu ),1.0d0 ) ) then
              lu = lu + macheps23 * max( abs( lu ),1.0d0 )
              l  = lu
          else
              call scalcu(n,bnnz,blin,bcol,bval,bdiag,l,idx,p,ueucn2,wd,
     +        memfail)
              chcnt = chcnt + 1
              
              if ( memfail ) then
                  go to 22
              end if
              
              ll = max( l, ll )
              ls = max( l + (d / ueucn2), ls )
              ll = max( ll, ls )
              l  = ls
          end if
          
          iter = iter + 1
          
C         Test whether the number of iterations is exhausted

          if ( iter .gt. maxit ) then
              msinfo = 5
              
              if ( iprintinn .ge. 5 ) then
                  write(*, 1090) 
                  write(10,1090) 
              end if
              
              go to 22
          end if
          
          go to 5
      end if
      
C     step 9: solve R^TRp = -g for p and calculate the squared 
C             Euclidian-norm of p
      
      do i = 1,n
          p(i) = -g(i)
      end do
      
      call lsssol(n,p)

C     Euclidian-norm of Rp = p^T R^TRp = p^T (-g) = - p^Tg
      
      rpeucn2 = 0.0d0
      do i = 1,n
          rpeucn2 = rpeucn2 - p(i) * g(i)
      end do
      
      peucn2 = 0.0d0
      do i = 1,n
          peucn2 = peucn2 + p(i)**2
      end do
      peucn = sqrt( peucn2 )
      
C     step 10: calculate z and tau, where tau * z is the approximation 
C              of the eigenvector associated with the smallest 
C              eigenvalue of B

      if ( peucn .lt. delta ) then

C        Calculate z

          call scalcz(n,wi,wd,z)
          
C         Calculate z Euclidian-norm

          tmp = 0.0d0
          do i = 1,n
              tmp = tmp + z(i)**2
          end do
          tmp = sqrt( tmp )
          
C         Divide z by its norm

          do i = 1,n
              z(i) = z(i) / tmp
          end do
          
C         Calculate the squared Euclidian-norm of the product Rz.         
C         Note that z^T R^T Rz = z^T (B + lI) z
          
          do i = 1,n
              wd(i) = 0.0d0
          end do
          
          do i = 1,bnnz
              lin = blin(i)
              col = bcol(i)
              
              if ( lin .eq. col ) then
                  wd(lin) = wd(lin) + (bval(i) + l) * z(col)
              else
                  wd(lin) = wd(lin) + bval(i) * z(col)
                  wd(col) = wd(col) + bval(i) * z(lin)
              end if
          end do
          
          rzeucn2 = 0.0d0
          do i = 1,n
              rzeucn2 = rzeucn2 + z(i) * wd(i)
          end do
          
C         Calculate tau

          ptz = 0.0d0
          do i = 1,n
              ptz = ptz + p(i) * z(i)
          end do
          
          if ( ptz .lt. 0.0d0 ) then
              sgn = -1.0d0
          else
              sgn =  1.0d0
          end if
          
          tmp = delta2 - peucn2
          
          tau = (ptz**2) + tmp
          tau = tmp / (ptz + sgn * sqrt( tau ))
          
      end if

C     Print informations (current iteration)

      if ( iprintinn .ge. 5 ) then
          write(*, 1000) iter
          write(*, 1010) ls,ll,lu,l
          write(*, 1020) peucn
          
          write(10,1000) iter
          write(10,1010) ls,ll,lu,l
          write(10,1020) peucn
          
          if ( peucn .lt. delta ) then
              write(*, 1030) sqrt( peucn2 + 2 * tau * ptz + tau ** 2 )
              write(10,1030) sqrt( peucn2 + 2 * tau * ptz + tau ** 2 )
          else
              write(*, 1030) peucn
              write(10,1030) peucn
          end if
          
          write(*, 1050) delta
          write(10,1050) delta
      end if
      
C     steps 11 and 12: update ll, lu and ls

      llant = ll
      luant = lu
      lsant = ls
      
      if ( peucn .lt. delta ) then
          lu = min( l, lu )
          ls = max( l - rzeucn2, ls )
      else
          ll = max( l, ll )
      end if

C     step 13: update ls when B + lI is not positive definite.
C              This was done right after the Cholesky decomposition 
C              failure

C     step 14: update ll

      ll = max( ll, ls )
      
C     step 15: convergence test

      if ( ( abs( l ) .le. eps ) .and. ( peucn .le. delta ) ) then
          msinfo = 1
          go to 21
      end if

C     step 16: second convergence test

      if ( abs( delta - peucn ) .le. sigma1*delta ) then
          msinfo = 2
      end if

C     step 17: convergence test for the hard case

      tmp = rpeucn2 + l * delta2
      tmp = max( sigma2, tmp )
      tmp = tmp * sigma1 * (2 - sigma1)
      
      if ( ( peucn .lt. delta ) .and. ( (rzeucn2*(tau**2) .le. tmp ) 
     +     .or. (lu - ll .le. eps) ) ) then     
          
          msinfo = msinfo + 3
          go to 20
      end if
      
      if ( msinfo .eq. 2 ) then
          go to 21
      end if
      
C     step 21: Calculate l to be used in the next iteration

      if ( ( abs( geucn ) .gt. eps ) .or. 
     +     ( ( l .le. ll + macheps23 * max( abs( ll ), 1.0d0 ) )
     +     .and. ( ls .le. ll ) ) ) then 
          
C         Solve R^T t = p for t and calculate t squared Euclidian-norm
         
          do i = 1,n
              t(i) = p(i)
          end do
          
          call lsssoltr('T',n,t)
          
          teucn2 = 0.0d0
          do i = 1,n
              teucn2 = teucn2 + t(i)**2
          end do
          
C         Update l using Newton's method update
        
          l = l + (peucn2 / teucn2) * ((peucn - delta) / delta)
      else
          l = ls
      end if
      
C     step 22: update iteration counter

      iter = iter + 1

C     Test whether the number of iterations is exhausted

      if ( iter .gt. maxit ) then
          msinfo = 5
          
          if ( iprintinn .ge. 5 ) then
              write(*, 1090) 
              write(10,1090) 
          end if
          
          go to 22
      end if      
      
C     step 23: start a new iteration

      if ( ( abs( llant-ll ) .le. 
     +     macheps23 * max( abs( ll ), 1.0d0 ) ) .and. 
     +     ( abs( luant-lu ) .le.
     +     macheps23 * max( abs( lu ), 1.0d0 ) ) .and. 
     +     ( abs( lsant-ls ) .le. 
     +     macheps23 * max( abs( ls ), 1.0d0 ) ) ) then
          
          ll = ll + macheps23 * max( abs( ll ), 1.0d0 )
          lu = lu - macheps23 * max( abs( lu ), 1.0d0 )
      end if
      
      go to 5

C     steps 18, 19 and 20:

C     The solution is given by p in 3 cases:
C     - if the first convergence criterion is satisfied;
C     - if only the second convergence criterion is satisfied;
C     - if both the second and the third convergence criteria are 
C       satisfied, but the squared Euclidian-norm of R(tau*z) is 
C       strictly bigger than l*(delta2 - peucn2)

C     The solution is given by p + tau*z when:
C     - just the third convergence criterion is satisfied;
C     - both the second and the third convergence criteria are 
C       satisfied, but the squared Euclidian-norm of R(tau*z) is smaller 
C       or equal to l*(delta2 - peucn2)

 20   tmp = (rzeucn2 * (tau**2)) - l * (delta2 - peucn2)
      
      if ( ( msinfo .eq. 3 ) .or. 
     +     ( ( msinfo .eq. 5 ) .and. ( tmp .le. 0.0d0 ) ) ) then
          
          peucn2 = 0.0d0
          do i = 1,n
              p(i) = p(i) + tau * z(i)
              peucn2 = peucn2 + p(i)**2
          end do
          peucn = sqrt( peucn2 )
          
          msinfo = 3
      else
          msinfo = 2
      end if
      
C     Print informations

 21   if ( iprintinn .ge. 5 ) then
          write(*, 1060) iter
          write(*, 1010) ls,ll,lu,l
          write(*, 1080) 
          write(*, 1040) peucn
          
          write(10,1060) iter
          write(10,1010) ls,ll,lu,l
          write(10,1080) 
          write(10,1040) peucn
      end if
      
 22   continue
      
C     Non-executable statements

 1000 format(/,10X,'More-Sorensen iteration: ',I7)
 1010 format(  10X,'ls = ',1P,D11.4,
     +         10X,'ll = ',1P,D11.4,
     +         10X,'lu = ',1P,D11.4,
     +         10X,'l  = ',1P,D11.4)
 1020 format(  10X,'Euclidian-norm of p: ',1P,D7.1)
 1030 format(  10X,'Euclidian-norm of p + z: ',1P,D7.1)
 1040 format(  10X,'Euclidian-norm of step: ',1P,D7.1)
 1050 format(  10X,'delta: ',1P,D7.1)
 1060 format(  10X,'Number of iterations: ',I7)
 1070 format(  10X,'Matrix is not positive definite!')

 1080 format(  10X,'Flag of More-Sorensen: ',
     +             'convergence criterion satisfied.')
 1090 format(  10X,'Flag of More-Sorensen: ',
     +             'maximum number of iterations achieved.')

      end

C     ******************************************************************
C     ******************************************************************
      
      subroutine scalcu(n,annz,alin,acol,aval,adiag,l,idx,u,ueucn2,wd,
     +memfail)
            
      implicit none

C     SCALAR ARGUMENTS
      logical memfail
      integer annz,idx,n
      double precision l,ueucn2
      
C     ARRAY ARGUMENTS
      integer acol(annz),adiag(n),alin(annz)
      double precision aval(annz),u(idx),wd(idx-1)
      
C     Solves a sparse linear system of the form (A + lI + ekek^Td)u = 0,
C     where A + lI is a definite positive matrix in R^{k x k}.
C     u is a vector of k positions with u(k) = 1.

C     LOCAL SCALARS
      integer col,i,lin,lssinfo
      
      if ( idx .eq. 1 ) then
          u(idx) = 1.0d0
          ueucn2 = 1.0d0
          return
      end if
      
      do i = 1,idx
          u(i) = 0.0d0
      end do
      
C     Permute columns and rows of A

      call lsspermind(annz,alin,acol)
      
C     Eliminate columns that have index greater than idx and define u
C     as idx-th column of A + lI
      
      do i = 1,annz
          col = acol(i)
          lin = alin(i)
          
          if ( ( col .eq. idx ) .and. ( lin .lt. idx ) ) then
              u(lin) = u(lin) - aval(i)
          else if ( ( lin .eq. idx ) .and. ( col .lt. idx ) ) then
              u(col) = u(col) - aval(i)
          end if
          
          if ( ( col .ge. idx ) .or. ( lin .ge. idx ) ) then
              acol(i) = col + n
              alin(i) = lin + n
          end if
      end do
      
C     Solve system (A + lI)x = u

      call lssini(.false.,.true.,.false.)
      
      call lsspermvec(n,adiag)
      
      do i = 1,idx-1
          wd(i) = l
      end do
      
      call lssafsol(idx-1,annz,alin,acol,aval,adiag,wd,u,lssinfo)
      
      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          memfail = .true.
          return

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          memfail = .true.
          return

      else if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          memfail = .true.
          return

      end if
      
      call lssunpermvec(n,adiag)
      
      call lssini(.false.,.false.,.true.)
      
C     Undo modifications in acol and alin
      
      do i = 1,annz
          col = acol(i)
          
          if ( col .gt. n ) then
              acol(i) = acol(i) - n
              alin(i) = alin(i) - n
          end if
      end do
      
      call lssunpermind(n,annz,alin,acol)
      
      u(idx) = 1.0d0
      
      ueucn2 = 0.0d0
      do i = 1,idx
          ueucn2 = ueucn2 + u(i)**2
      end do
      
      end
      
C     ******************************************************************
C     ******************************************************************      
      
      subroutine scalcz(n,rowind,rowval,z)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      integer rowind(n)
      double precision rowval(n),z(n)
      
C     Sparse implementation of the technique presented by Cline, Moler,
C     Stewart and Wilkinson to estimate the condition number of a upper
C     triangular matrix.

C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer col,i,k,rownnz
      double precision acs,acz,d,ek,rki,s,sm,w,wk,wkm
      double precision lssgetd

C     EXTERNAL FUNCTIONS
      external lssgetd

      call lsssetrow(n)
      
      do i = 1,n
          z(i) = 0.0d0
      end do
      
      acz = 0.0d0
      acs = 1.0d0
      ek  = 1.0d0
      
      do k = 1,n
          d = lssgetd(k)
          
          if ( abs( acs*z(k) ) .gt. macheps23 ) then
              ek = dsign(ek,-acs*z(k))
          end if
          
          if ( abs( ek - acs*z(k) ) .gt. abs( d ) ) then
              s   = abs( d ) / abs( ek - acs*z(k) ) 
              acs = acs * s
              ek  = s * ek
          end if
          
          wk  =  ek - acs * z(k)
          wkm = -ek - acs * z(k)
          s   = abs( wk )
          sm  = abs( wkm )
          
          if ( d .eq. 0.0d0 ) then 
              wk  = 1.0d0
              wkm = 1.0d0
          else
              wk  = wk / d
              wkm = wkm / d
          end if
          
          if ( k .eq. n ) then
              go to 10
          end if
          
          call lssgetrow(n,k,rownnz,rowind,rowval)
          
          sm = sm + acs * acz
          
          do i = 1,rownnz
              
              col = rowind(i)
              if ( col .gt. k ) then
                  rki    = rowval(i)
                  sm     = sm - abs( acs * z(col) )
                  sm     = sm + abs( acs * z(col) + wkm * rki )
                  acz    = acz - abs( z(col) )
                  acz    = acz + abs( z(col) + (wk/acs) * rki )
                  z(col) = z(col) + (wk/acs) * rki
              end if
          end do
          
          s = s + acs * acz
          
          if ( s .lt. 
     +         sm - macheps23 * max( abs( sm ), 1.0d0 ) ) then
              w = wkm - wk
              wk = wkm
              
              do i = 1,rownnz
                  
                  col = rowind(i)
                  if ( col .gt. k ) then
                      rki    = rowval(i)
                      acz    = acz - abs( z(col) )
                      acz    = acz + abs( z(col) + (w/acs) * rki )
                      z(col) = z(col) + (w/acs) * rki
                  end if
              end do
              
          end if
          acz = acz - abs( z(k+1) )
 10       z(k) = wk/acs
      end do
      
C     Divide z by its 1-norm to avoid overflow
      
      s = 0.0d0
      do i = 1,n
          s = s + abs( z(i) )
      end do
      
      do i = 1,n
          z(i) = z(i) / s
      end do
      
C     Solve Rz = y

      call lsssoltr(' ',n,z)

      end

C     ******************************************************************
C     ******************************************************************

C     moresor:
C
C     Method to minimize sparse quadratic functions subjected to 
C     ||w|| <= delta
C
C     minimize    psi(w) = 1/2 w^TBw + g^Tw 
C     subject to  ||w|| <= delta
C
C     Method described in "Computing a trust region step", by More and 
C     Sorensen.
C
C     The main ideia of this method is to find a positive scalar \mslamb 
C     that is a zero of the function
C
C     phi(\mslamb) = 1/||p|| - 1/delta,
C      
C     where p is the solution of the linear system
C
C     (B + \mslamb I)p = -g.
C
C     Note that symmetric matrix B, vector g and positive real number 
C     delta are presented in the minimization problem above. I is the 
C     identity matrix.
C
C     The method used to find the zero of that function is basically the 
C     Newton method to find roots.
C
C     On Entry
C
C     n        integer
C              dimension
C
C     g        double precision g(n)
C              vector used to define the quadratic function
C
C     bnnz     integer
C              number of nonzero elements of B
C
C     blin     integer blin(bnnz)
C              row indices of nonzero elements of B
C
C     bcol     integer bcol(bnnz)
C              column indices of nonzero elements of B
C
C     bval     double precision bval(bnnz)
C              nonzero elements of B, that is,
C              B(blin(i),bcol(i)) = bval(i). Since B is symmetric, just
C              one element B(i,j) or B(j,i) must appear in its sparse
C              representation. If more than one pair corresponding to
C              the same position of B appears in the sparse
C              representation, the multiple entries will be summed.
C              If any blin(i) or bcol(i) is out of range, the entry will
C              be ignored
C     
C     bdiag    integer bdiag(n)
C              indices of diagonal elements of B in blin, bcol and bval
C
C     delta    double precision
C              trust-region radius
C     
C     sigma1   double precision
C              allowed error for convergence criteria 1 and 2
C
C     sigma2   double precision
C              allowed error for convergence criteria 3 (hard case)
C
C     eps      double precision
C              allowed error
C
C     maxit    integer
C              maximum number of allowed iterations
C     
C     l        double precision
C              initial value for mslamb
C
C     On Return
C      
C     l        double precision
C              value that gives p as a solution to the minimization 
C              problem, because p is also solution to 
C              (B + l I)p = -g
C
C     pd       logical
C              set to true if the last Cholesky decomposition is 
C              successfull
C
C     p        double precision p(n)
C              solution to problem
C              minimize     psi(w)
C              subjected to ||w|| <= delta
C
C     chcnt    integer
C              number of Cholesky decompositions
C
C     memfail  logical
C              true iff linear solver failed because of lack of memory
C
C     msinfo   integer
C              stores which convergence criteria was satisfied:
C
C              0 = both g and B are null;
C
C              1 = first convergence criterion is satisfied;
C
C              2 = second convergence criterion is satisfied;
C
C              3 = third convergence criterion is satisfied
C
C              5 = maximum allowed number of iterations is achieved.

C     ******************************************************************
C     ******************************************************************

C     scalcu:
C
C     Solve a sparse linear system of the form (A + lI + ekek^Td)u = 0,
C     where A + lI is a definite positive matrix in R^{k x k}.
C     u is a vector of k positions with u(k) = 1.
C
C     On Entry
C      
C     n        integer
C              dimension of A
C
C     annz     integer
C              number of nonzero elements of A
C      
C     alin     integer alin(annz)
C              row indices of nonzero elements of A
C
C     acol     integer acol(annz)
C              column indices of nonzero elements of A
C
C     aval     integer aval(annz)
C              nonzero elements of A, that is,
C              A(alin(i),acol(i)) = aval(i). Since A is symmetric, just
C              one element A(i,j) or A(j,i) must appear in its sparse
C              representation. If more than one pair corresponding to
C              the same position of A appears in the sparse
C              representation, the multiple entries will be summed.
C              If any alin(i) or acol(i) is out of range, the entry will
C              be ignored
C     
C     adiag    integer adiag(n)
C              indices of diagonal elements of A in alin, acol and aval
C
C     l        double precision
C              used to compute A + lI
C      
C     idx      integer
C              index k 
C
C     On Return
C      
C     u        double precision u(n)
C              system solution
C
C     ueucn2   double precision
C              u squared Euclidian-norm
C
C     memfail  logical
C              true iff linear solver failed because of lack of memory

C     ******************************************************************
C     ******************************************************************      
      
C     scalcz:
C
C     Sparse implementation of the technique presented by Cline, Moler,
C     Stewart and Wilkinson to estimate the condition number of a matrix.
C     This technique is used by the More-Sorensen method to calculate an
C     approximation to the eigenvector associated to the smallest 
C     eigenvalue of a matrix B (\mslamb_1). 
C     In this technique, when \mslamb approaches -\mslamb_1, \|Rz\| 
C     approaches 0. This insures that z is an approximation to the 
C     wanted eigenvector.
C     Basically, it solves R^Ty = e, choosing e(k) as 1 or -1 (whatever 
C     gives maximum local growth of y). Then, it solves Rz = y.
C     Note that R is a sparse matrix given by P^T D^0.5 L^T P (obtained
C     applying subroutine MA27BD).
C      
C     On Entry
C
C     n        integer      
C              dimension of R
C      
C     rowind   integer rowind(n)
C     rowval   double precision rowval(n)
C              working arrays
C      
C     On Return
C      
C     z        double precision z(n)
C              approximation of the eigenvector of B associated to
C              \mslamb_1
C     ******************************************************************
C     ******************************************************************

      subroutine spgls(n,x,l,u,m,lambda,rho,equatn,linear,f,g,lamspg,
     +xp,fp,alpha,d,evalaldim,lsinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,lsinfo,m,n
      double precision alpha,f,fp,lamspg

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision d(n),g(n),l(n),lambda(m),rho(m),u(n),x(n),xp(n)

C     SUBROUTINE ARGUMENTS
      external evalaldim

C     This subroutine computes a line search in the Spectral Projected 
C     Gradient direction.
C
C     lsinfo:
C
C     0: Armijo satisfied
C     1: Small step with functional value similar to the current one
C     2: Unbounded objective function?
C     3: Too small backtracking step. Wrong gradient?

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/

C     LOCAL SCALARS
      integer i
      double precision dsupn,gtd,xsupn

C     ------------------------------------------------------------------
C     Compute search direction, directional derivative, dsupn, xsupn and 
C     first trial
C     ------------------------------------------------------------------

      gtd   = 0.0d0
      dsupn = 0.0d0
      xsupn = 0.0d0
      do i = 1,n
          d(i)  = - lamspg * g(i)
          xp(i) = x(i) + d(i)
          if ( xp(i) .lt. l(i) .or. xp(i) .gt. u(i) ) then
              xp(i) = max( l(i), min( xp(i), u(i) ) )
              d(i)  = xp(i) - x(i)
          end if
          gtd   = gtd + g(i) * d(i)
          dsupn = max( dsupn, abs( d(i) ) )
          xsupn = max( xsupn, abs( x(i) ) )
      end do

      if ( iprintinn .ge. 6 ) then
          write(* ,100) xsupn,lamspg,dsupn
          write(10,100) xsupn,lamspg,dsupn
      end if
      
      call evalaldim(n,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

      alpha = 1.0d0
      
      if ( iprintinn .ge. 6 ) then
          write(*, 110) alpha,fp,fcnt
          write(10,110) alpha,fp,fcnt
      end if
      
C     ==================================================================
C     Backtracking
C     ==================================================================

      call backtracking(n,x,m,lambda,rho,equatn,linear,f,xsupn,d,gtd,
     +dsupn,alpha,fp,xp,evalaldim,lsinfo,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     End of backtracking
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 100  format(/,5X,'SPG Line search (xsupn = ',1P,D7.1,1X,'SPGstep= ',
     +             1P,D7.1,1X,'dsupn = ',1P,D7.1,')')
 110  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)

      end
C     *****************************************************************
C     *****************************************************************

      subroutine tnls(nind,x,l,u,m,lambda,rho,equatn,linear,f,g,amax,d,
     +rbdnnz,rbdind,rbdtype,xp,fp,gp,lsinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,lsinfo,m,nind,rbdnnz
      double precision amax,f,fp

C     ARRAY ARGUMENTS
      integer rbdind(nind)
      character rbdtype(nind)
      logical equatn(m),linear(m)
      double precision d(nind),g(*),gp(*),l(nind),lambda(m),rho(m),
     +        u(nind),x(*),xp(*)

C     This subroutine computes a line search in direction d.
C
C     lsinfo:
C
C     At the first trial:
C
C     5: x + amax d is at the boundary and f(x + amax d) is smaller than f

C     Extrapolation:
C
C     2: Unbounded objective function?
C     4: beta-condition holds. No extrapolation is done
C     6: Maximum number of extrapolations reached
C     7: Similar consecutive projected points
C     8: Not-well-defined objective function
C     9: Functional value increases

C     Backtracking:
C
C     0: Armijo satisfied
C     1: Small step with functional value similar to the current one
C     2: Unbounded objective function?
C     3: Too small backtracking step. Wrong gradient?

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      logical boundary
      integer i
      double precision alpha,dsupn,gtd,xsupn

C     EXTERNAL SUBROUTINES
      external calcal

C     ==================================================================
C     ==================================================================

C     Test Armijo condition and beta condition. Decide between accept 
C     the first trial, extrapolate or backtrack. 

C     ==================================================================
C     ==================================================================

C     ------------------------------------------------------------------
C     Compute directional derivative, dsupn and xsupn
C     ------------------------------------------------------------------

      gtd = 0.0d0
      dsupn = 0.0d0
      xsupn = 0.0d0
      do i = 1,nind
          gtd = gtd + g(i) * d(i)
          dsupn = max( dsupn, abs( d(i) ) )
          xsupn = max( xsupn, abs( x(i) ) )
      end do

      if ( iprintinn .ge. 6 ) then
          write(*, 100) xsupn,amax,dsupn
          write(10,100) xsupn,amax,dsupn
      end if

C     ------------------------------------------------------------------
C     Compute first trial (projected point)
C     ------------------------------------------------------------------

      alpha = 1.0d0

      boundary = .false.
      do i = 1,nind
          xp(i) = x(i) + d(i)

          if ( xp(i) .lt. l(i) .or. xp(i) .gt. u(i) ) then
              boundary = .true.
              xp(i) = max( l(i), min( xp(i), u(i) ) )
          end if
      end do

      if ( amax .eq. 1.0d0 ) then
          boundary = .true.

          do i = 1,rbdnnz
              if ( rbdtype(i) .eq. 'L' ) then
                  xp(rbdind(i)) = l(rbdind(i))
              else if ( rbdtype(i) .eq. 'U' ) then
                  xp(rbdind(i)) = u(rbdind(i))
              end if
          end do
      end if

      call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

      if ( .not. boundary ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 110) fp,fcnt
              write(10,110) fp,fcnt
          end if

      else
          if ( iprintinn .ge. 6 ) then
              write(*, 120) fp,fcnt
              write(10,120) fp,fcnt
          end if
      end if

C     ------------------------------------------------------------------
C     The first trial is an interior point.
C     ------------------------------------------------------------------

      if ( .not. boundary ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 140)
              write(10,140)
          end if

C         Armijo condition holds.

          if ( fp .le. f + alpha * gamma * gtd ) then

              if ( iprintinn .ge. 6 ) then
                  write(*, 150)
                  write(10,150)
              end if

              go to 1000

          end if

C         Armijo condition does not hold. We will do backtracking.  

          if ( iprintinn .ge. 6 ) then
              write(* ,180)
              write(10,180)
          end if

          go to 2000

      end if

C     ------------------------------------------------------------------
C     The first trial is at the boundary.
C     ------------------------------------------------------------------

      if ( iprintinn .ge. 6 ) then
          write(*, 190) 
          write(10,190) 
      end if

C     Function value is smaller than at the current point. We will 
C     extrapolate.

      if ( fp .lt. f ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 200) 
              write(10,200) 
          end if

          go to 1000

      end if

C     Discard the projected point and consider x + amax d

      if ( iprintinn .ge. 6 ) then
          write(*, 210) 
          write(10,210) 
      end if

      alpha = amax

      do i = 1,nind
          xp(i) = x(i) + alpha * d(i)
      end do

      do i = 1,rbdnnz
          if ( rbdtype(i) .eq. 'L' ) then
              xp(rbdind(i)) = l(rbdind(i))
          else if ( rbdtype(i) .eq. 'U' ) then
              xp(rbdind(i)) = u(rbdind(i))
          end if
      end do

      call calcal(nind,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

      if ( iprintinn .ge. 6 ) then
          write(*, 130) alpha,fp,fcnt
          write(10,130) alpha,fp,fcnt
      end if

C     Function value is smaller than or equal to (or even just a little 
C     bit greater than) at the current point. Line search is over.

      if ( fp .lt. f + macheps23 * abs(f) ) then

          if ( iprintinn .ge. 6 ) then
              write(*, 220) 
              write(10,220) 
          end if

          call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
          if ( inform .lt. 0 ) return

          lsinfo = 5

          if ( iprintinn .ge. 6 ) then
              write(*, 900)
              write(10,900)
          end if

          return

      end if

C     Function value is greater than at the current point. We will 
C     do backtracking.

      if ( iprintinn .ge. 6 ) then
          write(*, 230) 
          write(10,230) 
      end if

      go to 2000

C     ==================================================================
C     ==================================================================

C     Extrapolation

C     ==================================================================
C     ==================================================================

 1000 continue

      call extrapolation(nind,x,l,u,m,lambda,rho,equatn,linear,g,xp,fp,
     +gp,d,alpha,amax,rbdnnz,rbdind,rbdtype,fmin,beta,etaext,maxextrap,
     +lsinfo,inform)

      return

C     ==================================================================
C     ==================================================================

C     Backtracking

C     ==================================================================
C     ==================================================================

 2000 continue

      call backtracking(nind,x,m,lambda,rho,equatn,linear,f,xsupn,d,gtd,
     +dsupn,alpha,fp,xp,calcal,lsinfo,inform)
      if ( inform .lt. 0 ) return

      call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     End of backtracking
C     ==================================================================

C     NON-EXECUTABLE STATEMENTS

 100  format(/,5X,'TN Line search (xsupn = ',1P,D7.1,', amax = ',
     +             1P,D7.1,', dsupn = ',1P,D7.1,')')
 110  format(  5X,'Unitary step    F = ',1P,D24.16,' FE = ',I7)
 120  format(  5X,'Projected point F = ',1P,D24.16,' FE = ',I7)
 130  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 140  format(  5X,'The first trial is an interior point.')
 150  format(  5X,'Armijo condition holds.')
 180  format(  5X,'Armijo condition does not hold. We will backtrack.')
 190  format(  5X,'The first trial is at the boundary.')
 200  format(  5X,'Function value at the boundary is smaller than at ',
     +            'the current point.',/,5X,'We will extrapolate.')
 210  format(  5X,'Discarding projected point. We will now consider x ',
     +            '+ amax d.')
 220  format(  5X,'Function value at the boundary is smaller than or ',
     +            'equal to than at the',/,5X,'current point. Line ',
     +            'search is over.')
 230  format(  5X,'Function value at the boundary is greater than at ',
     +            'the current point.',/,5X,'We will backtrack.')
 900  format(  5X,'Flag of TN Line search: First trial accepted.')

      end
      subroutine backtracking(dim,x,m,lambda,rho,equatn,linear,f,xsupn,
     +d,gtd,dsupn,alpha,fp,xp,evalaldim,btinfo,inform)

C     SCALAR ARGUMENTS
      integer btinfo,dim,inform,m
      double precision alpha,dsupn,f,fp,gtd,xsupn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision d(dim),lambda(m),rho(m),x(dim),xp(*)

C     SUBROUTINE ARGUMENTS
      external evalaldim

C     Backtracking with quadratic interpolation.
C
C     btinfo:
C
C     0: Armijo satisfied
C     1: Small step with functional value similar to the current one
C     2: Unbounded objective function?
C     3: Too small backtracking step. Wrong gradient?

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      logical smallstep
      integer i,interp
      double precision atmp

      interp = 0

 2010 continue

C     Test Armijo condition

      if ( fp .le. f + alpha * gamma * gtd ) then

C         Finish backtracking with the current point
   
          btinfo = 0

          if ( iprintinn .ge. 6 ) then
              write(*, 900)
              write(10,900)
          end if

          return

      end if

C     Test if we obtained a functional value similar to the current one
C     associated to a very small step
      
      if ( alpha * dsupn .le. macheps23 * xsupn .and.
     +     fp - f        .le. macheps23 * abs(f) ) then

C         Finish backtracking with the current point
          
          btinfo = 1
          
          if ( iprintinn .ge. 6 ) then
              write(*, 910)
              write(10,910)
          end if
          
          return
          
      end if

C     Test f going to -inf

      if ( fp .le. fmin ) then

C         Finish backtracking with the current point

          btinfo = 2

          if ( iprintinn .ge. 6 ) then
              write(*, 920)
              write(10,920)
          end if

          return

      end if

C     Compute new step

      interp = interp + 1

      atmp = ( - gtd * alpha ** 2 ) / 
     +       ( 2.0d0 * ( fp - f - alpha * gtd ) )

      if ( atmp .ge. sigma1 * alpha .and. 
     +     atmp .le. sigma2 * alpha ) then
          alpha = atmp
      else
          alpha = alpha / etaint
      end if

C     Compute new trial point

      do i = 1,dim
          xp(i) = x(i) + alpha * d(i)
      end do

      call evalaldim(dim,xp,m,lambda,rho,equatn,linear,fp,inform)
      if ( inform .lt. 0 ) return

C     Print information of this iteration

      if ( iprintinn .ge. 6 ) then
          write(*, 110) alpha,fp,fcnt
          write(10,110) alpha,fp,fcnt
      end if

C     Test whether the step is too small

      smallstep = .true.
      do i = 1,dim
          if ( alpha * abs( d(i) ) .gt. 
     +         macheps * max( 1.0d0, abs( x(i) ) ) ) then
              smallstep = .false.
          end if
      end do

      if ( interp .ge. mininterp .and. smallstep .and. 
     +     fp - f .le. macheps23 * abs(f) ) then

C         Finish backtracking with the current point

          btinfo = 3

          if ( iprintinn .ge. 6 ) then
              write(*, 930)
              write(10,930)
          end if
  
          return

      end if

      go to 2010

C     NON-EXECUTABLE STATEMENTS

 110  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 900  format(  5X,'Flag of backtracking: Armijo condition holds.')
 910  format(  5X,'Flag of backtracking: Small step with similar ',
     +            'functional value.')
 920  format(  5X,'Flag of backtracking: Unbounded objective function?')
 930  format(  5X,'Flag of backtracking: Too small backtracking step.')

      end
C     *****************************************************************
C     *****************************************************************

      subroutine extrapolation(nind,x,l,u,m,lambda,rho,equatn,linear,g,
     +xp,fp,gp,d,alpha,amax,rbdnnz,rbdind,rbdtype,fmin,beta,etaext,
     +maxextrap,extinfo,inform)

C     SCALAR ARGUMENTS
      integer extinfo,inform,m,maxextrap,nind,rbdnnz
      double precision alpha,amax,beta,etaext,fmin,fp

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      character rbdtype(nind)
      integer rbdind(nind)
      double precision d(nind),g(*),gp(*),l(nind),lambda(m),rho(m),
     +        u(nind),x(*),xp(*)
      
C     Performs extrapolation.
C
C     extinfo:
C
C     0: Success
C     2: Unbounded objective function?
C     4: beta-condition holds. No extrapolation is done
C     6: Maximum number of extrapolations reached
C     7: Similar consecutive projected points
C     8: Not-well-defined objective function
C     9: Functional value increases

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      logical projected,samep
      integer extrap,i
      double precision atmp,fbext,ftmp,gptd,gtd

C     LOCAL ARRAYS
      double precision xbext(nmax),xtmp(nmax)

      extinfo = 0

      extrap  = 0

C     Compute directional derivative

      gtd = 0.0d0
      do i = 1,nind
          gtd = gtd + g(i) * d(i)
      end do
      
      call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
      if ( inform .lt. 0 ) return
      
      gptd = 0.0d0
      do i = 1,nind
          gptd = gptd + gp(i) * d(i)
      end do

C     Beta condition holds. No extrapolation is done.

      if ( gptd .ge. beta * gtd ) then
          
          if ( iprintinn .ge. 6 ) then
              write(*, 110)
              write(10,110)
          end if
          
          extinfo = 4
          
          return
          
      end if

C     Beta condition does not holds. We will extrapolate.

      if ( iprintinn .ge. 6 ) then
          write(* ,120)
          write(10,120)
      end if
      
C     Save f and x before extrapolation to return in case of a
C     not-well-defined objective function at an extrapolated point

      fbext = fp

      do i = 1,nind
          xbext(i) = xp(i)
      end do

 1010 continue

C     Test f going to -inf

      if ( fp .le. fmin ) then

C         Finish the extrapolation with the current point

          if ( extrap .ne. 0 ) then
              call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
              if ( inform .lt. 0 ) return
          end if

          extinfo = 2

          if ( iprintinn .ge. 6 ) then
              write(*, 910)
              write(10,910)
          end if

          return

      end if

C     Test if the maximum number of extrapolations was exceeded

      if ( extrap .ge. maxextrap ) then

C         Finish the extrapolation with the current point

          if ( extrap .ne. 0 ) then
              call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
              if ( inform .lt. 0 ) return
          end if

          extinfo = 6

          if ( iprintinn .ge. 6 ) then
              write(*, 930)
              write(10,930)
          end if

          return

      end if

C     Chose new step 

      extrap = extrap + 1

      if ( alpha .lt. amax .and. etaext * alpha .gt. amax ) then
          atmp = amax
      else
          atmp = etaext * alpha
      end if

C     Compute new trial point

      do i = 1,nind
          xtmp(i) = x(i) + atmp * d(i)
      end do

      if ( atmp .eq. amax ) then
          do i = 1,rbdnnz
              if ( rbdtype(i) .eq. 'L' ) then
                  xtmp(rbdind(i)) = l(rbdind(i))
              else if ( rbdtype(i) .eq. 'U' ) then
                  xtmp(rbdind(i)) = u(rbdind(i))
              end if
          end do
      end if

C     Project

      if ( atmp .gt. amax ) then

          projected = .false.
          do i = 1,nind
              if ( xtmp(i) .lt. l(i) .or. xtmp(i) .gt. u(i) ) then
                  projected = .true.
                  xtmp(i) = max( l(i), min( xtmp(i), u(i) ) )
              end if
          end do

      end if

C     Test if this is not the same point as the previous one. This test 
C     is performed only when xtmp is in fact a projected point.

      if ( projected ) then

          samep = .true.
          do i = 1,nind
              if ( abs( xtmp(i) - xp(i) ) .gt. 
     +             macheps * max( 1.0d0, abs( xp(i) ) ) ) then
                  samep = .false.
              end if
          end do

          if ( samep ) then

C             Finish the extrapolation with the current point

              call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
              if ( inform .lt. 0 ) return

              extinfo = 7

              if ( iprintinn .ge. 6 ) then
                  write(*, 940)
                  write(10,940)
              end if

              return

          end if

      end if

C     Evaluate function

      call calcal(nind,xtmp,m,lambda,rho,equatn,linear,ftmp,inform)

C     If the objective function is not well defined in an extrapolated 
C     point, we discard all the extrapolated points and return to a 
C     safe region (to the last point before the extrapolation)

      if ( inform .lt. 0 ) then

          fp = fbext

          do i = 1,nind
              xp(i) = xbext(i)
          end do

          call csetp(nind,xp)

          call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
          if ( inform .lt. 0 ) return

          extinfo = 8

          if ( iprintinn .ge. 6 ) then
              write(*, 950)
              write(10,950)
          end if

          return

      end if

C     Print information of this iteration

      if ( iprintinn .ge. 6 ) then
          write(*, 100) atmp,ftmp,fcnt
          write(10,100) atmp,ftmp,fcnt
      end if

C     If the functional value decreases then set the current point and 
C     continue the extrapolation

      if ( ftmp .lt. fp ) then

          alpha = atmp

          fp = ftmp

          do i = 1,nind
              xp(i) = xtmp(i)
          end do

          go to 1010

      end if

C     If the functional value does not decrease then discard the last 
C     trial and finish the extrapolation with the previous point

      call csetp(nind,xp)

      call calcnal(nind,xp,m,lambda,rho,equatn,linear,gp,inform)
      if ( inform .lt. 0 ) return

      extinfo = 9

      if ( iprintinn .ge. 6 ) then
          write(*, 960)
          write(10,960)
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(  5X,'Alpha = ',1P,D7.1,' F = ',1P,D24.16,' FE = ',I7)
 110  format(  5X,'Beta condition also holds. ',
     +            'No extrapolation is done.')
 120  format(  5X,'Beta condition does not hold. We will extrapolate.')

 910  format(  5X,'Flag of Extrapolation: Unbounded objective ',
     +            'function?')
 930  format(  5X,'Flag of Extrapolation: Maximum of consecutive ',
     +            'extrapolations reached.')
 940  format(  5X,'Flag of Extrapolation: Very similar consecutive ',
     +            'projected points.')
 950  format(  5X,'Flag of Extrapolation: Not-well-defined objective ',
     +            'function in an extrapolated point.')
 960  format(  5X,'Flag of Extrapolation: Functional value increased ',
     +            'when extrapolating.')

      end
C     ******************************************************************
C     ******************************************************************

      subroutine newtonkkt(n,xo,lo,uo,m,lambdao,equatno,linearo,epsfeas,
     +epsopt,f,cnorm,nlnorm,iter,msqiter,accinfo,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer accinfo,inform,iter,m,msqiter,n
      double precision cnorm,epsfeas,epsopt,f,nlnorm

C     ARRAY ARGUMENTS
      logical equatno(m),linearo(m)
      double precision lo(n),lambdao(m),uo(n),xo(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     accinfo:
C
C     0: KKT system solved.
C     1: Ignored constraints were violated.
C     2: After correcting the Lagrange multipliers signal, optimality
C        was lost.
C     3: Maximum number of iterations reached.
C     4: Newton seems to be diverging.
C     5: Singular Jacobian.
C     6: Insufficient space to store the KKT linear system.
C     7: Insufficient double precision working space for linear solver.
C     8: Insufficient integer working space for linear solver.

C     PARAMETERS
      integer mtotmax,ntotmax
      parameter ( mtotmax  = mmax + 2 * nmax )
      parameter ( ntotmax  = mmax + 3 * nmax )

C     LOCAL SCALARS
      integer cind,col,divit,fun,hlnnz,i,itmp,j,jcnnz,k,lssinfo,lin,
     +        maxit,minp,mnop,mtot,nbds,nineq,ninn,nneigv,nnon,nsys,
     +        ntot,pind,sind,var,vind
      double precision cnorm0,cnormprev,cnormref,epsact,epsadd,nlnorm0,
     +        nlnormprev,nlnormref,pval,sgnvio,val

C     LOCAL ARRAYS
      logical equatn(mtotmax),linear(mmax)
      character constt(2*nmax),status(nmax)
      integer consti(2*nmax),hdiag(ntotmax+mtotmax),inn(ntotmax),
     +        inp(mtotmax),jcfun(jcnnzmax),slaind(mtotmax)
      double precision adddiag(ntotmax+mtotmax),b(ntotmax+mtotmax),
     +        l(ntotmax),lambda(mtotmax),nl(ntotmax),r(mtotmax),
     +        t(nmax),u(ntotmax),x(ntotmax)

C     ==================================================================
C     PRESENTATION
C     ==================================================================

      if ( iprintout .ge. 1 ) then
          write(* ,8000)
          write(10,8000)
      end if

C     ==================================================================
C     INITIALIZE
C     ==================================================================

      iter    =  0
      divit   =  0
      maxit   = 10
      msqiter =  0

      epsact = sqrt( epsfeas )
      epsadd = macheps12

      call lssini(sclsys,.true.,.false.)

C     ==================================================================
C     SET INITIAL POINT
C     ==================================================================

      do i = 1,n
          l(i) = lo(i)
          u(i) = uo(i)
          x(i) = max( lo(i), min( xo(i), uo(i) ) )
          if ( lo(i) .eq. uo(i) ) then
              write(*,*) 'There is a fixed variable: ',i,lo(i),xo(i)
          end if
      end do

      do j = 1,m
          lambda(j) = lambdao(j)
          equatn(j) = equatno(j)
          linear(j) = linearo(j)
      end do

C     ==================================================================
C     COMPUTE CONSTRAINTS AND GRADIENT OF THE LAGRANGIAN
C     ==================================================================

      call ssetp(n,x)

C     Objective function and constraints

      call sevalobjc(n,x,f,m,r,inform)
      if ( inform .lt. 0 ) return

C     Gradient of the Lagrangian (and save Jacobian of the constraints)

      call sevalnl(n,x,m,lambda,equatn,linear,nl,inform)
      if ( inform .lt. 0 ) return

C     ==================================================================
C     STRUCTURES FOR NEW VARIABLES AND (BOUND) CONSTRAINTS
C     SET SLACKS VARIABLES VALUES
C     ==================================================================

C     Relate constraints to slacks

      nineq = 0
      do j = 1,m
          if ( .not. equatn(j) ) then
              nineq     = nineq + 1
              sind      = n + nineq
              slaind(j) = sind
              x(sind)   = sqrt( 2.0d0 * max( 0.0d0, - r(j) ) )
          else
              slaind(j) = 0
          end if
      end do

      nbds = 0
      do i = 1,n
          if ( l(i) .gt. - 1.0d+20 ) then
              nbds = nbds + 1
              constt(nbds) = 'L'
              consti(nbds) =  i
              cind         = m + nbds
              sind         = n + nineq + nbds
              slaind(cind) = sind
              r(cind)      = l(i) - x(i)
              x(sind)      = sqrt( 2.0d0 * max( 0.0d0, - r(cind) ) )
              r(cind)      = r(cind) + 0.5d0 * x(sind) ** 2
              equatn(cind) = .true.
          end if

          if ( u(i) .lt. 1.0d+20 .and. u(i) .ne. l(i) ) then
              nbds = nbds + 1
              constt(nbds) = 'U'
              consti(nbds) =  i
              cind         = m + nbds
              sind         = n + nineq + nbds
              slaind(cind) = sind
              r(cind)      = x(i) - u(i)
              x(sind)      = sqrt( 2.0d0 * max( 0.0d0, - r(cind) ) )
              r(cind)      = r(cind) + 0.5d0 * x(sind) ** 2
              equatn(cind) = .true.
          end if
      end do

      mtot = m + nbds
      ntot = n + nineq + nbds

C     ==================================================================
C     MAIN LOOP
C     ==================================================================

 100  continue

C     ==================================================================
C     SET ACTIVE CONSTRAINTS AND VARIABLES
C     ==================================================================

C     A fixed variable will be fixed forever.

      minp = 0
      mnop = 0

      ninn = 0
      nnon = 0

C     Set active variables

      do i = 1,n
           if ( x(i) .le. l(i) ) then
               status(i) = 'L'
               x(i) = l(i)

               nnon   = nnon + 1
               inn(i) = ntot + 1 - nnon

           else if ( x(i) .ge. u(i) ) then

               status(i) = 'U'
               x(i) = u(i)

               nnon   = nnon + 1
               inn(i) = ntot + 1 - nnon

           else
              status(i) = 'F'
              ninn   = ninn + 1
              inn(i) = ninn
           end if
      end do

C     Set active (regular) constraints and their slacks

      do j = 1,m
          if ( equatn(j) ) then
C             Active equality constraint
              minp   = minp + 1
              inp(j) = minp
          else
              sind = slaind(j)

              if ( r(j) .ge. - epsact ) then
C                 Active inequality constraint and its slack
                  minp      = minp + 1
                  inp(j)    = minp

                  ninn      = ninn + 1
                  inn(sind) = ninn
              else
C                 Deactivate inequality constraint and its slack
                  mnop      = mnop + 1
                  inp(j)    = mtot + 1 - mnop

                  nnon      = nnon + 1
                  inn(sind) = ntot + 1 - nnon
              end if
          end if
      end do

C     Set active bound constraints and their slacks

      do i = 1,nbds
          cind = m + i
          vind = consti(i)
          sind = slaind(cind)

          if ( status(vind) .eq. 'F' ) then
              minp      = minp + 1
              inp(cind) = minp

              ninn      = ninn + 1
              inn(sind) = ninn

          else
              mnop      = mnop + 1
              inp(cind) = mtot + 1 - mnop

              nnon      = nnon + 1
              inn(sind) = ntot + 1 - nnon

              if ( constt(i) .eq. status(vind) ) then
                  x(sind) = 0.0d0
              else
                  x(sind) = sqrt( 2.0d0 * ( u(vind) - l(vind) ) )
              end if
          end if
      end do

      nsys = ninn + minp

C     ==================================================================
C     SAVE NORMS TO CHECK IMPROVEMENT
C     ==================================================================

      if ( iter .eq. 0 ) then
          nlnormprev = bignum
          cnormprev  = bignum
      else
          nlnormprev = nlnorm
          cnormprev  = cnorm
      end if

C     ==================================================================
C     COMPUTE OBJECTIVE FUNCTION AND CONSTRAINTS
C     ==================================================================

      call ssetp(n,x)

C     Objective function and constraints

      call sevalobjc(n,x,f,m,r,inform)
      if ( inform .lt. 0 ) return

C     Add slacks effect

C     Regular constraints

      do j = 1,m
          sind = slaind(j)
          if ( sind .ne. 0 ) then
              r(j) = r(j) + 0.5d0 * x(sind) ** 2
          end if
      end do

C     Bound constraints

      do i = 1,nbds
          cind = m + i
          sind = slaind(cind)
          vind = consti(i)
          if ( constt(i) .eq. 'L' ) then
              r(cind) = l(vind) - x(vind)
          else
              r(cind) = x(vind) - u(vind)
          end if

          r(cind) = r(cind) + 0.5d0 * x(sind) ** 2
      end do

C     Constraints norm

      cnorm = 0.0d0
      do j = 1,mtot
          cnorm = max( cnorm, abs( r(j) ) )
      end do

      if ( iter .eq. 0 ) then
          cnorm0 = cnorm
      end if

C     ==================================================================
C     COMPUTE FIRST DERIVATIVES
C     ==================================================================

C     Gradient of the objective function and Jacobian of the constraints

      if ( gjaccoded ) then

          call sevalgjac(n,x,t,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

C         First derivative related to the slacks of the regular
C         constraints

          do j = 1,m
              sind = slaind(j)

              if ( sind .ne. 0 ) then
                  jcnnz = jcnnz + 1

                  jcfun(jcnnz) = j
                  jcvar(jcnnz) = sind
                  jcval(jcnnz) = x(sind)
              end if
          end do

          call coo2csr(m,jcnnz,jcfun,jcvar,jcval,jclen,jcsta)

      else

C         Gradient of the objective function

          call sevalg(n,x,t,inform)
          if ( inform .lt. 0 ) return

C         Jacobian of regular constraints with the slacks effect

          k = 0

          do j = 1,m
              jcsta(j) = k + 1

              call sevaljac(n,x,j,jcvar(k+1),jcval(k+1),jclen(j),inform)
              if ( inform .lt. 0 ) return

              k = k + jclen(j)

C             First derivative related to the slacks of the regular
C             constraints

              sind = slaind(j)

              if ( sind .ne. 0 ) then
                  jclen(j) = jclen(j) + 1

                  jcvar(k+1) = sind
                  jcval(k+1) = x(sind)

                  k = k + 1
              end if
          end do

          jcnnz = k

      end if

C     Bound constraints with the slacks effect

      do i = 1,nbds
          cind = m + i
          sind = slaind(cind)
          vind = consti(i)

          jcsta(cind) = jcnnz + 1
          jclen(cind) = 2

          jcvar(jcnnz+1) = vind
          if ( constt(i) .eq. 'L' ) then
              jcval(jcnnz+1) = - 1.0d0
          else
              jcval(jcnnz+1) =   1.0d0
          end if

          jcvar(jcnnz+2) = sind
          jcval(jcnnz+2) = x(sind)

          jcnnz = jcnnz + jclen(cind)
      end do

C     Gradient of the Lagrangian

      do i = 1,ntot
          nl(i) = 0.0d0
      end do

C     Gradient of the objective function

      do i = 1,n
          nl(i) = nl(i) + t(i)
      end do

C     Effect of the Jacobian of regular constraints

      do j = 1,m
          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
          end do
      end do

C     Set Lagrange multipliers related to bound constraints

      do i = 1,nbds
          cind = m + i
          vind = consti(i)
          if ( status(vind) .ne. 'F' ) then
              if ( constt(i) .eq. status(vind) ) then
                  if ( constt(i) .eq. 'L' ) then
                      lambda(cind) =   nl(vind)
                  else
                      lambda(cind) = - nl(vind)
                  end if
              else
                  lambda(cind) = 0.0d0
              end if
          else
              lambda(cind) = 0.0d0
          end if
      end do

C     Effect of Jacobian of bound constraints

      do j = m + 1,m + nbds
          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
          end do
      end do

C     Gradient of the Lagrangian norm

      nlnorm = 0.0d0
      do i = 1,ntot
          nlnorm = max( nlnorm, abs( nl(i) ) )
      end do

      if ( iter .eq. 0 ) then
          nlnorm0 = nlnorm
      end if

C     ==================================================================
C     WRITE INFORMATION OF THE CURRENT POINT
C     ==================================================================

      if ( iprintout .ge. 1 ) then
          write(* ,8010) iter,f,cnorm,nlnorm
          write(10,8010) iter,f,cnorm,nlnorm
      end if

C     ==================================================================
C     TEST STOPPING CRITERIA
C     ==================================================================

      if ( cnorm .le. epsfeas .and. nlnorm .le. epsopt ) then
        ! THE POINT SATISFIES FEASIBILITY AND OPTIMALITY

        ! (Ignored constraints and Lagrange multipliers signal
        !  constraint must be checked)

          go to 400
      end if

      if ( iter .ge. maxit ) then
        ! MAXIMUM NUMBER OF ITERATIONS EXCEEDED

          accinfo = 3

          if ( iprintout .ge. 1 ) then
              write(* ,9030)
              write(10,9030)
          end if

          go to 500
      end if

      cnormref  = max(cnorm0 ,max(cnormprev ,max(epsfeas, 1.0d+01)))
      nlnormref = max(nlnorm0,max(nlnormprev,max(epsopt,  1.0d+01)))

      if ( cnorm .gt. cnormref  .or. nlnorm .gt. nlnormref .or.
     +     cnorm .eq. cnormprev .or. nlnorm .eq. nlnormprev ) then
        ! IT SEEMS TO BE DIVERGING

          divit = divit + 1

          if ( divit .ge. 3 ) then

              accinfo = 4

              if ( iprintout .ge. 1 ) then
                  write(* ,9040)
                  write(10,9040)
              end if

              go to 500

          end if

      else
          divit = 0
      end if

C     ==================================================================
C     DO AN ITERATION
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     COMPUTE SECOND DERIVATIVES
C     ==================================================================

C     Hessian of the Lagrangian

      call sevalhl(n,x,m,lambda,hlin,hcol,hval,hlnnz,inform)
      if ( inform .lt. 0 ) return

C     Second derivatives related to slacks of the regular constraints

      do j = 1,m
          sind = slaind(j)

          if ( sind .ne. 0 ) then
              hlnnz = hlnnz + 1

              hlin(hlnnz) = sind
              hcol(hlnnz) = sind
              hval(hlnnz) = lambda(j)
          end if
      end do

C     Second derivatives related to slacks of the bound constraints

      do i = 1,nbds
          cind = m + i
          sind = slaind(cind)

          hlnnz = hlnnz + 1

          hlin(hlnnz) = sind
          hcol(hlnnz) = sind
          hval(hlnnz) = lambda(cind)
      end do

C     ==================================================================
C     ASSEMBLE THE JACOBIAN OF THE KKT SYSTEM
C     ==================================================================

      hnnz = 0

      do i = 1,nsys
          hdiag(i) = 0
      end do

C     Hessian of the Lagrangian

      do i = 1,hlnnz
          if ( hlin(i) .ge. hcol(i) ) then

              lin = inn(hlin(i))
              col = inn(hcol(i))
              val = hval(i)

              if ( lin .le. ninn .and. col .le. ninn ) then
                  if ( val .ne. 0.0d0 ) then
C                     A(lin,col) = A(lin,col) + val
                      hnnz = hnnz + 1
                      hlin(hnnz) = lin
                      hcol(hnnz) = col
                      hval(hnnz) = val
                      if ( lin .eq. col ) hdiag(lin) = hnnz
                  end if
              end if

          end if
      end do

C     Jacobian of the constraints

      do j = 1,mtot
          do i = jcsta(j),jcsta(j) + jclen(j) - 1

              fun = inp(j)
              var = inn(jcvar(i))
              val = jcval(i)

              if ( var .le. ninn .and. fun .le. minp ) then
                  if ( val .ne. 0.0d0 ) then
C                     A(fun+ninn,var) = A(fun+ninn,var) + val
                      hnnz = hnnz + 1
                      hlin(hnnz) = fun + ninn
                      hcol(hnnz) = var
                      hval(hnnz) = val
                  end if
              end if

          end do
      end do

      do i = 1,nsys
          if ( hdiag(i) .eq. 0 ) then
              hnnz = hnnz + 1
              hlin(hnnz) = i
              hcol(hnnz) = i
              hval(hnnz) = 0.0d0

              hdiag(i) = hnnz
          end if
      end do

C     ==================================================================
C     ANALYSE SPARSITY PATTERN
C     ==================================================================

      call lssana(nsys,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          accinfo = 6
          go to 500

      end if

C     ==================================================================
C     SOLVE THE NEWTONIAN SYSTEM
C     ==================================================================

C     ==================================================================
C     COMPUTE REGULARIZATION
C     ==================================================================

 200  continue

      do i = 1,ninn
          adddiag(i) = epsadd
      end do

      do i = ninn + 1,nsys
          adddiag(i) = - macheps12
      end do

C     ==================================================================
C     FACTORIZE THE JACOBIAN OF THE NEWTONIAN SYSTEM
C     ==================================================================

      call lssfac(nsys,hnnz,hlin,hcol,hval,hdiag,adddiag,pind,pval,
     +nneigv,lssinfo)

      if ( lssinfo .eq. 0 .or. lssinfo .eq. 1 ) then

          if ( nneigv .ne. minp ) then
            ! WRONG INERTIA (SEE NOCEDAL AND WRIGHT)

C             Lemma 16.3 [pg. 447]: Suppose the the Jacobian of the
C             constraints has full rank and that the reduced Hessian
C             Z^T H Z is positive definite. Then the Jacobian of the
C             KKT system has ninn positive eigenvalues, minp negative
C             eigenvalues, and no zero eigenvalues.

C             Note that at this point we know that the matrix has no
C             zero eigenvalues. nneigv gives the number of negative
C             eigenvalues.

              epsadd = max( macheps12, epsadd * 10.0d0 )

              if ( iprintout .ge. 1 ) then
                  itmp = ninn+minp-nneigv
                  write(* ,8090) ninn,minp,itmp,nneigv,epsadd
                  write(10,8090) ninn,minp,itmp,nneigv,epsadd
              end if

              go to 200
          end if

      else if ( lssinfo .eq. 2 ) then
        ! SINGULAR JACOBIAN

          epsadd = max( macheps12, epsadd * 10.0d0 )

          if ( iprintout .ge. 1 ) then
              write(* ,8080) epsadd
              write(10,8080) epsadd
          end if

          if ( epsadd .le. 1.0d+20 ) then
              go to 200
          end if

          accinfo = 5

          if ( iprintout .ge. 1 ) then
              write(* ,9050)
              write(10,9050)
          end if

          go to 500

      else if ( lssinfo .eq. 6 ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          accinfo = 6
          go to 500

      else if ( lssinfo .eq. 7 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          accinfo = 7
          go to 500

      else ! if ( lssinfo .eq. 8 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          accinfo = 8
          go to 500

      end if

C     ==================================================================
C     SOLVE TRIANGULAR SYSTEMS
C     ==================================================================

C     SET RHS

      do i = 1,ntot
          if ( inn(i) .le. ninn ) then
              b(inn(i)) = - nl(i)
          end if
      end do

      do j = 1,mtot
          if ( inp(j) .le. minp ) then
              b(ninn+inp(j)) = - r(j)
          end if
      end do

C     SOLVE THE EQUATIONS

      call lsssol(nsys,b)

C     ==================================================================
C     UPDATE x AND lambda
C     ==================================================================

      do i = 1,ntot
          if ( inn(i) .le. ninn ) then
              x(i) = x(i) + b(inn(i))
          end if
      end do

      do j = 1,mtot
          if ( inp(j) .le. minp ) then
              lambda(j) = lambda(j) + b(ninn+inp(j))
          end if
      end do

C     ==================================================================
C     ITERATE
C     ==================================================================

      go to 100

C     ==================================================================
C     END OF MAIN LOOP
C     ==================================================================

C     ==================================================================
C     CHECK FEASIBILITY CONSIDERING ALL CONSTRAINTS
C     ==================================================================

 400  continue

      cnorm = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              cnorm = max( cnorm, abs( r(i) ) )
          else
              cnorm = max( cnorm, r(i) )
          end if
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,8020) cnorm
          write(10,8020) cnorm
      end if

      if ( cnorm .gt. epsfeas ) then
          accinfo = 1

          if ( iprintout .ge. 1 ) then
              write(* ,9010)
              write(10,9010)
          end if

          go to 500
      end if

C     ==================================================================
C     CHECK LAGRANGE MULTIPLIERS SIGNAL (RELATED TO BOUND CONSTRAINTS)
C     (Lagrange mutipliers of inequality constraints must be tested too.)
C     ==================================================================

      sgnvio = 0.0d0
      do i = m + 1,mtot
          sgnvio = max( sgnvio, - lambda(i) )
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,8030) sgnvio
          write(10,8030) sgnvio
      end if

      if ( sgnvio .eq. 0.0d0 ) then
          accinfo = 0

          if ( iprintout .ge. 1 ) then
              write(* ,9000)
              write(10,9000)
          end if

          go to 500
      end if

C     TRY TO CORRECT MULTIPLIERS SIGNAL

      if ( iprintout .ge. 1 ) then
          write(* ,8040)
          write(10,8040)
      end if

      call minsq(n,t,m,nbds,lambda,constt,consti,status,msqiter,inform)
      if ( inform .lt. 0 ) return

C     Compute signal violation

      sgnvio = 0.0d0
      do i = m + 1,mtot
         sgnvio = max( sgnvio, - lambda(i) )
      end do

C     Compute optimality

      do i = 1,n
          nl(i) = t(i)
      end do

      do i = n + 1,ntot
          nl(i) = 0.0d0
      end do

      do j = 1,mtot
          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
          end do
      end do

      nlnorm = 0.0d0
      do i = 1,ntot
          nlnorm = max( nlnorm, abs( nl(i) ) )
      end do

      if ( iprintout .ge. 1 ) then
          write(* ,8050) sgnvio,nlnorm
          write(10,8050) sgnvio,nlnorm
      end if

      if ( nlnorm .le. epsopt ) then
          accinfo = 0

          if ( iprintout .ge. 1 ) then
              write(* ,9000)
              write(10,9000)
          end if

          go to 500

      else
          accinfo = 2

          if ( iprintout .ge. 1 ) then
              write(* ,9020)
              write(10,9020)
          end if

          go to 500
      end if

C     ==================================================================
C     SET SOLUTION
C     ==================================================================

 500  continue

      do i = 1,n
          xo(i) = x(i)
      end do

      do j = 1,m
          lambdao(j) = lambda(j)
      end do

C     ==================================================================
C     NON-EXECUTABLE STATEMENTS
C     ==================================================================

 8000 format(/,' NEWTON-KKT scheme in action!')
 8010 format(/,' NEWTON-KKT Iteration',42X,' = ',5X,I6,
     +       /,' Objective function value',38X,' = ',1PD11.4,
     +       /,' Maximal violation of constraints',30X,' = ',1PD11.4,
     +       /,' Sup-norm of the gradient of the Lagrangian',20X,' = ',
     +           1PD11.4)
 8020 format(/,' Maximal violation of constraints (considering all',
     +         ' constraints) = ',1PD11.4)
 8030 format(/,' Maximal violation of Lagrange multipliers',
     +         ' non-negativity',6X,' = ',1PD11.4)
 8040 format(/,' GENCAN is being called to find the right',
     +         ' multipliers.')
 8050 format(/,' Maximal violation of Lagrange multipliers',
     +         ' non-negativity',6X,' = ',1PD11.4,
     +       /,' Sup-norm of the gradient of the Lagrangian',20X,' = ',
     +           1PD11.4)

 8080 format(/,' Singular Jacobian.',
     +         ' epsadd was increased to ',1PD11.4)
 8090 format(/,' Wrong Jacobian inertia. ',
     +       /,' Desired POS = ',I7,' NEG = ',I7,
     +       /,' Actual  POS = ',I7,' NEG = ',I7,
     +       /,' epsadd was increased to ',1PD11.4)

 9000 format(/,' Flag of NEWTON-KKT = KKT system solved!')
 9010 format(/,' Flag of NEWTON-KKT = Ignored constraints were',
     +         ' violated.')
 9020 format(/,' Flag of NEWTON-KKT = After correcting the Lagrange',
     +         ' multipliers signal,',/,' optimality was lost.')
 9030 format(/,' Flag of NEWTON-KKT = Maximum of iterations reached.')
 9040 format(/,' Flag of NEWTON-KKT = Newton seems to be diverging.')
 9050 format(/,' Flag of NEWTON-KKT = Singular Jacobian.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsq(n,t,m,nbds,lambda,constt,consti,status,iter,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,iter,m,n,nbds

C     ARRAY ARGUMENTS
      character constt(nbds),status(n)
      integer consti(nbds)
      double precision t(n),lambda(m+nbds)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     The problem to be solved is
C
C         min 0.5 || b + A^T lambda ||_2^2
C
C         subject to
C
C         lambda_i >= 0, i = m+1,...,m+nabds.
C
C     Columns a_i of A^T are multiplied by
C
C         an_i = 1 / max( 1, ||a_i||_infty ).
C
C     So, if we define D = diag(an_1,...,an_{m+nabds}), the problem
C     can be rewrittem as
C
C         min 0.5 || b + A^T D D^{-1} lambda ||_2^2
C
C         subject to
C
C         an_i lambda_i >= 0, i = m+1,...,m+nabds.
C
C     Subtituting an_i lambda_i by lambda_i the problem to be solved
C     becomes
C
C         min 0.5 || b + A^T D lambda ||_2^2
C
C         subject to
C
C         lambda_i >= 0, i = m+1,...,m+nabds.

C     PARAMETERS
      integer mtotmax,ntotmax
      parameter ( mtotmax  = mmax + 2 * nmax )
      parameter ( ntotmax  = mmax + 3 * nmax )

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      logical avoiddstmp,dum2(mtotmax),dum3(mtotmax),useustptmp
      integer cind,geninfo,i,j,k,maxit,nabds,vind
      double precision dum1(mtotmax),dum4(mtotmax),dum5,dum6,eps,msqf,
     +                 msqnlpsupn

C     LOCAL ARRAYS
      integer ind(mtotmax)
      double precision l(mtotmax),msqnl(mtotmax),scol(mtotmax),
     +        u(mtotmax),x(mtotmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

C     RHS

      do i = 1,n
          b(i) = t(i)
      end do

C     Matrix

      hnnz = 0

      do j = 1,m
          ind(j) = j

          do i = jcsta(j),jcsta(j) + jclen(j) - 1
              if ( jcvar(i) .le. n ) then
                  hnnz = hnnz + 1
                  hcol(hnnz) = j
                  hlin(hnnz) = jcvar(i)
                  hval(hnnz) = jcval(i)
              end if
          end do

          x(j) = lambda(j)
          l(j) = - 1.0d+20
          u(j) =   1.0d+20
      end do

      nabds = 0
      do j = 1,nbds
          cind = m + j
          vind = consti(j)

          if ( constt(j) .eq. status(vind) ) then
              nabds = nabds + 1
              k = m + nabds

              ind(k) = cind

              do i = jcsta(cind),jcsta(cind) + jclen(cind) - 1
                  if ( jcvar(i) .le. n ) then
                      hnnz = hnnz + 1
                      hcol(hnnz) = k
                      hlin(hnnz) = jcvar(i)
                      hval(hnnz) = jcval(i)
                  end if
              end do

              x(k) = lambda(cind)
              l(k) = 0.0d0
              u(k) = 1.0d+20
          end if
      end do

C     Dimensions

      nrows = n
      ncols = m + nabds

C     Columns scaling

      do j = 1,ncols
          scol(j) = 1.0d0
      end do

      do i = 1,hnnz
          scol(hcol(i)) = max( scol(hcol(i)), abs( hval(i) ) )
      end do

      do j = 1,ncols
          scol(j) = 1.0d0 / scol(j)
      end do

      do i = 1,hnnz
          hval(i) = hval(i) * scol(hcol(i))
      end do

      do i = 1,ncols
          x(i) = x(i) / scol(i)
      end do

C     Call the solver

      innercall = .true.

      avoiddstmp = avoidds
      useustptmp = useustp

      avoidds = .true.
      useustp = .true.

      eps   = 1.0d-16
      maxit = 200

      call gencan(ncols,x,l,u,0,dum1,dum2,dum3,dum4,0.0d0,eps,
     +maxit,iter,msqf,msqnl,msqnlpsupn,dum5,dum6,geninfo,inform)

      innercall = .false.

      avoidds = avoiddstmp
      useustp = useustptmp

C     Copy unscaled solution to lambda

      do i = 1,ncols
          lambda(ind(i)) = x(i) * scol(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsqf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision p(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          p(i) = b(i)
      end do

      do i = 1,hnnz
          p(hlin(i)) = p(hlin(i)) + x(hcol(i)) * hval(i)
      end do

      f = 0.0d0
      do i = 1,nrows
          f = f + p(i) ** 2
      end do

      f = 0.5d0 * f

      f = 1.0d+08 * f

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsqg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision p(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          p(i) = b(i)
      end do

      do i = 1,hnnz
          p(hlin(i)) = p(hlin(i)) + x(hcol(i)) * hval(i)
      end do

      do i = 1,ncols
          g(i) = 0.0d0
      end do

      do i = 1,hnnz
          g(hcol(i)) = g(hcol(i)) + p(hlin(i)) * hval(i)
      end do

      do i = 1,ncols
          g(i) = 1.0d+08 * g(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine minsqhp(n,x,p,hp,goth,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer inform,n

C     ARRAY ARGUMENTS
      double precision hp(n),p(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision tmp(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          tmp(i) = 0.0d0
      end do

      do i = 1,hnnz
          tmp(hlin(i)) = tmp(hlin(i)) + p(hcol(i)) * hval(i)
      end do

      do i = 1,ncols
          hp(i) = 0.0d0
      end do

      do i = 1,hnnz
          hp(hcol(i)) = hp(hcol(i)) + tmp(hlin(i)) * hval(i)
      end do

      do i = 1,ncols
          hp(i) = 1.0d+08 * hp(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      logical function minsqstop(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     COMMON SCALARS
      integer ncols,nrows

C     COMMON ARRAYS
      double precision b(nmax)

C     LOCAL SCALARS
      integer i
      double precision pnorm

C     LOCAL ARRAYS
      double precision p(nmax)

C     COMMON BLOCKS
      common /prodat/ b,ncols,nrows

      do i = 1,nrows
          p(i) = b(i)
      end do

      do i = 1,hnnz
          p(hlin(i)) = p(hlin(i)) + x(hcol(i)) * hval(i)
      end do

      pnorm = 0.0d0
      do i = 1,nrows
          pnorm = max( pnorm, abs( p(i) ) )
      end do

      minsqstop = .false.
      if ( pnorm .le. macheps12 ) then
          minsqstop = .true.
      end if

      end
C     ******************************************************************
C     ******************************************************************
      subroutine fparam(epsfeas,epsopt,iprint,ncomp)

      implicit none

C     SCALAR ARGUMENTS
      integer ncomp,iprint
      double precision epsfeas,epsopt

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     PARAMETERS
      integer nwords
      parameter ( nwords = 14 )

C     LOCAL SCALARS
      logical lss,scl
      integer i,ifirst,ikey,ilast,inum,j
      double precision dnum

C     LOCAL ARRAYS
      character * 80 line
      character * 10 keyword
      character *  4 lsssub
      character *  4 sclsub

C     DATA BLOCKS
      character * 1  addinfo(nwords),lower(26),upper(26)
      character * 10 dictionary(nwords)
      character * 38 description(nwords)

      data lower /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +            'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data upper /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      data dictionary( 1) /'INCREMENTA'/
      data dictionary( 2) /'HESSIAN-AP'/
      data dictionary( 3) /'TRUE-HESSI'/
      data dictionary( 4) /'PERFORM-AC'/
      data dictionary( 5) /'DIRECT-SOL'/
      data dictionary( 6) /'SCALE-LINE'/
      data dictionary( 7) /'DO-NOT-REM'/
      data dictionary( 8) /'ADD-SLACKS'/
      data dictionary( 9) /'DO-NOT-SCA'/
      data dictionary(10) /'IGNORE-OBJ'/
      data dictionary(11) /'FEASIBILIT'/
      data dictionary(12) /'OPTIMALITY'/
      data dictionary(13) /'OUTPUT-DET'/
      data dictionary(14) /'NCOMP-ARRA'/

      data description( 1) /'INCREMENTAL-QUOTIENTS-IN-CG           '/
      data description( 2) /'HESSIAN-APPROXIMATION-IN-CG           '/
      data description( 3) /'TRUE-HESSIAN-PRODUCT-IN-CG            '/
      data description( 4) /'PERFORM-ACCELERATION-STEP             '/
      data description( 5) /'DIRECT-SOLVER                         '/
      data description( 6) /'SCALE-LINEAR-SYSTEMS                  '/
      data description( 7) /'DO-NOT-REMOVE-FIXED-VARIABLES         '/
      data description( 8) /'ADD-SLACKS                            '/
      data description( 9) /'DO-NOT-SCALE-OBJECTIVE-AND-CONSTRAINTS'/
      data description(10) /'IGNORE-OBJECTIVE-FUNCTION             '/
      data description(11) /'FEASIBILITY-TOLERANCE                 '/
      data description(12) /'OPTIMALITY-TOLERANCE                  '/
      data description(13) /'OUTPUT-DETAIL                         '/
      data description(14) /'NCOMP-ARRAY                           '/

      data addinfo( 1) /' '/
      data addinfo( 2) /' '/
      data addinfo( 3) /' '/
      data addinfo( 4) /' '/
      data addinfo( 5) /' '/
      data addinfo( 6) /' '/
      data addinfo( 7) /' '/
      data addinfo( 8) /' '/
      data addinfo( 9) /' '/
      data addinfo(10) /' '/
      data addinfo(11) /'D'/
      data addinfo(12) /'D'/
      data addinfo(13) /'I'/
      data addinfo(14) /'I'/

C     EXTERNAL FUNCTIONS
      external lss,scl

C     BANNER
      if ( iprintctl(1) ) then
          write(*, 8000)
          write(10,8000)
      end if

C     OPENING THE SPECIFICATION FILE
      open(20,err=300,file='algencan.dat',status='old')

      if ( iprintctl(2) ) then
          write(*, 9005)
          write(10,9005)
      end if

C     MAIN LOOP

 100  continue

C     READING LINES
      read(20,fmt=1000,err=400,end=200) line

C     PROCESS LINES

C     Find first character
      i = 1
 110  if ( i .le. 80 .and. line(i:i) .eq. ' ' ) then
          i = i + 1
          go to 110
      end if

C     Skip blank lines
      if ( i .gt. 80 ) then
          go to 100
      end if

      ifirst = i

C     Skip comments
      if ( line(ifirst:ifirst) .eq. '*' .or.
     +     line(ifirst:ifirst) .eq. '#' ) then
          go to 100
      end if

C     Find the end of the keyword
      i = ifirst + 1
 120  if ( i .le. 80 .and. line(i:i) .ne. ' ' ) then
          i = i + 1
          go to 120
      end if

      ilast = i - 1

C     Obtain the first 10 characters and convert to upper-case
      keyword = '          '
      do i = 1,min( 10, ilast - ifirst + 1 )
          keyword(i:i) = line(ifirst+i-1:ifirst+i-1)
          do j = 1,26
              if ( keyword(i:i) .eq. lower(j) ) then
                  keyword(i:i) = upper(j)
              end if
          end do
      end do

C     Look up the keyword in the dictionary
      i = 1
 130  if ( i .le. nwords .and. keyword .ne. dictionary(i) ) then
          i = i + 1
          go to 130
      end if

C     Ignore unknown keywords
      if ( i .gt. nwords ) then
          if ( iprintctl(2) ) then
              write(*, 9020) line(ifirst:ilast)
              write(10,9020) line(ifirst:ilast)
          end if
          go to 100
      end if

      ikey = i

C     Read additional information if needed
      if ( addinfo(ikey) .ne. ' ' ) then

C         Skip blanks
          i = ilast + 1
 140      if ( i .le. 80 .and. line(i:i) .eq. ' ' ) then
              i = i + 1
              go to 140
          end if

C         Ignore keywords without the required information
          if ( i .gt. 80 ) then
              if ( iprintctl(2) ) then
                  write(*, 9030) description(ikey)
                  write(10,9030) description(ikey)
              end if
              go to 100
          end if

C         Read additional information
          if ( addinfo(ikey) .eq. 'I' ) then
              read(unit=line(i:80),fmt=2000) inum

          else if ( addinfo(ikey) .eq. 'D' ) then
              read(unit=line(i:80),fmt=3000) dnum
          end if

      end if

C     Process keyword
      if ( iprintctl(2) ) then
          if ( addinfo(ikey) .eq. ' ' ) then
              write(*, 9040) description(ikey)
              write(10,9040) description(ikey)
          else if ( addinfo(ikey) .eq. 'I' ) then
              write(*, 9041) description(ikey),inum
              write(10,9041) description(ikey),inum
          else if ( addinfo(ikey) .eq. 'D' ) then
              write(*, 9042) description(ikey),dnum
              write(10,9042) description(ikey),dnum
          end if
      end if

C     Set the corresponding algencan argument
      if ( ikey .eq.  1 ) then
          hptype = 'INCQUO'

      else if ( ikey .eq.  2 ) then
          hptype = 'HAPPRO'

      else if ( ikey .eq.  3 ) then
          if ( hlpcoded .or. truehl ) then
              hptype = 'TRUEHL'
          else
              if ( iprintctl(2) ) then
                  write(* ,9100) description(ikey)
                  write(10,9100) description(ikey)
              end if
          end if

      else if ( ikey .eq.  4 ) then
          if ( .not. truehl ) then
              if ( iprintctl(2) ) then
                  write(* ,9110) description(ikey)
                  write(10,9110) description(ikey)
              end if
          else if ( .not. lss(lsssub) ) then
              if ( iprintctl(2) ) then
                  write(* ,9120) description(ikey)
                  write(10,9120) description(ikey)
              end if
          else
              skipacc = .false.
              if ( iprintctl(2) ) then
                  write(* ,9060) lsssub
                  write(10,9060) lsssub
              end if
          end if

      else if ( ikey .eq.  5 ) then
          if ( .not. truehl ) then
              if ( iprintctl(2) ) then
                  write(* ,9110) description(ikey)
                  write(10,9110) description(ikey)
              end if
          else if ( .not. lss(lsssub) ) then
              if ( iprintctl(2) ) then
                  write(* ,9120) description(ikey)
                  write(10,9120) description(ikey)
              end if
          else
              avoidds = .false.
              if ( iprintctl(2) ) then
                  write(* ,9050) lsssub
                  write(10,9050) lsssub
              end if
          end if

      else if ( ikey .eq.  6 ) then
          if ( lss(lsssub) ) then
              if ( lsssub .eq. 'MA57' ) then
                  sclsys = .true.
                  if ( iprintctl(2) ) then
                      write(* ,9070) lsssub
                      write(10,9070) lsssub
                  end if
              else if ( scl(sclsub) ) then
                  sclsys = .true.
                  if ( iprintctl(2) ) then
                      write(* ,9080) sclsub
                      write(10,9080) sclsub
                  end if
              else
                  if ( iprintctl(2) ) then
                      write(* ,9130) description(ikey)
                      write(10,9130) description(ikey)
                  end if
              end if
          else
              if ( iprintctl(2) ) then
                  write(* ,9130) description(ikey)
                  write(10,9130) description(ikey)
              end if
          end if

      else if ( ikey .eq.  7 ) then
          rmfixv = .false.

      else if ( ikey .eq.  8 ) then
          slacks = .true.

      else if ( ikey .eq.  9 ) then
          scale = .false.

      else if ( ikey .eq. 10 ) then
          ignoref = .true.

      else if ( ikey .eq. 11 ) then
          epsfeas = dnum

      else if ( ikey .eq. 12 ) then
          epsopt = dnum

      else if ( ikey .eq. 13 ) then
          iprint = inum

      else if ( ikey .eq. 14 ) then
          ncomp = inum
      end if

C     IIERATE
      go to 100

C     END OF LOOP

C     TERMINATIONS

C     CLOSING SPECIFICATION FILE
 200  continue
      close(20)
      go to 500

C     NO SPECIFICATION FILE
 300  continue
      if ( iprintctl(2) ) then
          write(*, 9000)
          write(10,9000)
      end if
      go to 500

C     ERROR READING THE SPECIFICATION FILE
 400  continue
      if ( iprintctl(2) ) then
          write(*, 9010)
          write(10,9010)
      end if
      go to 500

C     PRINTING PARAMETERS VALUES
 500  continue
      if ( iprintctl(2) ) then
          write(* ,4000) hptype,avoidds,skipacc,sclsys,rmfixv,slacks,
     +                   scale,epsfeas,epsopt,iprint,ncomp
          write(10,4000) hptype,avoidds,skipacc,sclsys,rmfixv,slacks,
     +                   scale,epsfeas,epsopt,iprint,ncomp
      end if

C     NON-EXECUTABLE STATEMENTS

 1000 format(A80)
 2000 format(BN,I20)
 3000 format(BN,F24.0)
 4000 format(/,1X,'ALGENCAN PARAMETERS:',
     +       /,1X,'hptype   = ',     14X,A6,
     +       /,1X,'avoidds  = ',     19X,L1,
     +       /,1X,'skipacc  = ',     19X,L1,
     +       /,1X,'sclsys   = ',     19X,L1,
     +       /,1X,'rmfixv   = ',     19X,L1,
     +       /,1X,'slacks   = ',     19X,L1,
     +       /,1X,'scale    = ',     19X,L1,
     +       /,1X,'epsfeas  = ',8X,1P,D12.4,
     +       /,1X,'epsopt   = ',8X,1P,D12.4,
     +       /,1X,'iprint   = ',        I20,
     +       /,1X,'ncomp    = ',        I20)
 8000 format(/,1X,78('='),
     +       /,1X,'This is ALGENCAN 2.2.1.',
     +       /,1X,'ALGENCAN, an augmented Lagrangian method for ',
     +            'nonlinear programming, is part of',/,1X,'the TANGO ',
     +            'Project: Trustable Algorithms for Nonlinear ',
     +            'General Optimization.',/,1X,'See ',
     +            'http://www.ime.usp.br/~egbirgin/tango/ for details.',
     +       /,1X,78('='))
 9000 format(/,1X,'The optional specification file algencan.dat was ',
     +            'not found in the current',/,1X,'directory (this is ',
     +            'not a problem nor an error). The default values ',
     +            'for the',/,1X,'ALGENCAN parameters will be used.')
 9005 format(/,1X,'Specification file algencan.dat is being used.')
 9010 format(/,1X,'Error reading specification file algencan.dat.')
 9020 format(  1X,'Ignoring unknown keyword ',A38)
 9030 format(  1X,'Ignoring incomplete keyword ',A38)
 9040 format(1X,A38)
 9041 format(1X,A38,5X,I20)
 9042 format(1X,A38,1X,1P,D24.8)

 9050 format(1X,'(Subroutine ',A4,' from HSL will be used as a direct ',
     +          'solver for Newtonian',/,1X,'linear systems.)')
 9060 format(1X,'(Subroutine ',A4,' from HSL will be used as a direct ',
     +          'solver for the',/,1X,'acceleration-step linear ',
     +          'systems.)')
 9070 format(1X,'(Linear systems will be scaled using the embedded ',
     +          'scaling option of',/,1X,'subroutine ',A4,' from HSL.)')
 9080 format(1X,'(Subroutine ',A4,' from HSL will be used for scaling ',
     +          'linear systems.)')

 9100 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutines EVALH and EVALHC, ',
     +            'or, alternatively,',/,1X,'subroutine EVALHLP, to ',
     +            'be coded by the user. If you already coded them,',
     +       /,1X,'set array CODED in subrutine INIP appropiately.',/)
 9110 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutines EVALH and EVALHC, ',
     +            'or, alternatively,',/,1X,'subroutine EVALHL, to ',
     +            'be coded by the user. If you already coded them,',
     +       /,1X,'set array CODED in subrutine INIP appropiately.',/)
 9120 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutine MA27 or MA57 from ',
     +            'HSL to be provided by the',/,1X,'user. If you ',
     +            'have any of them, see the compilation instructions ',
     +            'for details.',/)
 9130 format(/,1X,'Warning: Ignoring keyword ',A38,'.',
     +       /,1X,'This option requires subroutine MA57 (that has an ',
     +            'embedded scaling option)',/,1X,'or, to be used in ',
     +            'connection with subroutine MA27, subroutine MC30 ',
     +            'or MC77',/,1X,'from HSL to be provided by the ',
     +            'user. If you have any of them, see the',/,1X,
     +            'compilation instructions for details.',/)

      end
C     *****************************************************************
C     *****************************************************************

      subroutine checkd(n,l,u,m,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision l(n),u(n)

C     This subrotutine checks the user supplied first and second
C     derivatives subroutines (evalg, evalh, evaljac and evalhc) for
C     computing the objective function gradient and Hessian and the
C     constraints gradients and Hessians, respectively.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

C     LOCAL SCALARS
      character answer
      integer i,j
      double precision drand,seed,smalll,smallu

C     LOCAL ARRAYS
      double precision x(nmax)

C     EXTERNAL FUNCTIONS
      external drand

      call vunsetp()

C     SET A RANDOM POINT 

      seed = 123456.0d0
      do i = 1,n
          smalll = max( l(i), - 10.0d0 )
          smallu = min( u(i),   10.0d0 )
          if ( .not. smalll .lt. smallu ) then
              smalll = l(i)
              smallu = u(i)
          end if
          x(i) = smalll + ( smallu - smalll ) * drand(seed)
      end do
 
      write(* ,100)
      write(10,100)

      do i = 1,n
          write(* ,110) i,x(i)
          write(10,110) i,x(i)
      end do

C     CHECK OBJECTIVE FUNCTION GRADIENT

      if ( .not. gcoded ) then
          write(* ,160) 'evalg'
          write(10,160) 'evalg'

          go to 1000
      end if

      write(* ,120)
      write(10,120)

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkg(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK CONSTRAINTS GRADIENTS

 1000 continue

      if ( .not. jaccoded ) then
          write(* ,160) 'evaljac'
          write(10,160) 'evaljac'

          go to 1020
      end if

      j = 1

 1010 if ( j .le. m ) then

          write(* ,130) j
          write(10,130) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              go to 1020

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkjac(n,x,j,inform)
              if ( inform .lt. 0 ) return
          end if

          j = j + 1

          go to 1010

      end if

C     CHECK HESSIAN OF THE OBJECTIVE FUNCTION

 1020 continue

      if ( .not. hcoded ) then
          write(* ,160) 'evalh'
          write(10,160) 'evalh'

          go to 1030
      end if

      write(* ,140)
      write(10,140)

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          return

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkh(n,x,inform)
          if ( inform .lt. 0 ) return
      end if

C     CHECK HESSIANS OF THE CONSTRAINTS

 1030  continue

      if ( .not. hccoded ) then
          write(* ,160) 'evalhc'
          write(10,160) 'evalhc'

          go to 1050
      end if

      j = 1

 1040 if ( j .le. m ) then

          write(* ,150) j
          write(10,150) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              return

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              return

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkhc(n,x,j,inform)
              if ( inform .lt. 0 ) return
          end if

          j = j + 1

          go to 1040

      end if

C     CHECK HESSIAN OF THE LAGRANGIAN

 1050 continue

      if ( .not. hlcoded ) then
          write(* ,160) 'evalhl'
          write(10,160) 'evalhl'

          go to 1060
      end if

      write(* ,*) 'Test of evalhl not implemented yet!'
      write(10,*) 'Test of evalhl not implemented yet!'

C     CHECK HESSIAN OF THE LAGRANGIAN TIMES A VECTOR

 1060 continue

      if ( .not. hlpcoded ) then
          write(* ,160) 'evalhlp'
          write(10,160) 'evalhlp'

          return
      end if

      write(* ,*) 'Test of evalhlp not implemented yet!'
      write(10,*) 'Test of evalhlp not implemented yet!'

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Derivatives will be tested at the random point: ')
 110  format(  1X,'x(',I6,') = ',1P,D15.8)
 120  format(/,1X,'Check the gradient of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 130  format(/,1X,'Check the gradient of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 140  format(/,1X,'Check the Hessian matrix of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 150  format(/,1X,'Check the Hessian of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 160  format(/,1X,'Skipping test of uncoded ',A7,' subroutine.')

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkg(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalg for 
C     computing the gradient of the objective function using central
C     finite differences with two different discretization steps.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer i
      double precision fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      double precision g(nmax)

      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

      write(* ,100)
      write(10,100)

      maxerr = 0.0d0

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )

          write(* ,110) i,g(i),gdiff1,gdiff2,tmp
          write(10,110) i,g(i),gdiff1,gdiff2,tmp

          maxerr = max( maxerr, tmp )

      end do

      write(* ,120) maxerr
      write(10,120) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of the objective function.',
     +       /,1X,'Index',13X,'evalg',2X,'Central diff (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkh(n,x,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalh for 
C     computing the Hessian of the objective function using central 
C     finite differences with two different discretization steps.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
          write(*, 100) nsmax,nsmax
          write(10,100) nsmax,nsmax

          return
      end if

C     Compute the gradient of the objective function at x

      call vevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

C     Compute the Hessian of the objective function at x and save in a
C     dense matrix

      call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

C     Test column by column

      write(* ,200)
      write(10,200)

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vevalg(n,x,gplus1,inform)
          if ( inform .lt. 0 ) return

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vevalg(n,x,gplus2,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp

          write(* ,210) j
          write(10,210) j

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
              if ( elem   .ne. 0.0d0 .or. 
     +             hdiff1 .ne. 0.0d0 .or. 
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      write(* ,220)
                      write(10,220)
                  end if
                  write(* ,230) i,elem,hdiff1,hdiff2,tmp
                  write(10,230) i,elem,hdiff1,hdiff2,tmp
              end if
              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( nullcol ) then
              write(* ,240)
              write(10,240)
          else
              write(* ,250) maxcoe(j)
              write(10,250) maxcoe(j)
          end if

      end do

      write(* ,*)
      write(10,*)

      do j = 1,n
          write(* ,260) j,maxcoe(j)
          write(10,260) j,maxcoe(j)
      end do

      write(* ,270) maxerr
      write(10,270) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKH uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of the objective function column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkjac(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evaljac for 
C     computing the gradients of the constraints using central finite 
C     differences with two different discretization steps.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      logical nullcol
      integer i,jcnnz
      double precision cminus,cplus,jacdiff1,jacdiff2,maxerr,step1,
     +        step2,tmp

C     LOCAL ARRAYS
      integer jcvar(nmax)
      double precision g(nmax),jcval(nmax)

C     COMPUTE THE GRADIENT OF THE CONSTRAINT AND SAVE IT INTO A DENSE 
C     VECTOR

      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     COMPARE WITH CENTRAL FINITE DIFFERENCES

      write(* ,100) ind
      write(10,100) ind

      maxerr = 0.0d0

      nullcol = .true.

      do i = 1,n
          tmp  = x(i)

          step1 = macheps13 * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step1
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jacdiff1 = ( cplus - cminus ) / ( 2.0d0 * step1 )

          step2 = macheps13 * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step2
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp

          jacdiff2 = ( cplus - cminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - jacdiff1 ), abs( g(i) - jacdiff2 ) )

          if ( g(i)     .ne. 0.0d0 .or. 
     +         jacdiff1 .ne. 0.0d0 .or. 
     +         jacdiff2 .ne. 0.0d0 ) then
              if ( nullcol ) then
                  nullcol = .false.
                  write(* ,110)
                  write(10,110)
              end if
              write(* ,120) i,g(i),jacdiff1,jacdiff2,tmp
              write(10,120) i,g(i),jacdiff1,jacdiff2,tmp
          end if

          maxerr = max( maxerr, tmp )
      end do

      if ( nullcol ) then
          write(* ,130)
          write(10,130)
      else
          write(* ,140) maxerr
          write(10,140) maxerr
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Gradient vector of constraints ',I5,'.')
 110  format(/,1X,'Index',11X,'evaljac',2X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 120  format(  1X,I5,4(3X,1P,D15.8))
 130  format(  1X,'All the elements of this gradient are null.')
 140  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************
      subroutine checkhc(n,x,ind,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n

C     ARRAY ARGUMENTS
      double precision x(n)

C     This subrotutine checks the user supplied subroutine evalhc for 
C     computing the Hessians of the constraints using finite
C     differences.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      logical nullcol
      integer i,j,hnnz,jcnnz
      double precision elem,hdiff1,hdiff2,maxerr,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2),jcvar(nsmax)
      double precision g(nsmax),gplus1(nsmax),gplus2(nsmax),
     +        H(nsmax,nsmax),hval(nsmax**2),jcval(nsmax),maxcoe(nsmax)

C     Check viability of the test

      if ( n .gt. nsmax ) then
          write(*, 100) nsmax,nsmax
          write(10,100) nsmax,nsmax

          return
      end if

C     Compute the gradient of constraint ind at x and save it in a 
C     dense vector

      call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,jcnnz
          g(jcvar(i)) = g(jcvar(i)) + jcval(i)
      end do

C     Compute the Hessian of constraint ind at x and save it in a 
C     dense matrix

      call vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,hnnz
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

      write(* ,200) ind
      write(10,200) ind

      maxerr = 0.0d0

      do j = 1,n

          tmp  = x(j)

C         Compute the gradient of constraint ind at xplus1 and 
C         save in a dense vector

          step1 = macheps12 * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus1(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus1(jcvar(i)) = jcval(i)
          end do

C         Compute the gradient of constraint ind at xplus2 and 
C         save in a dense vector

          step2 = macheps12 * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              gplus2(i) = 0.0d0
          end do

          do i = 1,jcnnz
              gplus2(jcvar(i)) = jcval(i)
          end do

          x(j) = tmp

          write(* ,210) j
          write(10,210) j

          maxcoe(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )

              if ( elem   .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or. 
     +             hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      write(* ,220)
                      write(10,220)
                  end if
                  write(* ,230) i,elem,hdiff1,hdiff2,tmp
                  write(10,230) i,elem,hdiff1,hdiff2,tmp
              end if

              maxcoe(j) = max( maxcoe(j), tmp )
          end do

          maxerr = max( maxerr, maxcoe(j) )

          if ( nullcol ) then
              write(* ,240)
              write(10,240)
          else
              write(* ,250) maxcoe(j)
              write(10,250) maxcoe(j)
          end if

      end do

      write(* ,*)
      write(10,*)

      do j = 1,n
          write(* ,260) j,maxcoe(j)
          write(10,260) j,maxcoe(j)
      end do

      write(* ,270) maxerr
      write(10,270) maxerr

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Subroutine CHECKHC uses dense matrices up to ',
     +            'dimension ',I6,' times ',I6,'. The Hessian ',
     +            'checking will be skipped.')

 200  format(/,1X,'Hessian matrix of constraint ',I5,' column by ',
     +            'column.')
 210  format(/,1X,'Column:  ',I6)
 220  format(/,1X,'Index',12X,'evalhc',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 230  format(  1X,I5,4(3X,1P,D15.8))
 240  format(  1X,'All the elements of this column are null.')
 250  format(  1X,'Maximum absolute error = ',1P,D15.8)
 260  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)
 270  format(/,1X,'Overall maximum absolute error = ',1P,D15.8)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine setp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine unsetp()

      implicit none

      end
C     ******************************************************************
C     ******************************************************************

      subroutine calcal(nind,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,nind
      double precision al

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Compute augmented Lagrangian

      call sevalal(nt,x,m,lambda,rho,equatn,linear,al,inform)
      if ( inform .lt. 0 ) return

C     Shrink x to the reduced space

      call shrink(nind,x)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcnal(nind,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,nind

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(*),rho(m),x(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Compute the gradient of the augmented Lagrangian

      call sevalnal(nt,x,m,lambda,rho,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

C     Shrink x and nal to the reduced space

      call shrink(nind,x)
      call shrink(nind,nal)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calchal(nind,x,m,lambda,rho,equatn,linear,hlin,hcol,
     +hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,nind

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer hcol(*),hlin(*)
      double precision lambda(m),hval(*),rho(m),x(*)

C     This subroutine computes the Hessian of the augmented Lagrangian
C     in the reduced space.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer col,i,k,lin

C     LOCAL ARRAYS
      integer wi(nmax)

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Compute the Hessian of the augmented Lagrangian

      call sevalhal(nt,x,m,lambda,rho,equatn,linear,hlin,hcol,hval,hnnz,
     +inform)
      if ( inform .lt. 0 ) return

C     Shrink x to the reduced space

      call shrink(nind,x)

C     Shrink representation of H
      
      do i = 1,nt
         wi(i) = 0
      end do

      do i = 1,nind
         wi(ind(i)) = i
      end do
      
      k = 0

      do i = 1,hnnz
          lin = wi(hlin(i))
          col = wi(hcol(i))

          if ( lin .ne. 0 .and. col .ne. 0 ) then
              k = k + 1
              hlin(k) = lin
              hcol(k) = col
              hval(k) = hval(i)
          end if
      end do

      hnnz = k

      end


C     ******************************************************************
C     ******************************************************************

      subroutine calchalp(nind,x,m,lambda,rho,equatn,linear,p,hp,gothl,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,nind

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(*),lambda(m),p(*),rho(m),x(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Complete p with zeroes

      do i = 1,nt - nind
          p(nind+i) = 0.0d0
      end do

C     Expand x and p to the full space

      call expand(nind,x)
      call expand(nind,p)

C     Compute the Hessian-vector product

      call sevalhalp(nt,x,m,lambda,rho,equatn,linear,p,hp,gothl,inform)
      if ( inform .lt. 0 ) return

C     Shrink x, p and hp to the reduced space

      call shrink(nind,x)
      call shrink(nind,p)
      call shrink(nind,hp)
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine capplyhpre(nind,m,rho,equatn,gotp,r,z)

      implicit none

C     SCALAR ARGUMENTS
      logical gotp
      integer m,nind

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision r(*),rho(m),z(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i

C     Complete r with zeroes

      do i = nind + 1,nt
          r(i) = 0.0d0
      end do

C     Expand r to the full space

      call expand(nind,r)

C     Solve P z = r

      call applyhpre(nt,m,rho,equatn,gotp,r,z)

C     Shrink r and z to the reduced space

      call shrink(nind,r)
      call shrink(nind,z)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine csetp(nind,x)

      implicit none

C     SCALAR ARGUMENTS
      integer nind

C     ARRAY ARGUMENTS
      double precision x(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = 1,nt - nind
          x(nind+i) = xcomplement(i)
      end do

C     Expand x to the full space

      call expand(nind,x)

C     Set point

      call ssetp(nt,x)

C     Shrink x

      call shrink(nind,x)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine shrink(nind,v)

      implicit none

C     This subroutine shrinks vector v from the full dimension space 
C     (dimension n) to the reduced space (dimension nind).

C     SCALAR ARGUMENTS
      integer nind

C     ARRAY ARGUMENTS
      double precision v(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = 1,nind
           indi = ind(i)
           if ( i .ne. indi ) then
               tmp     = v(indi)
               v(indi) = v(i)
               v(i)    = tmp
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine expand(nind,v)

      implicit none

C     This subroutine expands vector v from the reduced space 
C     (dimension nind) to the full space (dimension n).

C     SCALAR ARGUMENTS
      integer nind

C     ARRAY ARGUMENTS
      double precision v(*)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer nt

C     COMMON ARRAYS
      integer ind(nmax)
      double precision xcomplement(nmax)

C     COMMON BLOCKS
      common /rspace/ xcomplement,ind,nt
      save   /rspace/

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = nind,1,- 1
          indi = ind(i)
          if ( i .ne. indi ) then
              tmp     = v(indi)
              v(indi) = v(i)
              v(i)    = tmp
          end if
      end do
     
      end
 
C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

      call ssetp(n,x)

      if ( fccoded ) then
          call sevalobjcb(n,x,f,m,c,inform)

      else ! if ( fcoded .and. ( ccoded .or. m .eq. 0 ) ) then
          call sevalobjca(n,x,f,m,c,inform)
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjcb(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     COMPUTE OBJECTIVE FUNTION AND CONSTRAINTS

      call sevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalobjca(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     LOCAL SCALARS
      integer j

C     COMPUTE OBJECTIVE FUNTION

      call sevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

C     COMPUTE CONSTRAINTS

      do j = 1,m

C         COMPUTE THE j-TH CONSTRAINT
          call sevalc(n,x,j,c(j),inform)
          if ( inform .lt. 0 ) return

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalal(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

      if ( innercall ) then
          call minsqf(n,x,al,inform)
          return
      end if

      call ssetp(n,x)

      if ( fccoded ) then
          call sevalalb(n,x,m,lambda,rho,equatn,linear,al,inform)
      else
          call sevalala(n,x,m,lambda,rho,equatn,linear,al,inform)
      end if

      gotc = .true.

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalalb(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer j
      double precision f,p

C     COMPUTE OBJECTIVE FUNTION AND CONSTRAINTS

      call sevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

C     COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)

      al = f

      do j = 1,m

C         ADD P(c_j, rho_j, lambda_j)
          call evalp(c(j),rho(j),lambda(j),equatn(j),p)
          al = al + p

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalala(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer j
      double precision f,p

C     COMPUTE OBJECTIVE FUNTION

      call sevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

C     COMPUTES AL = f + sum_j P(c_j, rho_j, lambda_j)

      al = f

      do j = 1,m

C         COMPUTE j-TH CONSTRAINT
          call sevalc(n,x,j,c(j),inform)
          if ( inform .lt. 0 ) return

C         ADD P(c_j, rho_j, lambda_j)
          call evalp(c(j),rho(j),lambda(j),equatn(j),p)
          al = al + p

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnl(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

      if ( gjaccoded ) then
          call sevalnlb(n,x,m,lambda,equatn,linear,nl,inform)
      else
          call sevalnla(n,x,m,lambda,equatn,linear,nl,inform)
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnlb(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer i,j,jcnnz

C     LOCAL ARRAYS
      integer jcfun(jcnnzmax)

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND JACOBIAN
C     OF CONSTRAINTS

      call sevalgjac(n,x,nl,m,jcfun,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = nl(i)
      end do

C     CONVERT JACOBIAN OF CONSTRAINTS FROM COORDINATE FORMAT TO
C     COMPRESSED SPARSE ROW FORMAT

      call coo2csr(m,jcnnz,jcfun,jcvar,jcval,jclen,jcsta)

C     COMPUTE \nabla L = \nabla f + \sum_j lambda_j * \nabla c_j

      constrc = .false.

      do j = 1,m

          if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then

C             ADD lambda_j * \nabla c_j

              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
              end do

              if ( .not. linear(j) ) then
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      g(jcvar(i)) = g(jcvar(i)) + lambda(j) * jcval(i)
                  end do
              end if

              constrc = .true.

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnla(n,x,m,lambda,equatn,linear,nl,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer i,ind,j

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION

      call sevalg(n,x,nl,inform)
      if ( inform .lt. 0 ) return

      do i = 1,n
          g(i) = nl(i)
      end do

C     COMPUTE \nabla L = \nabla f + \sum_j lambda_j * \nabla c_j

      ind = 0

      constrc = .false.

      do j = 1,m

          if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then

              jcsta(j) = ind + 1

C             COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
              call sevaljac(n,x,j,jcvar(ind+1),jcval(ind+1),jclen(j),
     +        inform)
              if ( inform .lt. 0 ) return

              ind = ind + jclen(j)

C             ADD lambda_j * \nabla c_j
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
              end do

              if ( .not. linear(j) ) then
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      g(jcvar(i)) = g(jcvar(i)) + lambda(j) * jcval(i)
                  end do
              end if

              constrc = .true.

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnal(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

      if ( innercall ) then
          call minsqg(n,x,nal,inform)
          return
      end if

      if ( gjaccoded ) then
          call sevalnalb(n,x,m,lambda,rho,equatn,linear,nal,inform)
      else
          call sevalnala(n,x,m,lambda,rho,equatn,linear,nal,inform)
      end if

      gotc = .true.

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnalb(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer j
      double precision dum

C     COMPUTE CONSTRAINTS

      if ( .not. gotc ) then
          call sevalfc(n,x,dum,m,c,inform)
          if ( inform .lt. 0 ) return
      end if

      do j = 1,m
C         COMPUTE dP/dc
          call evaldpdy(c(j),rho(j),lambda(j),equatn(j),dpdc(j))
      end do

C     COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
      call sevalnlb(n,x,m,dpdc,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalnala(n,x,m,lambda,rho,equatn,linear,nal,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer j

      if ( .not. gotc ) then
          do j = 1,m
C             COMPUTE THE j-TH CONSTRAINT
              call sevalc(n,x,j,c(j),inform)
              if ( inform .lt. 0 ) return
          end do
      end if

      do j = 1,m
C         COMPUTE dP/dc
          call evaldpdy(c(j),rho(j),lambda(j),equatn(j),dpdc(j))
      end do

C     COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
      call sevalnla(n,x,m,dpdc,equatn,linear,nal,inform)
      if ( inform .lt. 0 ) return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalp(y,rho,lambda,equatn,p)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision lambda,p,rho,y

      if ( equatn ) then
          p  = y * ( lambda + 0.5d0 * rho * y )
      else
          if ( lambda + rho * y .ge. 0.0d0 ) then
              p = y * ( lambda + 0.5d0 * rho * y )
          else
              p = - 0.5d0 * lambda ** 2 / rho
          end if
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evaldpdy(y,rho,lambda,equatn,dpdy)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision y,rho,lambda,dpdy

      if ( equatn ) then
          dpdy = lambda + rho * y
      else
          dpdy = max( 0.0d0, lambda + rho * y )
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnalu(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical ignlin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

C     This subroutine computes the gradient of the Augmented Lagrangian 
C     function at a point xp, which is near to x, taking care of the
C     non-differentiability. The Augmented Lagrangian gradient must be
C     previously computed at x.

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

      if ( innercall ) then
          call minsqg(n,xp,nalp,inform)
          return
      end if

      call ssetp(n,xp)

      if ( gjaccoded ) then
          call ievalnalub(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +    inform)
      else
          call ievalnalua(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +    inform)
      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnalub(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical ignlin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer i,j,jcpnnz
      double precision dpdcp,dum

C     LOCAL ARRAYS
      integer jcpfun(jcnnzmax),jcplen(mmax),jcpsta(mmax),
     +        jcpvar(jcnnzmax)
      double precision cp(mmax),jcpval(jcnnzmax)

C     COMPUTE CONSTRAINTS AT xp

      if ( m .gt. 0 ) then
          call sevalfc(n,xp,dum,m,cp,inform)
          if ( inform .lt. 0 ) return
      end if

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AND JACOBIAN
C     OF CONSTRAINTS AT xp

      call sevalgjac(n,xp,nalp,m,jcpfun,jcpvar,jcpval,jcpnnz,inform)
      if ( inform .lt. 0 ) return

C     CONVERT JACOBIAN OF CONSTRAINTS FROM COODINATE FORMAT TO
C     COMPRESSED SPARSE ROW FORMAT

      call coo2csr(m,jcpnnz,jcpfun,jcpvar,jcpval,jcplen,jcpsta)

C     COMPUTE \nabla L = \nabla f + \sum_j dPdc * dcdx

      do j = 1,m

          if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +         .not. ( ignlin .and. linear(j) ) ) then

C             COMPUTE dP/dc

              dpdcp = lambda(j) + rho(j) * cp(j)

              if ( dpdcp .ne. 0.0d0 ) then

C                 ADD dPdc * dcdx

                  do i = jcpsta(j),jcpsta(j) + jcplen(j) - 1
                      nalp(jcpvar(i)) = 
     +                nalp(jcpvar(i)) + dpdcp * jcpval(i)
                  end do

              end if

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalnalua(n,xp,m,lambda,rho,equatn,linear,ignlin,nalp,
     +inform)

      implicit none

C     This subroutine computes the gradient of the Augmented Lagrangian 
C     function at a point xp, which is near to x, taking care of the
C     non-differentiability. The Augmented Lagrangian gradient must be
C     previously computed at x.

C     SCALAR ARGUMENTS
      logical ignlin
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer i,j,jcpnnz
      double precision cp,dpdcp

C     LOCAL ARRAYS
      integer jcpvar(nmax)
      double precision jcpval(nmax)

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AT xp

      call sevalg(n,xp,nalp,inform)
      if ( inform .lt. 0 ) return

C     COMPUTE \nabla L = \nabla f + \sum_j dPdc * dcdx

      do j = 1,m

          if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and.
     +         .not. ( ignlin .and. linear(j) ) ) then

C             COMPUTE THE i-TH CONSTRAINT
              call sevalc(n,xp,j,cp,inform)
              if ( inform .lt. 0 ) return

C             COMPUTE dP/dc
              dpdcp = lambda(j) + rho(j) * cp

              if ( dpdcp .ne. 0.0d0 ) then

C                 COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
                  call sevaljac(n,xp,j,jcpvar,jcpval,jcpnnz,inform)
                  if ( inform .lt. 0 ) return

C                 ADD dPdc * dcdx
                  do i = 1,jcpnnz
                      nalp(jcpvar(i)) = 
     +                nalp(jcpvar(i)) + dpdcp * jcpval(i)
                  end do

              end if

          end if

      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhal(n,x,m,lambda,rho,equatn,linear,hallin,halcol,
     +halval,halnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer halnnz,inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer hallin(*),halcol(*)
      double precision lambda(m),halval(*),rho(m),x(n)

C     This subroutine computes the Hessian of the augmented Lagrangian.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer i,lin,j,k,l,var

C     LOCAL ARRAYS
      integer stlin(nmax)
      double precision r(nmax)

      call sevalhl(n,x,m,dpdc,hallin,halcol,halval,halnnz,inform)
      if ( inform .lt. 0 ) return

      if ( m .eq. 0 ) return

C     PUT MATRIX INTO A ROW-LINKED LIST

      do i = 1,n
         r(i) = 0.0d0
      end do

      do i = 1,n
          stlin(i) = 0
      end do
      
      do i = 1,halnnz
          lin = hallin(i)
          k   = stlin(lin)
          stlin(lin) = i
          hallin(i)  = k
      end do

C     ADD \sum_j \rho_j * \nabla c_j \nabla c_j^t

      do j = 1,m
         
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             ADD \rho_j * \nabla c_j \nabla c_j^t

              do k = jcsta(j),jcsta(j) + jclen(j) - 1
               
                  var = jcvar(k)

C                 PUT ROW jcvar(k) INTO A DENSE VECTOR

                  lin = stlin(var)
 10               if ( lin .ne. 0 ) then
                      r(halcol(lin)) = r(halcol(lin)) + halval(lin)
                      lin = hallin(lin)
                      go to 10
                  end if
               
C                 ADD VALUE

                  do l = jcsta(j),jcsta(j) + jclen(j) - 1
                      if ( jcvar(l) .le. var ) then
                          r(jcvar(l)) = 
     +                    r(jcvar(l)) + rho(j) * jcval(l) * jcval(k)
                      end if
                  end do
               
C                 UPDATE VALUES IN HALVAL
               
                  lin = stlin(var)
 20               if ( lin .ne. 0 ) then
                      halval(lin) = r(halcol(lin))
                      r(halcol(lin)) = 0.0d0
                      lin = hallin(lin)
                      go to 20
                  end if
               
C                 INSERT NEW ELEMENTS IN HAL REPRESENTATION

                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      l = jcvar(i)
                      if ( r(l) .ne. 0.0d0 ) then
                          halnnz = halnnz + 1
                          halval(halnnz) = r(l)
                          halcol(halnnz) = l
                          hallin(halnnz) = stlin(var)
                          stlin(var) = halnnz
                          r(l) = 0.0d0
                      end if
                  end do
 
              end do

          end if

      end do
  
C     PUT MATRIX BACK INTO COORDINATE SQUEME

      do i = 1,n

          lin = stlin(i)
 30       if ( lin .ne. 0 ) then
              k = hallin(lin)
              hallin(lin) = i
              lin = k
              go to 30
          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine sevalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,gothl,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(n),lambda(m),p(n),rho(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

C     LOCAL SCALARS
      integer i,j
      double precision atp

      if ( innercall ) then
          call minsqhp(n,x,p,hp,gothl,inform)
          return
      end if

C     --------------------------------------------------------------
C     Hessian approximation
C     --------------------------------------------------------------

      if ( hptype .eq. 'HAPPRO' .and. constrc ) then

          call applyhapp(n,m,rho,equatn,gothl,p,hp)

C     --------------------------------------------------------------
C     Incremental quotients
C     --------------------------------------------------------------

      else if ( hptype .eq. 'INCQUO' .or. hptype .eq. 'HAPPRO' ) then

          call ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)
          if ( inform .lt. 0 ) return

C     --------------------------------------------------------------
C     True Hessian
C     --------------------------------------------------------------

      else if ( hptype .eq. 'TRUEHL' ) then

C         Compute Hessian of Lagrangian times p using dpdc 
C         instead of lambda

          call sevalhlp(n,x,m,dpdc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

C         Add rho A^T A

          do j = 1,m

              if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

                  atp = 0.0d0
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      atp = atp + jcval(i) * p(jcvar(i))
                  end do

                  atp = atp * rho(j)

                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
                  end do

              end if

          end do

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine ievalhalp(n,x,m,lambda,rho,equatn,linear,p,hp,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hp(n),lambda(m),p(n),rho(m),x(n)

C     Computes an approximation of the product of the Hessian of the 
C     Augmented Lagrangian times a vector using incremental quotients.

C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     LOCAL SCALARS
      integer i,j
      double precision atp,psupn,step,xsupn

C     LOCAL  ARRAYS
      double precision xp(nmax)

C     --------------------------------------------------------------
C     Set auxiliary point
C     --------------------------------------------------------------

      xsupn = 0.0d0
      psupn = 0.0d0
      do i = 1,n
          xsupn = max( xsupn, abs( x(i) ) )
          psupn = max( psupn, abs( p(i) ) )
      end do

      step = macheps12 * max( xsupn / psupn, 1.0d0 )

      do i = 1,n
          xp(i) = x(i) + step * p(i)
      end do

C     --------------------------------------------------------------
C     Compute gradient of the augmented Lagrangian at xp considering 
C     the same constraints considered at x and ignoring linear 
C     constraints
C     --------------------------------------------------------------

      call ievalnalu(n,xp,m,lambda,rho,equatn,linear,.true.,hp,inform)
      if ( inform .lt. 0 ) return

C     --------------------------------------------------------------
C     Compute gradients difference
C     --------------------------------------------------------------

      do i = 1,n
          hp(i) = ( hp(i) - g(i) ) / step
      end do

C     --------------------------------------------------------------
C     Add contribution of linear constraints
C     --------------------------------------------------------------

      do j = 1,m

          if ( ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) .and. 
     +         linear(j) ) then

C             Compute inner product <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             Add rho * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine coo2csr(m,nnz,alin,acol,aval,alen,asta)

      implicit none

C     SCALAR ARGUMENTS
      integer m,nnz

C     ARRAY ARGUMENTS
      integer acol(nnz),alen(m),alin(nnz),asta(m)
      double precision aval(nnz)

C     This subroutines converts a matrix from coordinate format to
C     compressed sparse row format.

C     LOCAL SCALARS
      integer i,j,col,coltmp,lin,lintmp
      double precision val,valtmp

      do i = 1,m
          alen(i) = 0
      end do

      do i = 1,nnz
          lin = alin(i)
          alen(lin) = alen(lin) + 1
      end do

      asta(1) = 1
      do i = 2,m
          asta(i) = asta(i-1) + alen(i-1)
      end do

      do i = 1,nnz

          val = aval(i)
          col = acol(i)
          lin = alin(i)

          alin(i) = - 1

 10       if ( lin .ge. 0 ) then

              j = asta(lin)
              asta(lin) = j + 1

              valtmp = aval(j)
              coltmp = acol(j)
              lintmp = alin(j)

              aval(j) = val
              acol(j) = col
              alin(j) = - 1

              val = valtmp
              col = coltmp
              lin = lintmp

              go to 10

          end if

      end do

      do i = 1,m
          asta(i) = asta(i) - alen(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      logical function sstop(n,x,m,lambda,rho,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

C     EXTERNAL FUNCTIONS
      logical minsqstop

      if ( innercall ) then
          sstop = minsqstop(n,x,inform)
          return
      end if

      end
C     *****************************************************************
C     *****************************************************************

      subroutine comphapp(n,m,rho,equatn)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision rho(m)

C     This subroutine computes an approximation H of the Hessian of the
C     Augmented Lagrangian following a very simple idea: "discard the 
C     second order terms and then correct the remaining matrix in order 
C     to satisfy a secant equation".
C
C     Hence, H takes the form
C
C     H = B + S + rho A^t A, 
C
C     where S is the spectral correction of (rho A^t A) and B is 
C     the BFGS correction of (S + rho A^t A). More specifically,
C
C     S = hlspg I,
C
C     where
C
C     hlspg = max(lspgmi, min(lspgma, s^t (y - rho A^t A s) / s^t s)),
C
C     D = S + rho A^t A,
C
C     and
C
C     B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ].
C
C     Note that this subroutine does not compute matrix H explicitly,
C     but computes some quantities that will be used later to compute 
C     the product of H by a vector p.
C
C     The quantities computed by this subroutine are:
C
C     (a) hlspg = s^t (y - rho A^t A s) / (s^t s)
C
C     (b) hds = D s = ( hlspg I + rho A^t A ) s, and
C
C     (c) hstds = <s,hds>.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/
C     COMMON SCALARS
      double precision hlspg,hstds

C     COMMON ARRAYS
      double precision hds(nmax)

C     COMMON BLOCKS
      common /happdata/ hlspg,hds,hstds
      save   /happdata/

C     LOCAL SCALARS
      integer i,j
      double precision ats

C     ------------------------------------------------------------------
C     Compute hds = rho A^t A s
C     ------------------------------------------------------------------

      do i = 1,n
          hds(i) = 0.0d0
      end do

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,s>
              ats = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  ats = ats + jcval(i) * s(jcvar(i))
              end do

              ats = ats * rho(j)

C             ADD rho * ats * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hds(jcvar(i)) = hds(jcvar(i)) + ats * jcval(i) 
              end do

          end if

      end do

      hstds = 0.0d0
      do i = 1,n
          hstds = hstds + s(i) * hds(i)
      end do

C     ------------------------------------------------------------------
C     Compute hlspg = s^t (y - rho A^t A s) / (s^t s)
C     ------------------------------------------------------------------

      if ( sty - hstds .le. 0.0d0 ) then
          hlspg = lspgmi
      else
          hlspg = max( lspgmi, min( (sty - hstds) / sts, lspgma ) )
      end if

      do i = 1,n
          hds(i) = hds(i) + hlspg * s(i)
      end do

      hstds = hstds + hlspg * sts

      end

C     *****************************************************************
C     *****************************************************************

      subroutine applyhapp(n,m,rho,equatn,goth,p,hp)

      implicit none

C     SCALAR ARGUMENTS
      logical goth
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision hp(n),p(n),rho(m)

C     This subroutine computes the product of the matrix computed by 
C     subroutine comphapp times vector p.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/
C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/
C     COMMON SCALARS
      double precision hlspg,hstds

C     COMMON ARRAYS
      double precision hds(nmax)

C     COMMON BLOCKS
      common /happdata/ hlspg,hds,hstds
      save   /happdata/

C     LOCAL SCALARS
      integer i,j
      double precision atp,c1,c2,ptds,pty

C     ------------------------------------------------------------------
C     Compute Hessian approximation
C     ------------------------------------------------------------------

      if ( .not. goth ) then

          goth = .true.
          call comphapp(n,m,rho,equatn)

      end if

C     ------------------------------------------------------------------
C     Compute ( hlspg I ) p
C     ------------------------------------------------------------------

      do i = 1,n
          hp(i) = hlspg * p(i)
      end do

C     ------------------------------------------------------------------
C     Add ( rho A^T A ) p
C     ------------------------------------------------------------------

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             ADD rho * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

C     ------------------------------------------------------------------
C     Add B p, 
C     where B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ]
C     ------------------------------------------------------------------

      if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

          pty = 0.0d0
          ptds = 0.0d0
          do i = 1,n
              pty = pty + p(i) * y(i)
              ptds = ptds + p(i) * hds(i)
          end do

          c1 = pty / sty
          c2 = ptds / hstds
          do i = 1,n
              hp(i) = hp(i) + c1 * y(i) - c2 * hds(i)
          end do

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine comphpre(n,m,rho,equatn)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision rho(m)

C     Consider the preconditioner 
C
C         P = Q + E + diag(rho A^t A) 
C
C     for matrix 
C
C         H = B + S + rho A^t A, 
C
C     where E is the spectral correction of diag(rho A^t A) and Q is 
C     the BFGS correction of (E + diag(rho A^t A)), while S and B are
C     the spectral and BFGS corrections of matrix (rho A^t A), 
C     respectively. 
C
C     This subroutine computes:
C
C     (a) pdiag = diag(rho A^t A),
C
C     (b) plspg such that E = plspg I, 
C
C     (c) psmdy = s - D^-1 y, where D = E + diag(rho A^t A), and
C
C     (d) the inner product psmdty = <psmdy,y>.
C
C     These quantities will be used later, in subroutine applyp, to
C     compute z = P^{-1} r.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     PARAMETERS

C     General constants

      double precision fmin

      parameter ( fmin       = - 1.0d+20 )

C     Line search constants

      integer maxextrap,mininterp
      double precision beta,etaint,etaext,gamma,sigma1,sigma2

      parameter ( maxextrap  =       100 )
      parameter ( mininterp  =         4 )

      parameter ( gamma      =   1.0d-04 )
      parameter ( beta       =     0.5d0 )
      parameter ( sigma1     =     0.1d0 )
      parameter ( sigma2     =     0.9d0 )
      parameter ( etaint     =     2.0d0 )
      parameter ( etaext     =     2.0d0 )

C     Safeguarding spectral step constants

      double precision lspgma,lspgmi

      parameter ( lspgma     =   1.0d+10 )
      parameter ( lspgmi     =   1.0d-10 )

C     Conjugate gradients constants

      integer maxcgitnp
      double precision epsnqmp,theta

      parameter ( theta      =   1.0d-06 )
      parameter ( epsnqmp    =   1.0d-08 )
      parameter ( maxcgitnp  =         5 )

C     BETRA constants

      logical extrp2,extrp4,extrp5
      integer msmaxit
      double precision mseps,msrho,mssig,phieps,trdelini,trdelmin,
     +        tralpha

      parameter ( extrp2     =   .true. )
      parameter ( extrp4     =   .true. )
      parameter ( extrp5     =   .true. )

      parameter ( msmaxit    =       20 )

      parameter ( mseps      =  1.0d-08 )    
      parameter ( mssig      =  1.0d-01 )
      parameter ( msrho      =  9.0d-01 )
      parameter ( phieps     =  1.0d-08 )
      parameter ( trdelini   =  1.0d+02 )
      parameter ( trdelmin   =  1.0d-08 )
      parameter ( tralpha    =  1.0d-01 )

C     GENCAN constants

      integer maxinnitnp
      double precision cggpnf,cgepsi,cgepsf,delmin,eta

      parameter ( maxinnitnp =         5 )
      parameter ( delmin     =   1.0d+04 )
      parameter ( eta        =   0.1d0   )
      parameter ( cggpnf     =   1.0d-04 )
      parameter ( cgepsi     =   1.0d-01 )
      parameter ( cgepsf     =   1.0d-08 )

C     ALGENCAN constants

      logical rhoauto,rhoiden
      integer maxoutitnp,maxoutit
      double precision lammax,lammin,rhofrac,rhomax,rhomult

      parameter ( maxoutitnp =        10 )
      parameter ( maxoutit   =        50 )

      parameter ( lammin     = - 1.0d+20 )
      parameter ( lammax     =   1.0d+20 )

      parameter ( rhoauto    =    .true. )
      parameter ( rhoiden    =    .true. )
      parameter ( rhofrac    =   0.5d0   )
      parameter ( rhomult    =   1.0d+01 )
      parameter ( rhomax     =   1.0d+20 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/

C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/
C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/
C     COMMON SCALARS
      double precision plspg,psmdyty

C     COMMON ARRAYS
      double precision pdiag(nmax),psmdy(nmax)

C     COMMON BLOCKS
      common /hpredata/ pdiag,psmdy,plspg,psmdyty
      save   /hpredata/

C     LOCAL SCALARS
      integer i,j
      double precision sttmp

C     ------------------------------------------------------------------
C     Compute diag( rho A^t A )
C     ------------------------------------------------------------------

      do i = 1,n
          pdiag(i) = 0.0d0
      end do

      do j = 1,m
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  pdiag(jcvar(i)) = 
     +            pdiag(jcvar(i)) + rho(j) * jcval(i) ** 2
              end do
          end if
      end do

C     ------------------------------------------------------------------
C     Compute plspg = s^t (y - diag( rho A^t A ) s) / (s^t s)
C     ------------------------------------------------------------------

      sttmp = 0.0d0
      do i = 1,n
          sttmp = sttmp + pdiag(i) * s(i) ** 2
      end do

      if ( sty - sttmp .le. 0.0d0 ) then
          plspg = lspgmi
      else
          plspg = max( lspgmi, min( ( sty - sttmp ) / sts, lspgma ) )
      end if

C     ------------------------------------------------------------------
C     Compute the BFGS correction Q of ( E + diag( rho A^t A ) )
C
C     Q = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
C         [ <s - D^-1 y, y> s s^t ] / (s^t y)^2,
C
C     where D = ( E + diag( rho A^t A ) )
C     ------------------------------------------------------------------
 
      if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

          psmdyty = 0.0d0
          do i = 1,n
              psmdy(i) = s(i) - y(i) / ( plspg + pdiag(i) )
              psmdyty = psmdyty + psmdy(i) * y(i)
          end do

      end if

      end

C     *****************************************************************
C     *****************************************************************

      subroutine applyhpre(n,m,rho,equatn,gotp,r,z)

      implicit none

C     SCALAR ARGUMENTS
      logical gotp
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m)
      double precision r(n),rho(m),z(n)

C     This subroutine computes the product of the inverse of the matrix
C     computed by subroutine comphpre times vector r, i.e., z = P^{-1} r.

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/
C     COMMON SCALARS
      double precision seucn,sts,sty,yeucn

C     COMMON ARRAYS
      double precision s(nmax),y(nmax)

C     COMMON BLOCKS
      common /sydata/ s,y,seucn,yeucn,sts,sty
      save   /sydata/
C     COMMON SCALARS
      logical sameface
      integer ittype

C     COMMON BLOCKS
      common /itedat/ sameface,ittype
      save   /itedat/
C     COMMON SCALARS
      double precision plspg,psmdyty

C     COMMON ARRAYS
      double precision pdiag(nmax),psmdy(nmax)

C     COMMON BLOCKS
      common /hpredata/ pdiag,psmdy,plspg,psmdyty
      save   /hpredata/

C     LOCAL SCALARS
      integer i
      double precision c1,c2,psmdytr,str

C     ------------------------------------------------------------------
C     Compute P
C     ------------------------------------------------------------------

      if ( .not. gotp ) then

          gotp = .true.
          call comphpre(n,m,rho,equatn)

      end if

C     ------------------------------------------------------------------
C     Compute ( E + diag( rho A^T A ) )^{-1} r
C     ------------------------------------------------------------------

      do i = 1,n
          z(i) = r(i) / ( plspg + pdiag(i) )
      end do

C     ------------------------------------------------------------------
C     Add Q^{-1} r, where
C
C     Q^{-1} = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
C              [ <s - D^-1 y, y> s s^t ] / (s^t y)^2
C
C     and D = ( E + diag( rho A^T A ) )
C     ------------------------------------------------------------------

      if ( sameface .and. sty .gt. macheps12 * seucn * yeucn ) then

          str = 0.0d0
          psmdytr = 0.0d0
          do i = 1,n
              str = str + s(i) * r(i)
              psmdytr = psmdytr + psmdy(i) * r(i)
          end do

          c1 = str / sty
          c2 = psmdytr / sty - psmdyty * str / sty ** 2

          do i = 1,n
              z(i) = z(i) + c1 * psmdy(i) + c2 * s(i) 
          end do

      end if

      end
C     SCALE OBJECTIVE FUNCTION AND CONSTRAINTS

C     ******************************************************************
C     ******************************************************************

      subroutine sinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

C     LOCAL SCALARS
      integer i,fun,j,jcnnz,nbds,neq
      double precision scmax

C     LOCAL ARRAYS
      integer jcfun(jcnnzmax),jcvar(jcnnzmax)
      double precision g(nmax),jcval(jcnnzmax)

      neq = 0
      do j = 1,m
          if ( equatn(j) ) neq = neq + 1
      end do

      nbds = 0
      do i = 1,n
          if ( l(i) .gt. - 1.0d+20 ) nbds = nbds + 1
          if ( u(i) .lt.   1.0d+20 ) nbds = nbds + 1
      end do

      if ( iprintctl(2) ) then
          write(* ,100) n,neq,m-neq,nbds
          write(10,100) n,neq,m-neq,nbds
      end if

      call tinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

C     Write classification line of final model

      if ( iprintctl(6) ) then
          open(50,file='class-tabline.out')
          write(50,400) n,neq,m-neq,nbds
          close(50)
      end if

C     Scaling

      usf = 1.0d0
      do j = 1,m
          usc(j) = 1.0d0
      end do

      if ( scale ) then

          if ( m .eq. 0 ) then
              sf = 1.0d0
              if ( iprintctl(2) ) then
                  write(* ,200) 1.0d0 / sf
                  write(10,200) 1.0d0 / sf
              end if

              return
          end if

          call tsetp(n,x)

          if ( gjaccoded ) then

              call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
              if ( inform .lt. 0 ) return

C             Scale constraints

              do j = 1,m
                  sc(j) = 1.0d0
              end do

              do i = 1,jcnnz
                  fun = jcfun(i)

                  sc(fun) = max( sc(fun), abs( jcval(i) ) )
              end do

          else

              call tevalg(n,x,g,inform)
              if ( inform .lt. 0 ) return

C             Scale constraints

              do j = 1,m
                  sc(j) = 1.0d0

                  call tevaljac(n,x,j,jcvar,jcval,jcnnz,inform)
                  if ( inform .lt. 0 ) return

                  do i = 1,jcnnz
                      sc(j) = max( sc(j), abs( jcval(i) ) )
                  end do
              end do

          end if

C         Scale objective function

          sf = 1.0d0
          do i = 1,n
              sf = max( sf, abs( g(i) ) )
          end do

C         Report scaling factors

          scmax = 0.0d0
          do j = 1,m
              scmax = max( scmax, sc(j) )
          end do

          if ( iprintctl(2) ) then
              write(* ,300) 1.0d0 / sf,1.0d0 / scmax
              write(10,300) 1.0d0 / sf,1.0d0 / scmax
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of variables               : ',I7,
     +       /,1X,'Number of equality constraints    : ',I7,
     +       /,1X,'Number of inequality constraints  : ',I7,
     +       /,1X,'Number of bound constraints       : ',I7)

 200  format(/,1X,'Objective function scale factor   : ',1P,D7.1,
     +       /,1X,'The scaling feature was mainly developed for ',
     +            'constrained problems. For',/,1X,'unconstrained and ',
     +            'bound-constrained problem, please, set the ',
     +            'optimality',/,1X,'tolerance (related to the ',
     +            'sup-norm of the projected gradient of the',/,1X,
     +            'objective function) with a convenient value.')

 300  format(/,1X,'Objective function scale factor   : ',1P,D7.1,
     +       /,1X,'Smallest constraints scale factor : ',1P,D7.1)

 400  format(  1X,I6,1X,I6,1X,I6,1X,I6)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

      if ( scale ) then
          do i = 1,m
              lambda(i) = lambda(i) * sf / sc(i)
          end do

          scale = .false.
      end if

      call tendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

      call tevalf(n,x,f,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) f = f / sf

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

      call tevalg(n,x,g,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,n
              g(i) = g(i) / sf
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

      call tevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,hnnz
              hval(i) = hval(i) / sf
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

      call tevalc(n,x,ind,c,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) c = c / sc(ind)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

      call tevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,jcnnz
              jcval(i) = jcval(i) / sc(ind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

      call tevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,hnnz
              hval(i) = hval(i) / sc(ind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhl(n,x,m,lambda,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

      if ( scale ) then
          call tevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call tevalhl(n,x,m,lambda,usf,usc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalhlp(n,x,m,lambda,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

      if ( scale ) then
          call tevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          call tevalhlp(n,x,m,lambda,usf,usc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer j

      call tevalfc(n,x,f,m,c,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          f = f / sf

          do j = 1,m
              c(j) = c(j) / sc(j)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine sevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical scale
      double precision sf,usf

C     COMMON ARRAYS
      double precision sc(mmax),usc(mmax)
 
C     COMMON BLOCK
      common /scadat/ sc,usc,sf,usf,scale
      save   /scadat/

C     LOCAL SCALARS
      integer i

      call tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
      if ( inform .lt. 0 ) return

      if ( scale ) then
          do i = 1,n
              g(i) = g(i) / sf
          end do

          do i = 1,jcnnz
              jcval(i) = jcval(i) / sc(jcfun(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ssetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

      call tsetp(n,x)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sunsetp()

      implicit none

      call tunsetp()

      end
C     ADD SLACKS

C     ******************************************************************
C     ******************************************************************

      subroutine tinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/

C     LOCAL SCALARS
      integer j
      double precision dum

C     LOCAL ARRAYS
      double precision c(mmax)

      call uinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

      if ( slacks ) then

          nws = n

          call usetp(nws,x)

          if ( fccoded ) then
              call uevalfc(nws,x,dum,m,c,inform)
              if ( inform .lt. 0 ) return

          else
              do j = 1,m
                  if ( .not. equatn(j) ) then
                      call uevalc(nws,x,j,c(j),inform)
                      if ( inform .lt. 0 ) return
                  end if
              end do
          end if

          do j = 1,m
              if ( equatn(j) ) then
                  slaind(j) = - 1
              else
                  equatn(j) = .true.

                  n         = n + 1
                  slaind(j) = n

                  l(n)      = - 1.0d+20
                  u(n)      =   0.0d0
                  x(n)      = max( l(n), min( c(j), u(n) ) )
              end if
          end do

          if ( n .eq. nws ) slacks = .false.

          if ( iprintctl(2) ) then
              write(* ,100) n - nws
              write(10,100) n - nws
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of added slack variables   : ',I7)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer j

      if  ( slacks ) then
          n = nws

          do j = 1,m
              if ( slaind(j) .ne. - 1 ) then
                  equatn(j) = .false.
              end if
          end do

          slacks = .false.
      end if

      call uendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

      if ( .not. slacks ) then
          call uevalf(n,x,f,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalf(nws,x,f,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer i

      if ( .not. slacks ) then
          call uevalg(n,x,g,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalg(nws,x,g,inform)
          if ( inform .lt. 0 ) return

          do i = nws + 1,n
              g(i) = 0.0d0
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

      if ( .not. slacks ) then
          call uevalh(n,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalh(nws,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer sind

      if ( .not. slacks ) then
          call uevalc(n,x,ind,c,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalc(nws,x,ind,c,inform)
          if ( inform .lt. 0 ) return

          sind = slaind(ind)
          if ( sind .ne. - 1 ) c = c - x(sind)
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer sind

      if ( .not. slacks ) then
          call uevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevaljac(nws,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          sind = slaind(ind)
          if ( sind .ne. - 1 ) then
              jcnnz = jcnnz +  1
              jcvar(jcnnz) = sind
              jcval(jcnnz) = - 1.0d0
          end if
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

      if ( .not. slacks ) then
          call uevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalhc(nws,x,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

      if ( .not. slacks ) then
          call uevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalhl(nws,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer i

      if ( .not. slacks ) then
          call uevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalhlp(nws,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

          do i = nws + 1,n
              hp(i) = 0.0d0
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer j,sind

      if ( .not. slacks ) then
          call uevalfc(n,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalfc(nws,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

          do j = 1,m
              sind = slaind(j)
              if ( sind .ne. - 1 ) c(j) = c(j) - x(sind)
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine tevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

C     LOCAL SCALARS
      integer i,j,sind

      if ( .not. slacks ) then
          call uevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          call uevalgjac(nws,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = nws + 1,n
              g(i) = 0.0d0
          end do

          do j = 1,m
              sind = slaind(j)
              if ( sind .ne. - 1 ) then
                  jcnnz = jcnnz +  1
                  jcfun(jcnnz) = j
                  jcvar(jcnnz) = sind
                  jcval(jcnnz) = - 1.0d0
              end if
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tsetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical slacks
      integer nws

C     COMMON ARRAYS
      integer slaind(mmax)

C     COMMON BLOCKS
      common /sladat/ slaind,nws,slacks
      save   /sladat/

      if ( .not. slacks ) then
          call usetp(n,x)

      else
          call usetp(nws,x)
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine tunsetp()

      implicit none

      call uunsetp()

      end
C     ADD REMOVED FIXED VARIABLES

C     ******************************************************************
C     ******************************************************************

      subroutine uinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(10),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vinip

      call vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,inform)
      if ( inform .lt. 0 ) return

C     Eliminate fixed variables (l=u) and save their values on y

      if ( rmfixv ) then

          yind(0) = n

          n = 0
          do i = 1,yind(0)
              if ( l(i) .lt. u(i) ) then
                  n = n + 1
                  yind(n) = i
                  ycor(i) = n
              else
                  y(i) = l(i)
                  ycor(i) = 0
              end if
          end do

          do i = 1,n
              x(i) = x(yind(i))
              l(i) = l(yind(i))
              u(i) = u(yind(i))
          end do

          if ( n .eq. yind(0) ) rmfixv = .false.

          if ( iprintctl(2) ) then
              write(* ,100) yind(0) - n
              write(10,100) yind(0) - n
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'Number of removed fixed variables : ',I7)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i,ind

C     EXTERNAL SUBROUTINES
      external vendp

C     Restore original x, l, u and n

      if ( rmfixv ) then
          do i = yind(0),1,-1
              ind = ycor(i)
              if ( ind .ne. 0 ) then
                  l(i) = l(ind)
                  u(i) = u(ind)
                  x(i) = x(ind)
              else
                  l(i) = y(i)
                  u(i) = y(i)
                  x(i) = y(i)
              end if
          end do

          n = yind(0)

          rmfixv = .false.
      end if

      call vendp(n,x,l,u,m,lambda,equatn,linear,inform)
      if ( inform .lt. 0 ) return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalf

      if ( .not. rmfixv ) then
          call vevalf(n,x,f,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalf: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalf(yind(0),y,f,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalg

      if ( .not. rmfixv ) then
          call vevalg(n,x,g,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalg: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalg(yind(0),y,g,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              g(i) = g(yind(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer col,i,j,lin

C     EXTERNAL SUBROUTINES
      external vevalh

      if ( .not. rmfixv ) then
          call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalh: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalh(yind(0),y,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,hnnz
              lin = ycor(hlin(i))
              col = ycor(hcol(i))
              if ( lin .ne. 0 .and. col .ne. 0 ) then
                  j = j + 1
                  hlin(j) = lin
                  hcol(j) = col
                  hval(j) = hval(i)
              end if
          end do

          hnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalc

      if ( .not. rmfixv ) then
          call vevalc(n,x,ind,c,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalc: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalc(yind(0),y,ind,c,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i,j,var

C     EXTERNAL SUBROUTINES
      external vevaljac

      if ( .not. rmfixv ) then
          call vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevaljac: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevaljac(yind(0),y,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,jcnnz
              var = ycor(jcvar(i))
              if ( var .ne. 0 ) then
                  j = j + 1
                  jcvar(j) = var
                  jcval(j) = jcval(i)
              end if
          end do

          jcnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer col,i,j,lin

C     EXTERNAL SUBROUTINES
      external vevalhc

      if ( .not. rmfixv ) then
          call vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalhc: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalhc(yind(0),y,ind,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,hnnz
              lin = ycor(hlin(i))
              col = ycor(hcol(i))
              if ( lin .ne. 0 .and. col .ne. 0 ) then
                  j = j + 1
                  hlin(j) = lin
                  hcol(j) = col
                  hval(j) = hval(i)
              end if
          end do

          hnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer col,i,j,lin

C     EXTERNAL SUBROUTINES
      external vevalhl

      if ( .not. rmfixv ) then
          call vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalhl: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalhl(yind(0),y,m,lambda,sf,sc,hlin,hcol,hval,hnnz,
     +    inform)
          if ( inform .lt. 0 ) return

          j = 0
          do i = 1,hnnz
              lin = ycor(hlin(i))
              col = ycor(hcol(i))
              if ( lin .ne. 0 .and. col .ne. 0 ) then
                  j = j + 1
                  hlin(j) = lin
                  hcol(j) = col
                  hval(j) = hval(i)
              end if
          end do

          hnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i

C     LOCAL ARRAYS
      double precision w(nmax)

C     EXTERNAL SUBROUTINES
      external vevalhlp

      if ( .not. rmfixv ) then
          call vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevalhlp: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          do i = 1,yind(0)
              w(i) = 0.0d0
          end do

          do i = 1,n
              w(yind(i)) = p(i)
          end do

          call vevalhlp(yind(0),y,m,lambda,sf,sc,w,hp,gothl,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              hp(i) = hp(yind(i))
          end do
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vevalfc

      if ( .not. rmfixv ) then
          call vevalfc(n,x,f,m,c,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevaljac: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalfc(yind(0),y,f,m,c,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine uevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i,j,var

C     EXTERNAL SUBROUTINES
      external vevalgjac

      if ( .not. rmfixv ) then
          call vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

      else
          if ( .not. yset ) then
              write(*,*) 'uevaljac: Opa!!!!!!!!!!!!!!!!!!!!!!!!!'
              do i = 1,n
                  y(yind(i)) = x(i)
              end do
          end if

          call vevalgjac(yind(0),y,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return

          do i = 1,n
              g(i) = g(yind(i))
          end do

          j = 0
          do i = 1,jcnnz
              var = ycor(jcvar(i))
              if ( var .ne. 0 ) then
                  j = j + 1
                  jcfun(j) = jcfun(i)
                  jcvar(j) = var
                  jcval(j) = jcval(i)
              end if
          end do

          jcnnz = j
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine usetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     LOCAL SCALARS
      integer i

C     EXTERNAL SUBROUTINES
      external vsetp

      if ( .not. rmfixv ) then
          call vsetp(n,x)
          return
      end if

      yset = .true.

      do i = 1,n
          y(yind(i)) = x(i)
      end do

      call vsetp(yind(0),y)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine uunsetp()

      implicit none

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical rmfixv,yset

C     COMMON ARRAYS
      integer ycor(nmax),yind(0:nmax)
      double precision y(nmax)

C     COMMON BLOCKS
      common /fixvar/ y,ycor,yind,yset,rmfixv
      save   /fixvar/

C     EXTERNAL SUBROUTINES
      external vunsetp

      yset = .false.

      call vunsetp()

      end
C     SAFE CALL TO THE USER PROVIDED SUBROUTINES

C     ******************************************************************
C     ******************************************************************

      subroutine vinip(n,x,l,u,m,lambda,equatn,linear,coded,checkder,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      logical checkder
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical coded(6),equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )

C     LOCAL SCALARS
      integer i

C     INTRINSEC FUNCTIONS
      intrinsic max,min

C     EXTERNAL SUBROUTINES
      external checkd

C     Avoid huge bounds

      do i = 1,n
          l(i) = max( l(i), - 1.0d+20 )
          u(i) = min( u(i),   1.0d+20 )
      end do

C     Project initial guess

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Check derivatives

      if ( checkder ) then
          call checkd(n,l,u,m,inform)
          if ( inform .lt. 0 ) return
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vendp(n,x,l,u,m,lambda,equatn,linear,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer i

      if ( iprintctl(2) ) then

C         Save solution

          open(20,file='solution.txt')

C         Point

          write(20,100)
          do i = 1,n
              write(20,200) i,x(i)
          end do

C         Lagrange multipliers

          if ( m .gt. 0 ) then
              write(20,300)
              do i = 1,m
                  write(20,200) i,lambda(i)
              end do
          end if

          close(20)

      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,'FINAL POINT:',//,2X,'INDEX',16X,'X(INDEX)')
 200  format(I7,1P,D24.16)
 300  format(/,'FINAL ESTIMATION OF THE LAGRANGE MULTIPLIERS: ',
     +       //,2X,'INDEX',11X,'LAMBDA(INDEX)')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalf(n,x,f,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalf,reperr

      call evalf(n,x,f,flag)

      efcnt = efcnt + 1
      fcnt  = fcnt  + 1

      if ( flag .ne. 0 ) then
          inform = - 90
          call reperr(inform)
          return
      end if

      if ( .not. IsANumber(f) ) then
          if ( iprintctl(3) ) then
              write(* ,100)
              write(* ,200) f
              write(10,100)
              write(10,200) f
          end if

          inform = - 90
          call reperr(inform)
          return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALF WARNING: The objective function value ',
     +            'computed by the user-supplied subroutine EVALF is ',
     +            '+Inf, -Inf or NaN.')

 200  format(/,1X,'Value: ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalg,ievalg,reperr

      if ( gcoded ) then
          call evalg(n,x,g,flag)

          egcnt = egcnt + 1

          if ( flag .ne. 0 ) then
              inform = - 92
              call reperr(inform)
              return
          end if

          do i = 1,n
              if ( .not. IsANumber(g(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,200) n,i,g(i)
                      write(10,100)
                      write(10,200) n,i,g(i)

                      inform = - 92
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ievalg(n,x,g,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALG WARNING: There is an element whose value is ',
     +            '+Inf, -Inf or NaN in the gradient of the objective ',
     +            'function computed by the user-supplied subroutine ',
     +            'EVALG.')

 200  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalg(n,x,g,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer j
      double precision fminus,fplus,step,tmp

C     INTRINSIC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external setp,vevalf

      do j = 1,n
          tmp  = x(j)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(j) = tmp + step
          call setp(n,x)
          call vevalf(n,x,fplus,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp - step
          call setp(n,x)
          call vevalf(n,x,fminus,inform)
          if ( inform .lt. 0 ) return

          g(j) = ( fplus - fminus ) / ( 2.0d0 * step )

          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalh(n,x,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalh,reperr

      call evalh(n,x,hlin,hcol,hval,hnnz,flag)

      ehcnt = ehcnt + 1

      if ( flag .ne. 0 ) then
          inform = - 94
          call reperr(inform)
          return
      end if

      do i = 1,hnnz
          if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +         hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. 
     +         hcol(i) .gt. hlin(i) ) then

              if ( iprintctl(3) ) then
                  write(* ,100)
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(10,100)
                  write(10,300) n,i,hlin(i),hcol(i),hval(i)
              end if

              hlin(i) = 1
              hcol(i) = 1
              hval(i) = 0.0d0
          end if

          if ( .not. IsANumber(hval(i)) ) then
              if ( iprintctl(3) ) then
                  write(* ,200) 
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(10,200) 
                  write(10,300) n,i,hlin(i),hcol(i),hval(i)

                  inform = - 94
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALH WARNING: There is an element out of range, ',
     +            'or in the upper triangle, of the Hessian of the ',
     +            'objetive function computed by the user-supplied ',
     +            'subroutine EVALH. It will be ignored.')

 200  format(/,1X,'VEVALH WARNING: There is an element whose value is ',
     +            '+Inf, -Inf or NaN in the Hessian of the objetive ',
     +            'function computed by the user-supplied subroutine ',
     +            'EVALH.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalc(n,x,ind,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision c

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalc,reperr

      call evalc(n,x,ind,c,flag)

      eccnt(ind) = eccnt(ind) + 1

      if ( flag .ne. 0 ) then
          inform = - 91
          call reperr(inform)
          return
      end if

      if ( .not. IsANumber(c) ) then
          if ( iprintctl(3) ) then
              write(* ,100) ind
              write(* ,200) c
              write(10,100) ind
              write(10,200) c

              inform = - 91
              call reperr(inform)
              return
          end if
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALC WARNING: The value of constraint ',I16,' ',
     +            'computed by the user-supplied subroutine EVALC is ',
     +            '+Inf, -Inf or NaN.')

 200  format(/,1X,'Value: ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision x(n),jcval(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evaljac,ievaljac,reperr

      if ( jaccoded ) then
          call evaljac(n,x,ind,jcvar,jcval,jcnnz,flag)

          ejccnt(ind) = ejccnt(ind) + 1

          if ( flag .ne. 0 ) then
              inform = - 93
              call reperr(inform)
              return
          end if

          do i = 1,jcnnz
              if ( jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

                  if ( iprintctl(3) ) then
                      write(* ,100) ind
                      write(* ,300) n,i,jcvar(i),jcval(i)
                      write(10,100) ind
                      write(10,300) n,i,jcvar(i),jcval(i)
                  end if

                  jcvar(i) = 1
                  jcval(i) = 0.0d0
              end if

              if ( .not. IsANumber(jcval(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,200) ind
                      write(* ,300) n,i,jcvar(i),jcval(i)
                      write(10,200) ind
                      write(10,300) n,i,jcvar(i),jcval(i)

                      inform = - 93
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ievaljac(n,x,ind,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALJAC WARNING: There is an element out of ',
     +            'range in the gradient of constraint ',I16,' ',
     +            'computed by the user-supplied subroutine EVALJAC. ',
     +            'It will be ignored.')

 200  format(/,1X,'VEVALJAC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the gradient of constraint ',
     +            I16,'computed by the user-supplied subroutine ',
     +            'EVALJAC.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Variable : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievaljac(n,x,ind,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer ind,inform,n,jcnnz

C     ARRAY ARGUMENTS
      integer jcvar(n)
      double precision jcval(n),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer j
      double precision cminus,cplus,step,tmp

C     INTRINSEC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external setp,vevalc

      jcnnz = 0

      do j = 1,n
          tmp  = x(j)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(j) = tmp + step
          call setp(n,x)
          call vevalc(n,x,ind,cplus,inform)
          if ( inform .lt. 0 ) return

          x(j) = tmp - step
          call setp(n,x)
          call vevalc(n,x,ind,cminus,inform)
          if ( inform .lt. 0 ) return

          jcvar(jcnnz + 1) = j
          jcval(jcnnz + 1) = ( cplus - cminus ) / ( 2.0d0 * step )

          if ( abs( jcval(jcnnz + 1) ) .gt. 0.0d0 ) then
              jcnnz = jcnnz + 1
          end if

          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalhc(n,x,ind,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,ind,n,hnnz

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhc,reperr

      call evalhc(n,x,ind,hlin,hcol,hval,hnnz,flag)

      ehccnt(ind) = ehccnt(ind) + 1

      if ( flag .ne. 0 ) then
          inform = - 95
          call reperr(inform)
          return
      end if

      do i = 1,hnnz
          if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +         hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. 
     +         hcol(i) .gt. hlin(i) ) then

              if ( iprintctl(3) ) then
                  write(* ,100) ind
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(10,100) ind
                  write(10,300) n,i,hlin(i),hcol(i),hval(i)
              end if

              hlin(i) = 1
              hcol(i) = 1
              hval(i) = 0.0d0
          end if

          if ( .not. IsANumber(hval(i)) ) then
              if ( iprintctl(3) ) then
                  write(* ,200) ind
                  write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                  write(10,200) ind
                  write(10,300) n,i,hlin(i),hcol(i),hval(i)

                  inform = - 95
                  call reperr(inform)
                  return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHC WARNING: There is an element out of range ',
     +            'or in the upper triangle of the Hessian of ',
     +            'constraint ',I16,' computed by the user-supplied ',
     +            'subroutine EVALHC. It will be ignored.')

 200  format(/,1X,'VEVALHC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the Hessian of constraint ',
     +            I16,' computed by the user-supplied subroutine ',
     +            'EVALHC.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhl,ievalhl,reperr

      if ( hlcoded ) then

          call evalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,flag)

          ehlcnt = ehlcnt + 1

          if ( flag .ne. 0 ) then
              inform = - 96
              call reperr(inform)
              return
          end if

          do i = 1,hnnz
              if ( hlin(i) .lt. 1 .or. hlin(i) .gt. n .or.
     +             hcol(i) .lt. 1 .or. hcol(i) .gt. n .or. 
     +             hcol(i) .gt. hlin(i) ) then

                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                      write(10,100)
                      write(10,300) n,i,hlin(i),hcol(i),hval(i)
                  end if

                  hlin(i) = 1
                  hcol(i) = 1
                  hval(i) = 0.0d0
              end if

              if ( .not. IsANumber(hval(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,200)
                      write(* ,300) n,i,hlin(i),hcol(i),hval(i)
                      write(10,200)
                      write(10,300) n,i,hlin(i),hcol(i),hval(i)

                      inform = - 96
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else if ( hcoded .and. ( hccoded .or. m .eq. 0 ) ) then
         call ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
         if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHL WARNING: There is an element out of range, ',
     +            'or in the upper triangle, of the',/,1X,'Hessian of ',
     +            'the Lagrangian computed by the user-supplied ',
     +            'subroutine EVALHL. It',/,1X,'will be ignored.')

 200  format(/,1X,'VEVALHL WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the',/,1X,'Hessian of the ',
     +            'Lagrangian computed by the user-supplied ',
     +            'subroutine EVALHL.')

 300  format(/,1X,'Dimension: ',I16,
     +       /,1X,'Position : ',I16,
     +       /,1X,'Row      : ',I16,
     +       /,1X,'Column   : ',I16,
     +       /,1X,'Value    : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz,inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      integer hlin(*),hcol(*)
      double precision hval(*),lambda(m),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )

C     LOCAL SCALARS
      integer col,con,hnnztmp,i,ind,itmp,j,lin,nextj,rnnz
      double precision val

C     LOCAL ARRAYS
      integer hcon(hnnzmax),pos(nmax),rind(nmax),stlin(nmax)
      double precision rval(nmax)

C     EXTERNAL SUBROUTINES
      external vevalh,vevalhc

C     ==================================================================
C     COMPUTE HESSIANS
C     ==================================================================

C     COMPUTE HESSIAN OF THE OBJECTIVE FUNCTION

      call vevalh(n,x,hlin,hcol,hval,hnnz,inform)
      if ( inform .lt. 0 ) return

C     For each element of the Hessian of the objective function,
C     set constraint index as zero
      do i = 1,hnnz
          hval(i) = hval(i) / sf
          hcon(i) = 0
      end do

      if ( m .eq. 0 ) return

C     COMPUTE HESSIANS OF THE CONSTRAINTS

      ind = 0

      do j = 1,m
C         Compute Hessian of constraint j
          call vevalhc(n,x,j,hlin(hnnz+ind+1),hcol(hnnz+ind+1),
     +    hval(hnnz+ind+1),hnnztmp,inform)
          if ( inform .lt. 0 ) return

C         For each element of the Hessian, set constraint as j
          do i = hnnz + ind + 1,hnnz + ind + hnnztmp
              hval(i) = hval(i) / sc(j)
              hcon(i) = j
          end do
         
          ind = ind + hnnztmp
      end do

      if ( ind .eq. 0 ) return

      hnnz = hnnz + ind

C     ==================================================================
C     SET ROW LINKED LISTS
C     ==================================================================

C     Initialize pointers to the first element of each row
      do i = 1,n
         stlin(i) = 0
      end do

C     Set row linked lists
      do i = 1,hnnz
         lin        = hlin(i)
         itmp       = stlin(lin)
         stlin(lin) = i
         hlin(i)    = itmp
      end do

C     ==================================================================
C     BUILD HESSIAN OF THE LAGRANGIAN ROW BY ROW
C     ==================================================================

C     Initialize array pos
      do i = 1,n
          pos(i) = 0
      end do

      do i = 1,n
C         Initialize the i-th row of the Hessian of the Lagrangian
          rnnz = 0

C         Process the i-th row of all the Hessians
          j = stlin(i)

 10       if ( j .ne. 0 ) then

C             Process element (i,hcol(j)) of the Hessian of constraint 
C             hcon(j) (Hessian of the objective function if hcon(j)=0)

              col = hcol(j)
              con = hcon(j)
              if ( con .eq. 0 ) then
                  val = hval(j)
              else
                  val = hval(j) * lambda(con)
              end if

              if ( pos(col) .ne. 0 ) then
                  rval(pos(col)) = rval(pos(col)) + val

              else
                  rnnz           = rnnz + 1
                  pos(col)       = rnnz
                  rind(pos(col)) = col
                  rval(pos(col)) = val
              end if

C             Get next element in the i-th row linked list
              j = hlin(j)
              go to 10
          end if

C         Clean array pos
          do j = 1,rnnz
              pos(rind(j)) = 0
          end do

C         Set i-th row of hl (over the i-th rows of the Hessians) 
C         and mark remaining elements to be deleted
          j = stlin(i)

 20       if ( j .ne. 0 ) then
              nextj = hlin(j)

              if ( rnnz .ne. 0 ) then
                  hlin(j) = i
                  hcol(j) = rind(rnnz)
                  hval(j) = rval(rnnz)
                  rnnz    = rnnz - 1
              else
                  hlin(j) = 0
              end if

              j = nextj
              go to 20
          end if

      end do

C     Eliminate remaining elements (marked with hlin(j)=0)
      j = 1

 30   if ( j .le. hnnz ) then
          if ( hlin(j) .eq. 0 ) then
              if ( j .ne. hnnz ) then
                  hlin(j) = hlin(hnnz)
                  hcol(j) = hcol(hnnz)
                  hval(j) = hval(hnnz)
              end if
              hnnz = hnnz - 1
          else
              j = j + 1
          end if

          go to 30
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUNCTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalhlp,ievalhlp,reperr

      if ( hlpcoded ) then
          call evalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,flag)

          ehlpcnt = ehlpcnt + 1

          if ( flag .ne. 0 ) then
              inform = - 97
              call reperr(inform)
              return
          end if

          do i = 1,n
              if ( .not. IsANumber(hp(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,200) n,i,hp(i)
                      write(10,100)
                      write(10,200) n,i,hp(i)

                      inform = - 97
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else if ( truehl ) then
          call ievalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALHLP WARNING: There is an element in the ',
     +            'product of the Hessian of the Lagrangian by a ',
     +            'computed by the user-supplied subroutine EVALHLP ',
     +            'whose value is +Inf, -Inf or NaN.')

 200  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalhlp(n,x,m,lambda,sf,sc,p,hp,gothl,inform)

      implicit none

C     SCALAR ARGUMENTS
      logical gothl
      integer inform,m,n
      double precision sf

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),sc(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer hnnz

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlin(hnnzmax)
      double precision hval(hnnzmax)

C     COMMON BLOCKS
      common /hdata/ hval,hlin,hcol,hnnz
      save   /hdata/

C     LOCAL SCALARS
      integer col,i,lin
      double precision val

C     EXTERNAL SUBROUTINES
      external vevalhl

      if ( .not. gothl ) then
          gothl = .true.

          call vevalhl(n,x,m,lambda,sf,sc,hlin,hcol,hval,hnnz,inform)
          if ( inform .lt. 0 ) return
      end if

      do i = 1,n
          hp(i) = 0.0d0
      end do

      do i = 1,hnnz
          lin = hlin(i)
          col = hcol(i)
          val = hval(i)

          hp(lin) = hp(lin) + p(col) * val 

          if ( lin .ne. col ) then
              hp(col) = hp(col) + p(lin) * val 
          end if
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vevalfc(n,x,f,m,c,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,m,n
      double precision f

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalfc,reperr

      call evalfc(n,x,f,m,c,flag)

      efccnt = efccnt + 1
      fcnt   = fcnt   + 1

      if ( flag .ne. 0 ) then
          inform = - 98
          call reperr(inform)
          return
      end if

      if ( .not. IsANumber(f) ) then
          if ( iprintctl(3) ) then
              write(* ,100)
              write(* ,300) f
              write(10,100)
              write(10,300) f

c             inform = - 98
c             call reperr(inform)
c             return
          end if
      end if

      do i = 1,m
          if ( .not. IsANumber(c(i)) ) then
              if ( iprintctl(3) ) then
                  write(* ,200)
                  write(* ,400) n,m,i,c(i)
                  write(10,200)
                  write(10,400) n,m,i,c(i)

c                 inform = - 98
c                 call reperr(inform)
c                 return
              end if
          end if
      end do

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALFC WARNING: The objective function value ',
     +            'computed by the user-supplied',/,1X,'subroutine ',
     +            'EVALFC is +Inf, -Inf or NaN.')

 200  format(/,1X,'VEVALFC WARNING: The value of a constraint ',
     +            'computed by the user-supplied',/,1X,'subroutine ',
     +            'EVALFC is +Inf, -Inf or NaN.')

 300  format(/,1X,'Value: ',1P,D24.16)

 400  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Number of constraints : ',I16,
     +       /,1X,'Constraint            : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine vevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      character * 6 hptype
      logical ccoded,fccoded,fcoded,gcoded,gjaccoded,innercall,jaccoded,
     +        hcoded,hccoded,hlcoded,hlpcoded,truehl,ignoref,avoidds,
     +        skipacc,sclsys,useustp

C     COMMON BLOCKS
      common /algparam/ fcoded,gcoded,hcoded,ccoded,jaccoded,hccoded,
     +                  hlcoded,hlpcoded,fccoded,gjaccoded,truehl,
     +                  ignoref,avoidds,skipacc,sclsys,innercall,
     +                  useustp,hptype
      save   /algparam/
C     COMMON SCALARS
      integer efcnt,efccnt,egcnt,egjccnt,ehcnt,ehlcnt,ehlpcnt,fcnt

C     COMMON ARRAYS
      integer eccnt(mmax),ehccnt(mmax),ejccnt(mmax)

C     COMMON BLOCKS
      common /counters/ eccnt,ehccnt,ejccnt,efcnt,efccnt,egcnt,egjccnt,
     +                  ehcnt,ehlcnt,ehlpcnt,fcnt
      save   /counters/
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


C     LOCAL SCALARS
      integer flag,i

C     EXTERNAL FUCNTIONS
      logical IsANumber
      external IsANumber

C     EXTERNAL SUBROUTINES
      external evalgjac,ievalgjac,reperr

      if ( gjaccoded ) then
          call evalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag)

          egjccnt = egjccnt + 1

          if ( flag .ne. 0 ) then
              inform = - 99
              call reperr(inform)
              return
          end if

          do i = 1,jcnnz
              if ( jcfun(i) .lt. 1 .or. jcfun(i) .gt. m .or.
     +             jcvar(i) .lt. 1 .or. jcvar(i) .gt. n ) then

                  if ( iprintctl(3) ) then
                      write(* ,100)
                      write(* ,300) n,m,i,jcfun(i),jcvar(i),jcval(i)
                      write(10,100)
                      write(10,300) n,m,i,jcfun(i),jcvar(i),jcval(i)
                  end if

                  jcfun(i) = 1
                  jcvar(i) = 1
                  jcval(i) = 0.0d0
              end if

              if ( .not. IsANumber(jcval(i)) ) then
                  if ( iprintctl(3) ) then
                      write(* ,200)
                      write(* ,300) n,m,i,jcfun(i),jcvar(i),jcval(i)
                      write(10,200)
                      write(10,300) n,m,i,jcfun(i),jcvar(i),jcval(i)

                      inform = - 99
                      call reperr(inform)
                      return
                  end if
              end if
          end do

      else
          call ievalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)
          if ( inform .lt. 0 ) return
      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'VEVALGJAC WARNING: There is an element out of ',
     +            'range in the gradient of the objective function or ',
     +            'in the Jacobian of the constraints computed by the ',
     +            'user-supplied subroutine EVALGJAC. It will be ',
     +            'ignored.')

 200  format(/,1X,'VEVALGJAC WARNING: There is an element whose value ',
     +            'is +Inf, -Inf or NaN in the',/,1X,'gradient of the ',
     +            'objective function or in the Jacobian of the ',
     +            'constraints',/,1X,'computed by the user-supplied ',
     +            'subroutine EVALGJAC.')

 300  format(/,1X,'Dimension of the space: ',I16,
     +       /,1X,'Number of constraints : ',I16,
     +       /,1X,'Position              : ',I16,
     +       /,1X,'Constraint            : ',I16,
     +       /,1X,'Variable              : ',I16,
     +       /,1X,'Value                 : ',1P,D24.16)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine ievalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer inform,jcnnz,m,n

C     ARRAY ARGUMENTS
      integer jcfun(*),jcvar(*)
      double precision g(n),jcval(*),x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

C     LOCAL SCALARS
      integer i,j
      double precision fminus,fplus,step,tmp

C     LOCAL ARRAYS
      double precision cminus(mmax),cplus(mmax)

C     INTRINSIC FUNCTIONS
      intrinsic abs,max

C     EXTERNAL SUBROUTINES
      external setp,vevalfc

      jcnnz = 0

      do i = 1,n
          tmp  = x(i)

          step = macheps13 * max( 1.0d0, abs( tmp ) )

          x(i) = tmp + step
          call setp(n,x)
          call vevalfc(n,x,fplus,m,cplus,inform)
          if ( inform .lt. 0 ) return

          x(i) = tmp - step
          call setp(n,x)
          call vevalfc(n,x,fminus,m,cminus,inform)
          if ( inform .lt. 0 ) return

          do j = 1,m
              jcfun(jcnnz + 1) = j
              jcvar(jcnnz + 1) = i
              jcval(jcnnz + 1) = ( cplus(j) - cminus(j) ) /
     +                           ( 2.0d0 * step )

              if ( abs( jcval(jcnnz + 1) ) .gt. 0.0d0 ) then
                  jcnnz = jcnnz + 1
              end if
          end do

          g(i) = ( fplus - fminus ) / ( 2.0d0 * step )

          x(i) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine reperr(inform)

C     SCALAR ARGUMENTS
      integer inform

C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/


      if ( iprintctl(3) ) then

          if ( inform .eq. -90 ) then
              write(* ,100) 'EVALF' 
              write(10,100) 'EVALF' 

          else if ( inform .eq. -91 ) then
              write(* ,100) 'EVALC' 
              write(10,100) 'EVALC' 

          else if ( inform .eq. -92 ) then
              write(* ,100) 'EVALG' 
              write(10,100) 'EVALG' 

          else if ( inform .eq. -93 ) then
              write(* ,100) 'EVALJAC' 
              write(10,100) 'EVALJAC' 

          else if ( inform .eq. -94 ) then
              write(* ,100) 'EVALH' 
              write(10,100) 'EVALH' 

          else if ( inform .eq. -95 ) then
              write(* ,100) 'EVALHC' 
              write(10,100) 'EVALHC' 

          else if ( inform .eq. -96 ) then
              write(* ,100) 'EVALHL' 
              write(10,100) 'EVALHL' 

          else if ( inform .eq. -97 ) then
              write(* ,100) 'EVALHLP' 
              write(10,100) 'EVALHLP' 

          else if ( inform .eq. -98 ) then
              write(* ,100) 'EVALFC' 
              write(10,100) 'EVALFC' 

          else if ( inform .eq. -99 ) then
              write(* ,100) 'EVALGJAC' 
              write(10,100) 'EVALGJAC' 
          end if

      end if

C     NON-EXECUTABLE STATEMENTS

 100  format(/,1X,'*** There was an error in the user supplied ',
     +            'subroutine ',A10,' ***',/)

      end

C     ******************************************************************
C     ******************************************************************

      logical function IsANumber(x)

      implicit none

C     SCALAR ARGUMENTS
      double precision x

C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

      IsANumber = .true.

      if ( .not. abs( x ) .le. bignum ) then
          IsANumber = .false.
      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vsetp(n,x)

      implicit none

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     EXTERNAL SUBROUTINES
      external setp

      gotc = .false.

      call setp(n,x)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine vunsetp()

      implicit none

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical constrc,gotc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision c(mmax),dpdc(mmax),g(nmax),jcval(jcnnzmax)
c     double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     COMMON BLOCKS
      common /gdata/ g,c,dpdc,jcval,jcvar,jcsta,jclen,constrc,gotc
c     common /gdata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      save   /gdata/


C     EXTERNAL SUBROUTINES
      external unsetp

      gotc = .false.

      call unsetp()

      end
C     ******************************************************************
C     ******************************************************************

      double precision function drand(ix)

C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision ix

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= ix/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (ix-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (ix.lt.0) ix= ix + p
      drand= ix*4.656612875d-10

      return

      end
C *******************************************************************
C COPYRIGHT (c) 1999 CCLRC Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 13 September 1999
C 01/11/00  Entries in IW initialized to zero in MA57O/OD to avoid copy
C           of unassigned variables by MA57E/ED.
C           AINPUT and IINPUT reset in call to MA57E/ED.
C 06/02/01  Default values for ICNTL(12) and ICNTL(13) changed.
C           Control for direct addressing in solve changed to be
C           on number of rows and number columns in block pivot.
C           Several comments changed as consequence.
C           INFO(31) added to record number of block pivots.
C           Subroutines MA57X/XD and MA57Y/YD added for efficiency when
C           only one rhs (equivalent to MA57Q/QD and MA57R/RD resp).
C 04/07/01  Use of MC41 changed to use of MC71.
C 26/10/01  Printing controls corrected to ensure ICNTL(5) is used and
C           unit number always checked for being positive before
C           printing. Text and comments changed to reflect that D
C           inverse is held in factors and text for solution changed
C           from Right-hand side to solution.
C           Option of choosing two 1 x 1 pivots when 2 x 2 fails
C           removed.
C           MC47B/BD given remaining length in KEEP to avoid compresses
C 20/12/01  INFO(1) initialized to zero in MA57E/ED
C 06/12/02  The test for convergence of iterative refinement changed to
C           avoid any problem with comparisons of numbers held in
C           registers.
C 25/03/03  MC50 (AMD with dense row protection) and MA27 (minimum
C           degree) added. Invoked by ICNTL(6) equal to 2 and 3,
C           respectively. Routines MA57H/HD, MA57V/VD, and MA57Z/ZD
C           have been added to duplicate routines MA27H/HD, MA27G/GD,
C           and MA27U/UD from MA57 and MC50B/BD is another internal
C           routine of MA57. ICNTL(14) has been added to control
C           density of rows regarded as dense by the MC50 and MA27
C           orderings.
C 24/05/04  Statment functions in MA57U/UD replaced by in-line code.

C 12th July 2004 Version 1.0.0. Version numbering added.

C 20/07/04  Several changes incorporated for HSL 2004 code.
C           Removed unused INT,ABS from MA57U/UD
C           INFO(32), INFO(33), and INFO(34) added
C           INFO(32): no. of zeros in the triangle of the factors
C           INFO(33): no. of zeros in the rectangle of the factors
C           INFO(34): no. of zero columns in rectangle of the factors
C           Static pivoting available (controlled by CNTL(4), CNTL(5))
C           Scaling using symmetrized MC64 (ICNTL(15))
C           Links to METIS_NODEND ordering


C 31st July 2004 Version 2.0.0 established at HSL 2004 release.

C 1st Sept  2004 Version 2.1.0. Default changed to static pivoting off.
C 10th Sept 2004 Version 2.2.0. Defaults for ICNTL(6), ICNTL(9) and
C           CNTL(5) changed. Scaling factors (optionally) printed.
C  4th Nov  2004 Version 2.2.1. Change to assembly of reals in MA57O/OD
C           leading to more efficient code at suggestion of Stephane
C           Pralet.
C 13th Dec  2004 Version 2.3.0. Several minor changes after field
C           testing.
C           Scale factors (RINFO(16) and RINFO(17) set to 1
C           if scaling not used.
C           Option to handle dense columns invoked for METIS ordering.
C           Value of SCHNAB(1) set to 1. to allow Schnabel-Eskow to
C           work on matrix with a rows of zeros.
C           Some diagnostic printing and STOP statements removed from
C           MC50.

C 2nd March 2005  Version 3.0.0.  A new option has been added for
C           ordering the matrix.  If ICNTL(6) is equal to 5 then the
C           ordering chosen depends on the matrix characteristics.
C           At the moment the choices are MC50 or METIS.
C           INFO(36) is set to ordering used.
C           A minor chnage has been made to the pivot control to reduce
C           the amount of researching on failed pivots (resetting of
C           KR). FD05 dependence changed to FD15.
C 15th June 2005  Version 3.0.1.  Setting of ALENB in MA57B/BD moved
C           before first error exit to avoid undefined variable
C           if error invoked.  INFO(1) initialized to zero in call to
C           MA57C/CD.
C 1 December 2006. Version 3.0.2. Comments adjusted to meet the
C           72-character limit.

C 3 August 2007  Version 3.1.0.
C           The new version of MC47 (that incorporates an updated
C           version of MC50 is used).

C 19 September 2007  Version 3.2.0
C           New option added (ICNTL(16)) that allows the removal of
C           blocks of small entries to the end of the factorization.
C           Is particularly powerful when matrix is severely rank
C           deficient.  Numerous other mainly cosmetic changes that
C           don't affect interface.

      SUBROUTINE MA57ID(CNTL, ICNTL)
C****************************************************************
      DOUBLE PRECISION    CNTL(5)
      INTEGER             ICNTL(20)
      INTEGER I
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C===============================================
C===============================================
      CNTL(1)   = 0.01D0
      CNTL(2)   = 1.0D-20
      CNTL(3)   = 0.5D0
      CNTL(4) = ZERO
      CNTL(5) = ZERO
      ICNTL(1)  = 6
      ICNTL(2)  = 6
      ICNTL(3)  = 6
      ICNTL(4)  = -1
      ICNTL(5)  = 2
      ICNTL(6)  = 5
      ICNTL(7)  = 1
      ICNTL(8)  = 0
      ICNTL(9)  = 10
      ICNTL(10) = 0
      ICNTL(11) = 16
      ICNTL(12) = 16
      ICNTL(13) = 10
      ICNTL(14) = 100
      ICNTL(15) = 1
      ICNTL(16) = 0
      DO 110 I=17,20
        ICNTL(I) = 0
  110 CONTINUE
      RETURN
      END
      SUBROUTINE MA57AD(N,NE,IRN,JCN,LKEEP,KEEP,IWORK,ICNTL,INFO,RINFO)
      INTEGER N,NE,IRN(NE),JCN(NE),IWORK(5*N),LKEEP,KEEP(LKEEP),
     *        ICNTL(20),INFO(40)
      DOUBLE PRECISION RINFO(20)
C**** Still to be updated
      INTRINSIC MIN
      EXTERNAL MA57GD,MC47ID,MC47BD,MA57VD,MA57HD,MA57JD,MA57KD,
     *         MA57LD,MA57MD,MA57ND
      INTEGER I,IL,IN,IPE,IRNPRM,COUNT,FILS,FRERE,HOLD,IFCT,INVP,IPS,
     +        IW,IWFR,K,LDIAG,LP,LW,LROW,MAP,EXPNE,
     +        MP,NCMPA,NEMIN,NODE,NST,NSTEPS,NV,PERM,
     +        IW1,IW2,IW3,IW4,IW5,NSTK,ND,NELIM,NZE,ALENB,
     +        J,JJ,J1,J2,SIZE22,OXO
      INTEGER METOPT(8),METFTN,ICNTL6,INF47(10),ICNT47(10)
      DOUBLE PRECISION ZERO,THRESH,AVNUM,MC47FI,RINF47(10)
      PARAMETER (ZERO=0.0D0)
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      DO 10 I = 1,40
        INFO(I) = 0
   10 CONTINUE
      DO 11 I = 1,20
        RINFO(I) = ZERO
   11 CONTINUE
      IF (N.LT.1)  GO TO 20
      IF (NE.LT.0) GO TO 30
      IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
      IF (ICNTL(6).EQ.1) THEN
        DO 12 I = 1,N
          IWORK(I) = 0
   12   CONTINUE
        DO 14 I=1,N
          K = KEEP(I)
          IF (K.LE.0 .OR. K.GT.N) GO TO 80
          IF (IWORK(K).NE.0) GO TO 80
          IWORK(K) = I
   14   CONTINUE
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE(MP,99980) N,NE,(ICNTL(I),I=1,7),ICNTL(12),ICNTL(15)
99980 FORMAT (//'Entering analysis phase (MA57AD) with ...'/
     1      'N         Order of matrix                     =',I12/
     2      'NE        Number of entries                   =',I12/
     6      'ICNTL(1)  Stream for errors                   =',I12/
     7      ' --- (2)  Stream for warnings                 =',I12/
     8      ' --- (3)  Stream for monitoring               =',I12/
     9      ' --- (4)  Stream for statistics               =',I12/
     1      ' --- (5)  Level of diagnostic printing        =',I12/
     2      ' --- (6)  Flag for input pivot order          =',I12/
     2      ' --- (7)  Numerical pivoting control (st est) =',I12/
     2      ' --- (12) Node amalgamation parameter         =',I12/
     2      ' --- (15) Scaling control (storage estimate)  =',I12)
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(3(I6,A,2I8,A)))') ' Matrix entries:',
     +        (I,': (',IRN(I),JCN(I),')',I=1,K)
        IF (K.LT.NE) WRITE (MP,'(A)') '     . . .'
        IF (ICNTL(6).EQ.1) THEN
          K = MIN(10,N)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(A,10I6:/(7X,10I6))') ' KEEP =', (KEEP(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(7X,A)') '     . . .'
        END IF
      END IF
      IW1 = 1
      IW2 = IW1 + N
      IW3 = IW2 + N
      IW4 = IW3 + N
      IW5 = IW4 + N
      FILS  = IW1
      FRERE = IW2
      ND    = IW3
      NELIM = IW4
      NV    = IW5
      PERM = 1
      NSTEPS = PERM + N
      EXPNE  = NSTEPS + 1
      HOLD   = EXPNE + 1
      LROW = HOLD + 40
      NODE = LROW + N
      NSTK = NODE + N
      MAP  = NSTK + N
      IRNPRM = MAP + MAX(N,NE)
      INVP  = NODE
      IW    = NODE
      IPE   = LROW
      IFCT  = MAP
      IPS   = MAP
      COUNT = NSTK
      KEEP(HOLD) = 0
      ICNTL6 = ICNTL(6)
      IF (ICNTL(6).GT.5) ICNTL6 = 5
      IF (ICNTL6.EQ.4 .OR. ICNTL6.EQ.5) THEN
        METFTN    = 1
        METOPT(1) = 0
        KEEP(IPE)   = 1
        KEEP(IPE+1) = 2
        KEEP(IFCT)  = 1
        CALL METIS_NODEND(1,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                    KEEP(NSTK),KEEP(PERM))
        IF (KEEP(PERM).EQ.-1) THEN
          IF (ICNTL6 .EQ. 4) GO TO 90
          ICNTL6 = 2
        ENDIF
      ENDIF
      IF (ICNTL6.NE.1) THEN
        CALL MC47ID(ICNT47)
        IF (ICNTL6 .NE. 3) THEN
          CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +                KEEP(IW),IWFR,ICNTL,INFO)
          IF (ICNTL6.EQ.5) THEN
            IF (ICNTL(7).EQ.2) THEN
              AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
              IF (N.GE.50000) THEN
                ICNTL6 = 4
                GO TO 97
              ENDIF
              IF (N.LE.30000) THEN
                ICNTL6 = 2
                IF (AVNUM.GT.100.0) ICNTL6 = 4
                GO TO 97
              ENDIF
              IF (N.GT.30000 .AND. N.LT.50000) THEN
                IF (AVNUM.GT.46.0) THEN
                  ICNTL6 = 4
                ELSE
                  ICNTL6 = 2
                ENDIF
                GO TO 97
              ENDIF
            ELSE
              AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
              OXO = 0
              J2 = IWFR - 1
              SIZE22 = 0
              DO 100 J = N,1,-1
                J1 = KEEP(IPE+J-1)
                DO  99 JJ = J1,J2
                  IF (KEEP(IFCT+JJ-1).GT.J) GO TO 101
   99           CONTINUE
                SIZE22 = SIZE22 + 1
                J2 = J1-1
  100         CONTINUE
  101         IF (SIZE22 .GT. 0) THEN
                DO 98 I = 1,NE
                  IF (IRN(I) .LE. N-SIZE22
     *          .AND. JCN(I) .LE. N-SIZE22) THEN
                      AVNUM = FLOAT(IWFR+N-SIZE22-1)/FLOAT(N)
                      GO TO 96
                  ENDIF
   98           CONTINUE
                OXO = 1
                AVNUM = FLOAT(IWFR-1)/FLOAT(N)
              ENDIF
   96         IF (N .GE. 100000) THEN
                IF (AVNUM.GT.5.42) THEN
                  ICNTL6 = 4
                ELSE
                  ICNTL6 = 2
                ENDIF
                GO TO 97
              ENDIF
              IF (OXO.EQ.1) THEN
                IF (FLOAT(N-SIZE22)/FLOAT(SIZE22) .GT .1.8D0) THEN
                  ICNTL6 = 2
                ELSE
                  ICNTL6 = 4
                ENDIF
                GO TO 97
              ENDIF
              LW = LKEEP-IFCT+1
              CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
     +                    KEEP(IFCT),IWORK(NV),
     +                    KEEP(INVP),KEEP(PERM),IWORK(IW1),
     +                    IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +                    ICNT47,INF47,RINF47)
              INFO(13) = INF47(2)
              ICNTL6 = 2
              NEMIN    = ICNTL(12)
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
              NST = KEEP(NSTEPS)
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))
              KEEP(EXPNE) = IWORK(IW5)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)
              IF (FLOAT(INFO(5))/FLOAT(NE) .LT. 10.0) THEN
                GO TO 93
              ELSE
                MC47FI = FLOAT(INFO(5))/FLOAT(NE)
        CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +              KEEP(IW),IWFR,ICNTL,INFO)
                KEEP(IPE+N) = IWFR
                METFTN    = 1
                METOPT(1) = 0
                IF (N.LT.50) GO TO 92
                DO 91 I = 1,N
                  IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
                    METOPT(1) = 1
                    METOPT(2) = 3
                    METOPT(3) = 1
                    METOPT(4) = 2
                    METOPT(5) = 0
                    METOPT(6) = 1
                    METOPT(7) = 200
                    METOPT(8) = 1
                    GO TO 92
                  ENDIF
   91           CONTINUE
   92     CALL METIS_NODEND(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                      KEEP(NSTK),KEEP(PERM))
        CALL MA57JD(N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
     +              KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
                LW = 2*NE
        CALL MA57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
     +              KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
                INFO(13) = NCMPA
                NEMIN = ICNTL(12)
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
                NST = KEEP(NSTEPS)
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))
                KEEP(EXPNE) = IWORK(IW5)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)
                IF (FLOAT(INFO(5))/FLOAT(NE).LT.MC47FI) THEN
                  ICNTL6 = 4
                  GO TO 93
                ELSE
                  ICNTL6=2
        CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +              KEEP(IW),IWFR,ICNTL,INFO)
                  GO TO 97
                ENDIF
              ENDIF
            ENDIF
          ENDIF
   97     IF (ICNTL6.EQ.4) THEN
            KEEP(IPE+N) = IWFR
            METFTN    = 1
            METOPT(1) = 0
            IF (N.LT.50) GO TO 103
            DO 102 I = 1,N
              IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
                METOPT(1) = 1
                METOPT(2) = 3
                METOPT(3) = 1
                METOPT(4) = 2
                METOPT(5) = 0
                METOPT(6) = 1
                METOPT(7) = 200
                METOPT(8) = 1
                GO TO 103
              ENDIF
  102       CONTINUE
  103     CALL METIS_NODEND(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                      KEEP(NSTK),KEEP(PERM))
            GO TO 111
          ENDIF
          LW = LKEEP-IFCT+1
          IF (ICNTL6 .EQ. 0) ICNT47(4) = -1
          CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
     +                KEEP(IFCT),IWORK(NV),
     +                KEEP(INVP),KEEP(PERM),IWORK(IW1),
     +                IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +                ICNT47,INF47,RINF47)
          INFO(13) = INF47(2)
        ELSE
          LW = LKEEP-IFCT+1
        CALL MA57VD(N,NE,IRN,JCN,KEEP(IFCT),LW,KEEP(IPE),IWORK(IW1),
     *              IWORK(IW2),IWFR,ICNTL,INFO)
          THRESH = FLOAT(ICNTL(14))/100.0
        CALL MA57HD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,IWORK(NV),
     *              IWORK(IW1),IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +              2139062143,INFO(13),THRESH)
          DO 110 I = 1,N
            IF (IWORK(NV+I-1).NE.0) GO TO 110
            IN = I
  105       IL = IN
            IN = - KEEP(IPE+IL-1)
            IF (IWORK(NV+IN-1).EQ.0) GO TO 105
            KEEP(IPE+I-1) = -IN
  110     CONTINUE
        ENDIF
      ENDIF
  111 IF (ICNTL6.EQ.1 .OR. ICNTL6.EQ.4) THEN
        CALL MA57JD(N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
     +              KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
        LW = 2*NE
        CALL MA57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
     +              KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
        INFO(13) = NCMPA
      END IF
      NEMIN = ICNTL(12)
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
      NST = KEEP(NSTEPS)
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))
      KEEP(EXPNE) = IWORK(IW5)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)
   93 INFO(36) = ICNTL6
      ALENB    = 1
      IF (ICNTL(7).EQ.4) ALENB = ALENB + N + 5
      IF (ICNTL(15).EQ.1) ALENB = ALENB + N
      INFO(9)  = MAX(INFO(9)+ALENB,ALENB+KEEP(EXPNE)+1)
      INFO(11) = MAX(INFO(11)+ALENB,ALENB+KEEP(EXPNE)+1)
      INFO(10) = MAX(INFO(10),KEEP(EXPNE)+N+5)
      INFO(12) = MAX(INFO(12),KEEP(EXPNE)+N+5)
      IF (ICNTL(15).EQ.1) THEN
        INFO(9) = MAX(INFO(9),ALENB+3*KEEP(EXPNE)+3*N)
        INFO(11) = MAX(INFO(11),ALENB+3*KEEP(EXPNE)+3*N)
        INFO(10) = MAX(INFO(10),3*KEEP(EXPNE)+5*N+1)
        INFO(12) = MAX(INFO(12),3*KEEP(EXPNE)+5*N+1)
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        NZE = KEEP(EXPNE)
        WRITE (MP,99999) INFO(1),NZE,
     *                  (INFO(I),I=3,13),INFO(36),(RINFO(I),I=1,2)
99999 FORMAT (/'Leaving analysis phase (MA57AD) with ...'/
     1    'INFO(1)  Error indicator                      =',I12/
     2    'Number of entries in matrix with diagonal     =',I12/
     2    'INFO(3)  Number of out-of-range indices       =',I12/
     2    'INFO(4)  Number of off-diagonal duplicates    =',I12/
     2    'INFO(5)  Forecast real storage for factors    =',I12/
     3    '----(6)  Forecast integer storage for factors =',I12/
     3    '----(7)  Forecast maximum front size          =',I12/
     4    '----(8)  Number of nodes in assembly tree     =',I12/
     5    '----(9)  Size of FACT without compress        =',I12/
     6    '----(10) Size of IFACT without compress       =',I12/
     5    '----(11) Size of FACT with compress           =',I12/
     5    '----(12) Size of IFACT with compress          =',I12/
     5    '----(13) Number of compresses                 =',I12/
     5    '----(36) Ordering strategy used by code       =',I12/
     9    'RINFO(1) Forecast additions for assembly      =',1P,D12.5/
     9    'RINFO(2) Forecast ops for elimination         =',1P,D12.5)
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')  'Permutation array:',
     +                  (KEEP(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +        'Number of entries in rows of permuted matrix:',
     +        (KEEP(LROW+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NZE)
        IF (LDIAG.GE.4) K = NZE
        WRITE (MP,'(/A/(5I12))')
     *        'Column indices of permuted matrix:',
     *                           (KEEP(IRNPRM+I-1),I=1,K)
        IF (K.LT.NZE) WRITE (MP,'(16X,A)') '     . . .'
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')
     +    'Tree nodes at which variables eliminated:',
     +    (KEEP(NODE+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(5I12))') 'Map array:',
     *                               (KEEP(I),I=MAP,MAP+K-1)
        IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
      END IF
      RETURN
   20 INFO(1) = -1
      INFO(2) = N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'N has value ',INFO(2)
      RETURN
   30 INFO(1) = -2
      INFO(2) = NE
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'NE has value',INFO(2)
       RETURN
   40 INFO(1) = -15
      INFO(2) = LKEEP
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'LKEEP has value    ',INFO(2),
     +    'Should be at least ',5*N+NE+MAX(N,NE)+42
       RETURN
   80 INFO(1) = -9
      INFO(2) = I
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A/A,I10,A)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'Invalid permutation supplied in KEEP',
     +    'Component',INFO(2),' is faulty'
      RETURN
   90 INFO(1) = -18
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'MeTiS ordering requested but MeTiS not linked'
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA57BD(N, NE, A, FACT, LFACT, IFACT, LIFACT,
     * LKEEP, KEEP, PPOS, ICNTL, CNTL, INFO, RINFO)
      INTEGER N,NE,LFACT,LIFACT,LKEEP
      DOUBLE PRECISION A(NE),FACT(LFACT)
      DOUBLE PRECISION RINFO(20)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(20), IFACT(LIFACT)
      INTEGER   INFO(40), KEEP(LKEEP), PPOS(N)
      INTEGER EXPNE,HOLD,I,IRNPRM,K,LDIAG,LLFACT,LP,LROW,MAP,MM1,MM2,MP
      INTEGER J,JJ,KK,ISCALE,NUM,NE64,IDUP,IMAT,IPT,JLOOP,JNEW,NN,ISING
      INTEGER NSTEPS,NODE,NSTK,PERM,INEW,ALENB,BIGA
      DOUBLE PRECISION ONE,ZERO,RINF,FD15AD,FCT,SMAX,SMIN,REPS
      PARAMETER (ONE = 1.0D0, ZERO=0.0D0)
      INTRINSIC MIN
      EXTERNAL MA57OD,MA57UD,FD15AD,MC34AD,MC64WD
      RINF = FD15AD('H')
      REPS = FD15AD('E')
      LP     = ICNTL(1)
      MP     = ICNTL(3)
      LDIAG  = ICNTL(5)
C??
      IF (N.LE.0)  GO TO 25
      IF (NE.LT.0) GO TO 30
      IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
      IF (ICNTL(7).LT.1 .OR. ICNTL(7).GT.4) GO TO 35
      NSTEPS = KEEP(N+1)
      EXPNE  = KEEP(N+2)
      PERM = 1
      HOLD = PERM + N + 2
      LROW = HOLD + 40
      NODE = LROW + N
      NSTK = NODE + N
      MAP  = NSTK + N
      IRNPRM = MAP + MAX(NE,N)
      BIGA = LFACT
      LLFACT = LFACT - 1
      IF (ICNTL(15).EQ.1) THEN
        ISCALE = LLFACT - N + 1
        LLFACT = ISCALE - 1
      ENDIF
      IF (ICNTL(7).EQ.4) THEN
        LLFACT = LLFACT - N - 5
        MM1 = LLFACT+6
        MM2 = LLFACT+1
      ELSE
        MM1 = 1
        MM2 = 1
      ENDIF
      ALENB = 1
      IF (ICNTL(7).EQ.4)  ALENB = ALENB + N + 5
      IF (ICNTL(15).EQ.1) ALENB = ALENB + N
      IF (LLFACT.LT.EXPNE+1)   GO TO 85
      IF (LIFACT.LT.EXPNE+N+5)  GO TO 95
      IF (ICNTL(15).EQ.1)  THEN
        IF (LFACT .LT. ALENB + 3*EXPNE  + 3*N) GO TO 85
        IF (LIFACT .LT. 3*EXPNE + 5*N + 1) GO TO 95
      ENDIF
C*****************************
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999)
99999 FORMAT (//'Entering factorization phase (MA57BD) with ...')
        IF (KEEP(HOLD).GT.0) WRITE (MP,99998)
99998 FORMAT ('Re-entry call after call to MA57ED')
        WRITE (MP,99997) N,NE,EXPNE,(ICNTL(I),I=1,5),ICNTL(7),ICNTL(8),
     +         ICNTL(11),ICNTL(15),LFACT, LIFACT, NSTEPS,
     +         CNTL(1), CNTL(2), CNTL(4), CNTL(5)
99997 FORMAT ('N       Order of input matrix               =',I12/
     2        'NE      Entries in input matrix             =',I12/
     2        '        Entries in input matrix (inc diags) =',I12/
     6        'ICNTL(1)  Stream for errors                 =',I12/
     7        ' --- (2)  Stream for warnings               =',I12/
     8        ' --- (3)  Stream for monitoring             =',I12/
     9        ' --- (4)  Stream for statistics             =',I12/
     1        ' --- (5)  Level of diagnostic printing      =',I12/
     1        ' --- (7)  Numerical pivoting control        =',I12/
     1        ' --- (8)  Restart or discard factors        =',I12/
     1        ' --- (11) Block size for Level 3 BLAS       =',I12/
     1        ' --- (15) Scaling control (1 on)            =',I12/
     4        'LFACT   Size of real working space          =',I12/
     5        'LIFACT  Size of integer working space       =',I12/
     7        '        Number nodes in assembly tree       =',I12/
     9        'CNTL(1) Value of threshold parameter        =',D12.5/
     9        'CNTL(2) Threshold for zero pivot            =',D12.5/
     9        'CNTL(4) Control for value of static pivots  =',D12.5/
     9        'CNTL(5) Control for number delayed pivots   =',D12.5)
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        IF (NE.GT.0) THEN
          WRITE (MP,'(/A/(3(I6,A,1P,D16.8,A)))') 'Matrix entries:',
     +     (I,': (',A(I),')',I=1,K)
          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .'
        END IF
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')  'Permutation array:',
     +                    (KEEP(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +          'Number of entries in rows of permuted matrix:',
     +          (KEEP(LROW+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +    'Tree nodes at which variables eliminated:',
     +    (KEEP(NODE+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NSTEPS)
        IF (LDIAG.GE.4) K = NSTEPS
        IF (K.GT.0) WRITE (MP,'(/A/(5I12))')
     +     'Number of assemblies at each tree node:',
     +     (KEEP(NSTK+I-1),I=1,K)
        IF (K.LT.NSTEPS) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(5I12))') 'Map array:',
     *                               (KEEP(I),I=MAP,MAP+K-1)
        IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,EXPNE)
        IF (LDIAG.GE.4) K = EXPNE
        WRITE (MP,'(/A/(5I12))')
     *          'Column indices of permuted matrix:',
     *                             (KEEP(IRNPRM+I-1),I=1,K)
        IF (K.LT.EXPNE) WRITE (MP,'(16X,A)') '     . . .'
      ENDIF
      IF (KEEP(HOLD) .GT. 0) GO TO 22
C***************************************************
C***************************************************
C?? For the moment to handle missing diagonals
      DO 19 K = 1,EXPNE
        FACT(LLFACT-EXPNE+K) = ZERO
   19 CONTINUE
      FACT(BIGA) = ZERO
      DO 20 K = 1,NE
        FACT(BIGA) = MAX(FACT(BIGA),ABS(A(K)))
        FACT(KEEP(MAP+K-1)+LLFACT-EXPNE) = A(K)
   20 CONTINUE
      RINFO(18) = FACT(BIGA)
      DO 21 K = 1,EXPNE
        IFACT(LIFACT-EXPNE+K) = KEEP(IRNPRM+K-1)
   21 CONTINUE
      DO 23 I = 1,N
        PPOS(KEEP(PERM+I-1)) = I
   23 CONTINUE
      IF (ICNTL(15).EQ.1) THEN
        IPT = 1
        IDUP = IPT+N+1
        IMAT = IDUP+N
        ISING = IMAT + MAX(NE,EXPNE)
        DO 4444 I = 1,N
          IFACT(IDUP+I-1) = 0
 4444   CONTINUE
C9999   CONTINUE
        IFACT(IPT) = 1
        KK = 1
        K = 1
        DO 3333 J = 1,N
          DO 2222 JJ = 1,KEEP(LROW+J-1)
            I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
            IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
              FACT(IFACT(IDUP+I-1)) =
     &          FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
            ELSE
              IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
                IFACT(IDUP+I-1) = KK
                FACT(KK) = FACT(LLFACT-EXPNE+K)
                IFACT(IMAT-1+KK) = I
                KK = KK+1
              ENDIF
            ENDIF
            K = K + 1
 2222     CONTINUE
          IFACT(IPT+J) = KK
 3333   CONTINUE
        CALL MC34AD(N,IFACT(IMAT),IFACT(IPT),.TRUE.,FACT,KEEP(PERM))
        NE64 = IFACT(IPT+N)-1
        DO 75 J = 1,N
          FCT = ZERO
          DO 60 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            FACT(K) = ABS(FACT(K))
            IF (FACT(K).GT.FCT) FCT = FACT(K)
   60     CONTINUE
          FACT(NE64+2*N+J) = FCT
          IF (FCT.NE.ZERO) THEN
            FCT = LOG(FCT)
          ELSE
            FCT = RINF/N
          ENDIF
          DO 70 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
CCC
              FACT(K) = FCT - LOG(FACT(K))
   70     CONTINUE
   75   CONTINUE
        CALL MC64WD(N,NE64,IFACT(IPT),IFACT(IMAT),FACT,KEEP(PERM),NUM,
     &      IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
     &      IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
     &      FACT(NE64+1),FACT(NE64+N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
              FACT(NE64+N+J) = FACT(NE64+N+J) - LOG(FACT(NE64+2*N+J))
CCC
   80     CONTINUE
          DO 5555 I=1,N
            FACT(ISCALE+PPOS(I)-1) =
     &        SQRT(EXP(FACT(NE64+I)+FACT(NE64+N+I)))
 5555     CONTINUE
        ELSE
        K = 0
        DO 3501 I = 1,N
          IF (KEEP(PERM+I-1).LT.0) THEN
            PPOS(I) = -PPOS(I)
            IFACT(ISING+I-1) = 0
          ELSE
            K = K + 1
            IFACT(ISING+I-1) = K
          ENDIF
 3501   CONTINUE
        DO 3502 I = 1,N
          KEEP(PERM+ABS(PPOS(I))-1) = I
 3502   CONTINUE
        DO 3503 I = 1,N
          IFACT(IDUP+I-1) = 0
 3503   CONTINUE
        IFACT(IPT) = 1
        KK = 1
        K = 1
        JNEW = 0
        NN = N
        DO 3505 J = 1,N
          IF (PPOS(J).LT.0) THEN
            NN = NN - 1
            K = K + KEEP(LROW+J-1)
            GO TO 3505
          ENDIF
          JNEW = JNEW + 1
          DO 3504 JJ = 1,KEEP(LROW+J-1)
            I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
            IF (PPOS(I).GT.0) THEN
              IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
                FACT(IFACT(IDUP+I-1)) =
     &            FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
              ELSE
                IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
                  IFACT(IDUP+I-1) = KK
                  FACT(KK) = FACT(LLFACT-EXPNE+K)
                  IFACT(IMAT-1+KK) = IFACT(ISING+I-1)
                  KK = KK+1
                ENDIF
              ENDIF
            ENDIF
            K = K + 1
 3504     CONTINUE
          IFACT(IPT+JNEW) = KK
 3505   CONTINUE
      NE64 = IFACT(IPT+NN)-1
        CALL MC34AD(NN,IFACT(IMAT),IFACT(IPT),.TRUE.,FACT,KEEP(PERM))
        NE64 = IFACT(IPT+NN)-1
        DO 3508 J = 1,NN
          FCT = ZERO
          DO 3506 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            FACT(K) = ABS(FACT(K))
            IF (FACT(K).GT.FCT) FCT = FACT(K)
 3506     CONTINUE
          FACT(NE64+2*N+J) = FCT
CCC
            FCT = LOG(FCT)
          DO 3507 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
              FACT(K) = FCT - LOG(FACT(K))
 3507     CONTINUE
 3508   CONTINUE
        CALL MC64WD(NN,NE64,IFACT(IPT),IFACT(IMAT),FACT,KEEP(PERM),NUM,
     &      IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
     &      IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
     &      FACT(NE64+1),FACT(NE64+N+1))
        DO 3509 J = 1,NN
CCC
              FACT(NE64+N+J) = FACT(NE64+N+J) - LOG(FACT(NE64+2*N+J))
              FACT(NE64+N+J) = ZERO
 3509     CONTINUE
          K=0
          DO 3510 I=1,N
            IF (PPOS(I).LT.0) THEN
              K = K + 1
              FACT(ISCALE-PPOS(I)-1) = ZERO
            ELSE
              FACT(ISCALE+PPOS(I)-1) =
     &          SQRT(EXP(FACT(NE64+I-K)+FACT(NE64+N+I-K)))
            ENDIF
 3510     CONTINUE
          DO 3516 I = 1,N
            KEEP(PERM+ABS(PPOS(I))-1) = I
 3516     CONTINUE
          K = 1
          DO 3514 JJ = 1,N
            J = PPOS(JJ)
            IF (J.GT.0) THEN
              DO 3511 JLOOP = 1,KEEP(LROW+JJ-1)
                I = IFACT(LIFACT-EXPNE+K)
                INEW = KEEP(PERM+I-1)
                IF (PPOS(INEW).LT.0)
     &            FACT(ISCALE+I-1) = MAX(FACT(ISCALE+I-1),
     &                 ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+J-1))
                K = K + 1
 3511         CONTINUE
            ELSE
              DO 3512 JLOOP = 1,KEEP(LROW+JJ-1)
                I = IFACT(LIFACT-EXPNE+K)
                INEW = KEEP(PERM+I-1)
                IF (I .NE. -J)  THEN
                FACT(ISCALE-J-1) =
     &              MAX(FACT(ISCALE-J-1),
     &              ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+I-1))
                ENDIF
                K = K + 1
 3512         CONTINUE
            ENDIF
 3514     CONTINUE
          DO 3513 I = 1,N
            INEW = KEEP(PERM+I-1)
            IF (PPOS(INEW) .LT. 0) THEN
              PPOS(INEW) = - PPOS(INEW)
              IF (FACT(ISCALE+I-1) .EQ. ZERO) THEN
                FACT(ISCALE+I-1) = ONE
              ELSE
                FACT(ISCALE+I-1) = ONE/FACT(ISCALE+I-1)
              ENDIF
            ENDIF
 3513     CONTINUE
        ENDIF
C8888     CONTINUE
          SMAX = FACT(ISCALE)
          SMIN = FACT(ISCALE)
          DO 5566 I = 1,N
            SMAX = MAX(SMAX,FACT(ISCALE+I-1))
            SMIN = MIN(SMIN,FACT(ISCALE+I-1))
 5566     CONTINUE
          RINFO(16) = SMIN
          RINFO(17) = SMAX
          K = 1
          FACT(BIGA) = ZERO
          DO 6666 JJ = 1,N
            J = PPOS(JJ)
            DO 7777 JLOOP = 1,KEEP(LROW+JJ-1)
              I = IFACT(LIFACT-EXPNE+K)
              FACT(LLFACT-EXPNE+K) =
     &          FACT(ISCALE+I-1)*FACT(LLFACT-EXPNE+K)*FACT(ISCALE+J-1)
              FACT(BIGA) = MAX(FACT(BIGA), ABS(FACT(LLFACT-EXPNE+K)))
              K = K + 1
 7777       CONTINUE
 6666     CONTINUE
      ELSE
        RINFO(16) = ONE
        RINFO(17) = ONE
      ENDIF
C**********************************
C**********************************
   22 CALL MA57OD(N, EXPNE, FACT, LLFACT, IFACT, LIFACT, KEEP(LROW),
     *            PPOS,
     *            NSTEPS, KEEP(NSTK), KEEP(NODE), FACT(MM1),
     *            FACT(MM2),
     *            KEEP(PERM),
     *            CNTL, ICNTL,
     *            INFO, RINFO, KEEP(HOLD), FACT(BIGA))
      IF (INFO(1).EQ.10 .OR. INFO(1).EQ.11) THEN
        IF (LDIAG.GT.2 .AND. MP.GE.0)  THEN
          IF (INFO(1).EQ.10) WRITE (MP,99982) INFO(1)
99982 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'Factorization suspended because of lack of real space'/
     1  'INFO (1) = ',I3)
          IF (INFO(1).EQ.11) WRITE (MP,99983) INFO(1)
99983 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'Factorization suspended because of lack of integer space'/
     1  'INFO (1) = ',I3)
        ENDIF
        RETURN
      ENDIF
      DO 24 I = 1,N
        KEEP(PERM+PPOS(I)-1) = I
   24 CONTINUE
        INFO(17) = ALENB + INFO(17)
        INFO(19) = ALENB + INFO(19)
      IF (ICNTL(15).EQ.1) THEN
        INFO(17) = MAX(INFO(17),ALENB + 3*EXPNE+3*N)
        INFO(19) = MAX(INFO(19),ALENB + 3*EXPNE+3*N)
        INFO(18) = MAX(INFO(18),3*EXPNE+5*N+1)
        INFO(20) = MAX(INFO(18),3*EXPNE+5*N+1)
      ENDIF
      IF (INFO(1).LT.0) RETURN
      GO TO 100
C************************
C************************
   25 INFO(1) = -1
      INFO(2) =  N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'N has value ',INFO(2)
      RETURN
   30 INFO(1) = -2
      INFO(2) = NE
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'NE has value',INFO(2)
      RETURN
   40 INFO(1) = -15
      INFO(2) = LKEEP
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'LKEEP has value    ',INFO(2),
     +    'Should be at least ',5*N+NE+MAX(N,NE)+42
      RETURN
   35 INFO(1) = -10
      INFO(2) = ICNTL(7)
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'ICNTL(7) has value',ICNTL(7)
      RETURN
   85 INFO(1) = -3
      INFO(2) = LFACT
      INFO(17) = ALENB + EXPNE + 1
      IF (ICNTL(15).EQ.1) INFO(17) = ALENB + 3*EXPNE + 3*N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'Insufficient real space in FACT, LFACT = ',INFO(2)
      RETURN
   95 INFO(1) = -4
      INFO(2) = LIFACT
      INFO(18) = EXPNE+N+5
      IF (ICNTL(15).EQ.1) INFO(18) = 3*EXPNE + 5*N + 1
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'Insufficient integer space in IFACT, LIFACT = ',INFO(2)
      RETURN
C****************
C****************
 100  IF (LDIAG.LE.2 .OR. MP.LT.0) RETURN
      WRITE (MP,99980) INFO(1), INFO(2),
     *    (INFO(I),I=14,25),INFO(28),INFO(29)
      WRITE (MP,99984) (INFO(I),I=31,35),RINFO(3), RINFO(4),
     *                 RINFO(5), RINFO(18)
99980 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'INFO (1)                                      =',I12/
     2  ' --- (2)                                      =',I12/
     3  ' --- (14) Number of entries in factors        =',I12/
     4  ' --- (15) Real storage for factors            =',I12/
     5  ' --- (16) Integer storage for factors         =',I12/
     6  ' --- (17) Min LFACT with compresses           =',I12/
     7  ' --- (18) Min LIFACT with compresses          =',I12/
     8  ' --- (19) Min LFACT without compresses        =',I12/
     9  ' --- (20) Min LIFACT without compresses       =',I12/
     *  ' --- (21) Order of largest frontal matrix     =',I12/
     1  ' --- (22) Number of 2x2 pivots                =',I12/
     2  ' --- (23) Number of delayed pivots            =',I12/
     3  ' --- (24) Number of negative eigenvalues      =',I12/
     4  ' --- (25) Rank of factorization               =',I12/
     5  ' --- (28) Number compresses on real data      =',I12/
     6  ' --- (29) Number compresses on integer data   =',I12)
      IF (ICNTL(15).EQ.1) WRITE (MP,99985) RINFO(16),RINFO(17)
99985 FORMAT (
     1  'RINFO(16) Minimum value of scaling factor     =  ',1PD10.3/
     2  '-----(17) Maximum value of scaling factor     =  ',1PD10.3)
99984 FORMAT (
     7  ' --- (31) Number of block pivots in factors   =',I12/
     7  ' --- (32) Number of zeros factors triangle    =',I12/
     7  ' --- (33) Number of zeros factors rectangle   =',I12/
     7  ' --- (34) Number of zero cols factors rect    =',I12/
     7  ' --- (35) Number of static pivots             =',I12/
     1  'RINFO(3)  Operations during node assembly     =  ',1PD10.3/
     2  '-----(4)  Operations during node elimination  =  ',1PD10.3/
     3  '-----(5)  Extra operations because of BLAS    =  ',1PD10.3/
     3  '-----(18) Largest modulus of entry in matrix  =  ',1PD10.3)
      IF (INFO(27).GT.0) WRITE (MP,99981) INFO(27),RINFO(14),RINFO(15)
99981 FORMAT (/'Matrix modification performed'/
     1  'INFO (27) Step at which matrix first modified =',I12/
     2  'RINFO(14) Maximum value added to diagonal     =  ',1PD10.3/
     2  'RINFO(15) Smallest pivot in modified matrix   =  ',1PD10.3)
      CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
      IF (ICNTL(15).NE.1) RETURN
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
     +                    (FACT(ISCALE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      END
      SUBROUTINE MA57CD(JOB,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,W,
     *                  LW,IW1,ICNTL,INFO)
      INTEGER JOB,N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20),INFO(40)
      INTRINSIC MIN
      EXTERNAL MA57QD,MA57RD,MA57SD,MA57TD,MA57UD,MA57XD,MA57YD
      DOUBLE PRECISION SCALE,ONE
      PARAMETER (ONE = 1.0D0)
      INTEGER I,J,K,LDIAG,LLW,LP,MP,ISCALE
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      INFO(1) = 0
      IF (N.LE.0) THEN
        INFO(1) = -1
        INFO(2) = N
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'N has value',N
        GOTO 500
      ENDIF
      IF (NRHS.LT.1) THEN
        INFO(1) = -16
        INFO(2) = NRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of NRHS =',NRHS,' is less than 1'
        GOTO 500
      ENDIF
      IF (LRHS.LT.N) THEN
        INFO(1) = -11
        INFO(2) = LRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of LRHS =',LRHS,' is less than N=',N
        GOTO 500
      ENDIF
      IF (LW.LT.N*NRHS) THEN
        INFO(1) = -17
        INFO(2) = N*NRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of LW =',LW,' is less than', N*NRHS
        GOTO 500
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999) JOB,N,(ICNTL(I),I=1,5),LFACT,LIFACT,NRHS,
     +         LRHS,LW,ICNTL(13)
99999 FORMAT(/'Entering solution phase (MA57CD) with ...'/
     +    'JOB       Control on coefficient matrix       =',I12/
     +    'N         Order of matrix                     =',I12/
     6    'ICNTL(1)  Stream for errors                   =',I12/
     7    ' --- (2)  Stream for warnings                 =',I12/
     8    ' --- (3)  Stream for monitoring               =',I12/
     9    ' --- (4)  Stream for statistics               =',I12/
     1    ' --- (5)  Level of diagnostic printing        =',I12/
     +    'LFACT     Length of array FACT                =',I12/
     +    'LIFACT    Length of array IFACT               =',I12/
     +    'NRHS      Number of right-hand sides          =',I12/
     +    'LRHS      Leading dimension of RHS array      =',I12/
     +    'LW        Leading dimension of work array     =',I12/
     +    'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12)
        CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
        IF (ICNTL(15).EQ.1) THEN
          ISCALE = LFACT-N
          K = MIN(10,N)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
     +                        (FACT(ISCALE+I-1),I=1,K)
          IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        ENDIF
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        DO 10 J = 1,NRHS
          WRITE(MP,'(/A,I10)') 'Right-hand side',J
          WRITE (MP,'((1P,5D13.3))') (RHS(I,J),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
   10   CONTINUE
      END IF
      LLW = LW/NRHS
      IF (ICNTL(15).EQ.1) THEN
        ISCALE = LFACT-N
        DO 5555 I = 1, N
          SCALE = FACT(ISCALE+I-1)
          IF (JOB.GE.4) SCALE = ONE/FACT(ISCALE+I-1)
          DO 4444 J = 1, NRHS
            RHS(I,J) = SCALE*RHS(I,J)
 4444     CONTINUE
 5555   CONTINUE
      ENDIF
      IF (JOB.LE.2) THEN
        IF (NRHS.EQ.1) THEN
          CALL MA57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ELSE
          CALL MA57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ENDIF
        IF (JOB.EQ.2) GO TO 15
        IF (NRHS.EQ.1) THEN
          CALL MA57YD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ELSE
          CALL MA57RD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ENDIF
      ENDIF
      IF (JOB.EQ.3)
     *  CALL MA57SD(FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *              W,LLW,ICNTL)
      IF (JOB.GE.4)
     *  CALL MA57TD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *              W,LLW,IW1,ICNTL)
   15 IF (ICNTL(15).EQ.1) THEN
        ISCALE = LFACT-N
        DO 6666 I = 1, N
          SCALE = FACT(ISCALE+I-1)
          IF (JOB.EQ.2) SCALE = ONE/FACT(ISCALE+I-1)
          DO 7777 J = 1, NRHS
            RHS(I,J) = SCALE*RHS(I,J)
 7777     CONTINUE
 6666   CONTINUE
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,'(//A)')
     *       'Leaving solution phase (MA57CD) with ...'
        DO 20 J = 1,NRHS
          WRITE(MP,'(/A,I10)') 'Solution       ',J
          WRITE (MP,'(1P,5D13.3)') (RHS(I,J),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
   20   CONTINUE
      ENDIF
  500 RETURN
      END
      SUBROUTINE MA57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,
     +        NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = 1
      IWPOS = 4
      DO 270 IBLK = 1,IFACT(3)
        IW1(IBLK) = IWPOS
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 11 J = 1,NRHS
              W(I,J) = RHS(II,J)
   11       CONTINUE
   10     CONTINUE
          DO 12 J = 1,NRHS
            CALL DTPSV('L','N','U',NROWS,FACT(APOS),W(1,J),1)
   12     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2
          IF (NCOLS.GT.NROWS) CALL DGEMM('N','N',NCOLS-NROWS,NRHS,NROWS,
     +                                  ONE,FACT(APOS),NCOLS-NROWS,
     +                                  W,LW,ONE,W(NROWS+1,1),LW)
          APOS = APOS + NROWS* (NCOLS-NROWS)
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 36 J = 1,NRHS
              RHS(II,J) = W(I,J)
   36       CONTINUE
   35     CONTINUE
        ELSE
        J1 = IWPOS
        J2 = IWPOS + NROWS - 1
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          DO 101 II = 1,NRHS
            W1 = RHS(ABS(IFACT(J1)),II)
            K = APOS
            DO 100 J = J1+1,J2
              IRHS = ABS(IFACT(J))
              RHS(IRHS,II) = RHS(IRHS,II) - FACT(K)*W1
              K = K + 1
  100       CONTINUE
  101     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS
          DO 135 II = 1,NRHS
            K = APOS
            W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)),II)
            DO 133 J = J1,J2
              IRHS = ABS(IFACT(J))
              RHS(IRHS,II) = RHS(IRHS,II) + W1*FACT(K)
              K = K + 1
  133       CONTINUE
  135     CONTINUE
          APOS = K
  136   CONTINUE
      END IF
      IWPOS = IWPOS + NCOLS
  270 CONTINUE
      END
      SUBROUTINE MA57RD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,J1,J2,K,KK,LROW,NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = IFACT(1)
      APOS2 = IFACT(2)
      DO 380 IBLK = IFACT(3),1,-1
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 5 I = NROWS + 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 3 J = 1,NRHS
              W(I,J) = RHS(II,J)
    3       CONTINUE
    5     CONTINUE
          DO 10 IPIV = NROWS,1,-1
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            APOS = APOS - (NROWS+1-IPIV)
            DO 9 J = 1,NRHS
              W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9       CONTINUE
   10     CONTINUE
          JPIV = -1
          DO 20 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              DO 19 J = 1,NRHS
                W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19         CONTINUE
              IF (JPIV.EQ.1) APOS2 = APOS2 - 1
              JPIV = -JPIV
            END IF
   20     CONTINUE
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMM('T','N',NROWS,NRHS,K,ONE,
     +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
     +                           W(NROWS+1,1),LW,ONE,W,LW)
          DO 22 J = 1,NRHS
            CALL DTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE
        ELSE
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1
          JPIV = -1
          DO 210 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            LROW = NROWS + 1 - IPIV
            IF (IRHS.GT.0) THEN
              APOS = APOS - LROW
              DO 65 J = 1,NRHS
                RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65         CONTINUE
            ELSE
              IF (JPIV.EQ.-1) THEN
                IRHS1 = -IFACT(IWPOS+IPIV-2)
                IRHS2 = -IRHS
                APOS = APOS - LROW - LROW - 1
                DO 68 J = 1,NRHS
                  W1 = RHS(IRHS1,J)*FACT(APOS) +
     +                 RHS(IRHS2,J)*FACT(APOS2)
                  RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
     +                           RHS(IRHS2,J)*FACT(APOS+LROW+1)
                  RHS(IRHS1,J) = W1
   68           CONTINUE
                APOS2 = APOS2 - 1
              END IF
              JPIV = -JPIV
            END IF
  210     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2
          KK = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 218 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = KK
              DO 215 J = J1,J2
                W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  215         CONTINUE
              RHS(IRHS,II) = W1
  218       CONTINUE
            KK = K
  220     CONTINUE
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            DO 240 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = APOS + 1
              DO 230 J = J1,J2
                W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  230         CONTINUE
              RHS(IRHS,II) = W1
  240       CONTINUE
            J1 = J1 - 1
  260     CONTINUE
        END IF
  380 CONTINUE
      END
      SUBROUTINE MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),ICNTL(20)
      INTRINSIC MIN,SIGN
      CHARACTER*72 LINE
      INTEGER APOS,APOS2,IBLK,ILINE,IROW,IWPOS,J,JPIV,J1,J2,K,
     +        LDIAG,LEN,MP,NBLK,NCOLS,NROWS
      CHARACTER*1 PM(-2:2)
      DATA PM/'*','-','.','+','.'/
      DOUBLE PRECISION ZERO,TINY,FD15AD
      PARAMETER (ZERO=0.0D0)
      EXTERNAL FD15AD
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      TINY = FD15AD('T')
      APOS2 = IFACT(1)
      NBLK = IFACT(3)
      IF (LDIAG.EQ.3) NBLK = MIN(1,NBLK)
      LEN = 12
      IF (LDIAG.EQ.5) LEN = 1
      IF (LEN.EQ.12) THEN
        IF (NBLK.EQ.IFACT(3)) THEN
          WRITE (MP,'(/A)')
     +      'For each block, the following information is provided:'
        ELSE
          WRITE (MP,'(/A,A)') 'For the first block only,',
     +      ' the following information is provided:'
        END IF
      END IF
      IF (LEN.EQ.12) WRITE (MP,'(A)')
     +    '   1. Block number, number of rows, number of columns',
     +    '   2. List of indices for the pivot, each negated if part of'
     +    ,'      a 2x2 pivot',
     +    '   3. The factorized block pivot',
     +    '      It has the form',
     +    '            -1  T',
     +    '        L  D   L ',
     +    '                         -1    T',
     +    '      and is printed as D and L  packed together.',
     +    '   4. List of indices for the non-pivot columns',
     +    '   5. The non-pivot part as rectangular block by rows'
      IWPOS = 4
      APOS = 1
      DO 300 IBLK = 1,NBLK
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        WRITE (MP,'(/4(A,I6))') 'Block pivot',IBLK,' with',NROWS,
     +        ' rows and', NCOLS,' columns'
        IF (LEN.EQ.12) WRITE (MP,'(6I12)')
     +                       (IFACT(K),K=IWPOS,IWPOS+NROWS-1)
        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (PM(SIGN(1,IFACT(K))),
     +      K=IWPOS,IWPOS+NROWS-1)
        JPIV = 0
        DO 30 IROW = 1,NROWS
          IF (JPIV.EQ.1) THEN
            JPIV = 0
          ELSE
            IF (IFACT(IWPOS+IROW-1).LT.0) JPIV = 1
          END IF
          ILINE = 1
          DO 10 J = 1,IROW - 1
            WRITE (LINE(ILINE:ILINE+LEN-1),'(A)') ' '
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
   10     CONTINUE
          DO 20 J = IROW,NROWS
            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +          '(1P,D12.4)') FACT(APOS)
            IF (LEN.EQ.1) THEN
               IF (FACT(APOS).EQ.ZERO) THEN
                  WRITE (LINE(ILINE:ILINE),'(A)') '.'
               ELSE
                  WRITE (LINE(ILINE:ILINE),'(A)') '*'
               END IF
            END IF
            APOS = APOS + 1
            IF (J.EQ.IROW+1) THEN
              IF (JPIV.EQ.1) THEN
                IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +              '(1P,D12.4)') FACT(APOS2)
                IF (LEN.EQ.1) THEN
                    IF (FACT(APOS2).EQ.ZERO) THEN
                       WRITE (LINE(ILINE:ILINE),'(A)') '.'
                    ELSE
                       WRITE (LINE(ILINE:ILINE),'(A)') '*'
                    END IF
                END IF
                APOS2 = APOS2 + 1
              END IF
            END IF
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
   20     CONTINUE
          IF (ILINE.GT.1) THEN
            LINE(ILINE:) = ' '
            WRITE (MP,'(A)') LINE
          END IF
   30   CONTINUE
        IWPOS = IWPOS + NROWS
        IF (LEN.EQ.12) WRITE (MP,'(6I12)') (IFACT(K),K=IWPOS,
     +      IWPOS+NCOLS-NROWS-1)
        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (PM(SIGN(1,IFACT(K))),
     +      K=IWPOS,IWPOS+NCOLS-NROWS-1)
        IWPOS = IWPOS + NCOLS - NROWS
        DO 280 IROW = 1,NROWS
          J1 = NROWS
          J2 = NCOLS
          ILINE = 1
          DO 110 J = J1 + 1,J2
            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +          '(1P,D12.4)') FACT(APOS)
            IF (LEN.EQ.1) THEN
               IF (FACT(APOS).EQ.ZERO) THEN
                  WRITE (LINE(ILINE:ILINE),'(A)') '.'
               ELSE
                  WRITE (LINE(ILINE:ILINE),'(A)') '*'
               END IF
            END IF
            APOS = APOS + 1
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
  110     CONTINUE
          IF (ILINE.GT.1) THEN
            LINE(ILINE:) = ' '
            WRITE (MP,'(A)') LINE
          END IF
  280   CONTINUE
  300 CONTINUE
      END
      SUBROUTINE MA57SD(FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,ICNTL)
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = 1
      APOS2 = IFACT(1)
      IWPOS = 4
      DO 380 IBLK = 1,IFACT(3)
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 10 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 9 J = 1,NRHS
              W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9       CONTINUE
            APOS = APOS + (NROWS+1-IPIV)
   10     CONTINUE
          JPIV = 1
          DO 20 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              DO 19 J = 1,NRHS
                W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19         CONTINUE
              IF (JPIV.EQ.-1) APOS2 = APOS2 + 1
              JPIV = -JPIV
            END IF
   20     CONTINUE
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE
        ELSE
          JPIV = 1
          DO 210 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.GT.0) THEN
              DO 65 J = 1,NRHS
                RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65         CONTINUE
              APOS = APOS + NROWS - IPIV + 1
            ELSE
              IF (JPIV.EQ.1) THEN
                IRHS1 = -IRHS
                IRHS2 = -IFACT(IWPOS+IPIV)
                DO 68 J = 1,NRHS
                  W1 = RHS(IRHS1,J)*FACT(APOS) +
     +                 RHS(IRHS2,J)*FACT(APOS2)
                  RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
     +                           RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
                  RHS(IRHS1,J) = W1
   68           CONTINUE
                APOS2 = APOS2 + 1
              END IF
              JPIV = -JPIV
              APOS = APOS + NROWS - IPIV + 1
            END IF
  210     CONTINUE
        END IF
        IWPOS = IWPOS + NCOLS
        APOS = APOS + NROWS*(NCOLS-NROWS)
  380 CONTINUE
      END
      SUBROUTINE MA57TD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,
     +        IWPOS,J,J1,J2,K,KK,NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = IFACT(1)
      IWPOS = 4
      DO 10 I = 1,IFACT(3)-1
        IW1(I) = IWPOS
        IWPOS = IWPOS + ABS(IFACT(IWPOS))+2
   10 CONTINUE
      IW1(IFACT(3)) = IWPOS
      DO 380 IBLK = IFACT(3),1,-1
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 5 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 3 J = 1,NRHS
              W(I,J) = RHS(II,J)
    3       CONTINUE
    5     CONTINUE
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMM('T','N',NROWS,NRHS,K,ONE,
     +                           FACT(APOS),K,
     +                           W(NROWS+1,1),LW,ONE,W,LW)
          APOS = APOS-(NROWS*(NROWS+1))/2
          DO 22 J = 1,NRHS
            CALL DTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE
        ELSE
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1
          KK = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 218 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = KK
              DO 215 J = J1,J2
                W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  215         CONTINUE
              RHS(IRHS,II) = W1
  218       CONTINUE
            KK = K
  220     CONTINUE
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            DO 240 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = APOS + 1
              DO 230 J = J1,J2
                W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  230         CONTINUE
              RHS(IRHS,II) = W1
  240       CONTINUE
            J1 = J1 - 1
  260     CONTINUE
        END IF
  380 CONTINUE
      END
      SUBROUTINE MA57DD(JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,
     *                  RHS,X,RESID,W,IW,ICNTL,CNTL,INFO,RINFO)
      INTEGER JOB,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),JCN(NE),LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT)
      DOUBLE PRECISION RHS(N),X(N),RESID(N),W(N,*)
      INTEGER IW(N),ICNTL(20)
      DOUBLE PRECISION CNTL(5)
      INTEGER INFO(40)
      DOUBLE PRECISION RINFO(20)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.0D0)
      DOUBLE PRECISION COND(2),CTAU,DXMAX,ERROR,OLDOMG(2),OLDOM2,
     *                 OMEGA(2),OM2,TAU
      INTEGER I,ICNTLC(20),ITER,J,K,KASE,KK,LDIAG,LP,MP,KEEP71(5)
      LOGICAL LCOND(2)
      INTRINSIC MIN
      EXTERNAL MA57CD,MA57UD,FD15AD,MC71AD
      DOUBLE PRECISION EPS,FD15AD
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      INFO(1) = 0
      IF (N.LE.0) THEN
        INFO(1) = -1
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'N has value',N
        INFO(2) = N
        GOTO 500
      ENDIF
      IF (NE.LT.0) THEN
        INFO(1) = -2
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'NE has value',NE
        INFO(2) = NE
        GOTO 500
      ENDIF
      IF (ICNTL(9).LT.1) THEN
        INFO(1) = -13
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'ICNTL(9) has value',ICNTL(9)
        INFO(2) = ICNTL(9)
        GOTO 500
      ENDIF
      IF (JOB.LT.0 .OR. JOB.GT.4 .OR. (ICNTL(9).GT.1 .AND.
     *    (JOB.NE.0 .AND. JOB.NE.2)))  THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'JOB has value',JOB
        IF (ICNTL(9).GT.1 .AND. LDIAG.GT.0 .AND. LP.GE.0)
     +    WRITE (LP,'(A,I3)') 'and ICNTL(9) =',ICNTL(9)
        GOTO 500
      ENDIF
      IF (NE.EQ.0) THEN
        IF (JOB.NE.3) THEN
          DO 8 I = 1,N
            RESID(I) = ZERO
  8       CONTINUE
        ENDIF
        DO 9 I = 1,N
          X(I) = ZERO
  9     CONTINUE
        INFO(30)=0
        DO 10 I = 6,13
          RINFO(I) = ZERO
 10     CONTINUE
        GO TO 500
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999) JOB,N,NE,(ICNTL(I),I=1,5),LFACT,LIFACT,
     +   ICNTL(9),ICNTL(10),ICNTL(13),CNTL(3)
99999 FORMAT(/'Entering iterative refinement solution phase ',
     +  '(MA57DD) with ...'/
     +  'JOB       Control for coefficient matrix      =',I12/
     +  'N         Order of matrix                     =',I12/
     +  'NE        Number of entries in matrix         =',I12/
     6  'ICNTL(1)  Stream for errors                   =',I12/
     7  ' --- (2)  Stream for warnings                 =',I12/
     8  ' --- (3)  Stream for monitoring               =',I12/
     9  ' --- (4)  Stream for statistics               =',I12/
     1  ' --- (5)  Level of diagnostic printing        =',I12/
     +  'LFACT     Length of array FACT                =',I12/
     +  'LIFACT    Length of array IFACT               =',I12/
     +  'ICNTL(9)  Number steps iterative refinement   =',I12/
     +  'ICNTL(10) Control for error analysis          =',I12/
     +  'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12/
     +  'CNTL(3)   Convergence test for IR             =',1P,D12.4)
        CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE(MP,'(/A)') 'Right-hand side'
        WRITE (MP,'((4X, 1P,5D13.3))') (RHS(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
      END IF
      DO 15 I=1,5
        ICNTLC(I) = ICNTL(I)
   15 CONTINUE
      ICNTLC(13) = ICNTL(13)
      ICNTLC(15) = ICNTL(15)
      ICNTLC(3) = -1
      IF (JOB.LE.2) THEN
        IF (JOB .LE. 1) THEN
          DO 14 I = 1,N
            X(I) = RHS(I)
            RESID(I) = RHS(I)
   14     CONTINUE
          CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,X,N,W,N,IW,
     +                ICNTLC,INFO)
        ELSE
          DO 13 I = 1,N
            RESID(I) = RHS(I)
   13     CONTINUE
        ENDIF
        IF (ICNTL(9).EQ.1) THEN
          DO 16 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 16
            RESID(J) = RESID(J) - A(KK)*X(I)
            IF (I.NE.J) RESID(I) = RESID(I) - A(KK)*X(J)
   16     CONTINUE
          IF (JOB.EQ.0) GO to 340
        ELSE
          DO 18 I = 1,N
            W(I,1) = ZERO
            W(I,3) = ZERO
   18     CONTINUE
          DO 17 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 17
            RESID(J) = RESID(J) - A(KK)*X(I)
            W(J,1) = W(J,1) + ABS(A(KK)*X(I))
            W(J,3) = W(J,3) + ABS(A(KK))
            IF (I.NE.J) THEN
              RESID(I) = RESID(I) - A(KK)*X(J)
              W(I,1) = W(I,1) + ABS(A(KK)*X(J))
              W(I,3) = W(I,3) + ABS(A(KK))
            ENDIF
   17     CONTINUE
        DXMAX = ZERO
        DO 221 I = 1,N
          DXMAX = MAX(DXMAX,ABS(X(I)))
  221   CONTINUE
      EPS = FD15AD('E')
        CTAU = 1000.*EPS
          OMEGA(1) = ZERO
          OMEGA(2) = ZERO
          DO 231 I = 1,N
            TAU = (W(I,3)*DXMAX+ABS(RHS(I)))*N*CTAU
            IF ((W(I,1)+ABS(RHS(I))).GT.TAU) THEN
              OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
     +                   (W(I,1)+ABS(RHS(I))))
              IW(I) = 1
            ELSE
              IF (TAU.GT.ZERO) THEN
                OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
     +                     (W(I,1)+W(I,3)*DXMAX))
              END IF
              IW(I) = 2
            END IF
  231     CONTINUE
          OM2 = OMEGA(1) + OMEGA(2)
          ITER = 0
          IF (OM2.LE.EPS) THEN
            GO TO 270
          ENDIF
          DO 251 I = 1,N
            W(I,2) = X(I)
  251     CONTINUE
          OLDOMG(1) = OMEGA(1)
          OLDOMG(2) = OMEGA(2)
          OLDOM2 = OM2
        ENDIF
      ENDIF
      DO 260 ITER = 1,ICNTL(9)
        CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,RESID,N,W,N,IW,
     +              ICNTLC,INFO)
        DO 141 I = 1,N
          X(I) = X(I) + RESID(I)
  141   CONTINUE
        IF (JOB.LT.4 .AND. ICNTL(9).EQ.1) GO TO 340
        IF (ICNTL(9).EQ.1) THEN
          DO 151 I = 1,N
            RESID(I) = RHS(I)
  151     CONTINUE
          DO 181 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 181
            RESID(J) = RESID(J) - A(KK)*X(I)
            IF (I.NE.J) RESID(I) = RESID(I) - A(KK)*X(J)
  181     CONTINUE
          GO TO 340
        ELSE
          DO 153 I = 1,N
            RESID(I) = RHS(I)
            W(I,1) = ZERO
  153     CONTINUE
          DO 183 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 183
            RESID(J) = RESID(J) - A(KK)*X(I)
            W(J,1) = W(J,1) + ABS(A(KK)*X(I))
            IF (I.NE.J) THEN
              RESID(I) = RESID(I) - A(KK)*X(J)
              W(I,1) = W(I,1) + ABS(A(KK)*X(J))
            ENDIF
  183     CONTINUE
        ENDIF
        DXMAX = ZERO
        DO 220 I = 1,N
          DXMAX = MAX(DXMAX,ABS(X(I)))
  220   CONTINUE
        OMEGA(1) = ZERO
        OMEGA(2) = ZERO
        DO 230 I = 1,N
          TAU = (W(I,3)*DXMAX+ABS(RHS(I)))*N*CTAU
          IF ((W(I,1)+ABS(RHS(I))).GT.TAU) THEN
            OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
     +                 (W(I,1)+ABS(RHS(I))))
            IW(I) = 1
          ELSE
            IF (TAU.GT.ZERO) THEN
              OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
     +                   (W(I,1)+W(I,3)*DXMAX))
            END IF
            IW(I) = 2
          END IF
  230   CONTINUE
        OM2 = OMEGA(1) + OMEGA(2)
        IF ((OM2+ONE).LE.ONE) THEN
          GO TO 270
        ENDIF
        IF (OM2.GT.OLDOM2*CNTL(3)) THEN
          IF (OM2.GT.OLDOM2) THEN
            OMEGA(1) = OLDOMG(1)
            OMEGA(2) = OLDOMG(2)
            DO 240 I = 1,N
              X(I) = W(I,2)
  240       CONTINUE
          END IF
          GO TO 270
        ELSE
          DO 250 I = 1,N
            W(I,2) = X(I)
  250     CONTINUE
          OLDOMG(1) = OMEGA(1)
          OLDOMG(2) = OMEGA(2)
          OLDOM2 = OM2
        END IF
  260 CONTINUE
      INFO(1) = -8
      IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
 9170 FORMAT ('Error return from MA57D/DD because of ','nonconvergence',
     +       ' of iterative refinement'/'Error INFO(1) = ',I2,'  with',
     +       ' ICNTL','(9) = ',I10)
  270 RINFO(6)  = OMEGA(1)
      RINFO(7)  = OMEGA(2)
      RINFO(8) = ZERO
      DO 271 I=1,N
        RINFO(8) = MAX(RINFO(8),W(I,3))
  271 CONTINUE
      RINFO(9) = DXMAX
      RINFO(10) = ZERO
      DO 272 I=1,N
        RINFO(10) = MAX(RINFO(10),ABS(RESID(I)))
  272 CONTINUE
      IF (RINFO(8)*RINFO(9).NE.ZERO)
     *RINFO(10) = RINFO(10)/(RINFO(8)*RINFO(9))
      INFO(30) = ITER
      IF (INFO(1).LT.0) GO TO 340
      IF (ICNTL(10).LE.0) GO TO 340
      LCOND(1) = .FALSE.
      LCOND(2) = .FALSE.
      ERROR    = ZERO
      DO 280 I = 1,N
        IF (IW(I).EQ.1) THEN
          W(I,1) = W(I,1) + ABS(RHS(I))
          W(I,2) = ZERO
          LCOND(1) = .TRUE.
        ELSE
          W(I,2) = W(I,1) + W(I,3)*DXMAX
          W(I,1) = ZERO
          LCOND(2) = .TRUE.
        END IF
  280 CONTINUE
      DO 330 K = 1,2
        IF (LCOND(K)) THEN
          KASE = 0
          DO 310 KK = 1,40
            CALL MC71AD(N,KASE,W(1,3),COND(K),W(1,4),IW,KEEP71)
            IF (KASE.EQ.0) GO TO 320
            IF (KASE.EQ.1) THEN
              CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,W(1,3),
     *                    N,W(1,4),N,IW,ICNTLC,INFO)
              DO 290 I = 1,N
                W(I,3) = W(I,K)*W(I,3)
  290         CONTINUE
            END IF
            IF (KASE.EQ.2) THEN
              DO 300 I = 1,N
                W(I,3) = W(I,K)*W(I,3)
  300         CONTINUE
              CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,W(1,3),N,
     *                    W(1,4),N,IW,ICNTLC,INFO)
            END IF
  310     CONTINUE
          INFO(1) = -14
          IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9160)
 9160 FORMAT ('Error return from MA57D/DD because of ','error in MC71',
     +       'A/AD'/'Error not calculated')
  320     IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
          ERROR = ERROR + OMEGA(K)*COND(K)
        ELSE
          COND(K) = ZERO
        ENDIF
  330 CONTINUE
      RINFO(11)  = COND(1)
      RINFO(12)  = COND(2)
      RINFO(13)  = ERROR
 340  IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE(MP,99980) INFO(1)
99980 FORMAT (/'Leaving iterative refinement solution phase ',
     +  '(MA57DD) with ...'/
     1      'INFO (1)                                      =',I12/)
        IF (INFO(1).LT.0) GO TO 500
        IF (ICNTL(9).GT.1) THEN
          WRITE(MP,99981) INFO(30),(RINFO(I),I=6,10)
99981     FORMAT(
     1     'INFO(30)  Number steps iterative ref   =',I10/
     1     'RINFO(6)  Backward errors  (OMEGA(1))  =',1PD10.3/
     2     '-----(7)  Backward errors  (OMEGA(2))  =',1PD10.3/
     3     '-----(8)  Infinity norm of matrix      =',1PD10.3/
     4     '-----(9)  Infinity norm of solution    =',1PD10.3/
     5     '-----(10) Norm of scaled residuals     =',1PD10.3)
          IF (ICNTL(10).GT.0) WRITE(MP,99979) (RINFO(I),I=11,13)
99979       FORMAT (
     1       'RINFO(11) Condition number (COND(1))   =',1PD10.3/
     1       'RINFO(12) Condition number (COND(2))   =',1PD10.3/
     1       'RINFO(13) Error in solution            =',1PD10.3)
          WRITE(MP,'(/A,I10)') 'Residual'
          K=MIN(N,10)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(1P,5D13.3)') (RESID(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
        ELSE
          IF (JOB.GE.1 .AND. JOB.LE.3) THEN
            WRITE(MP,'(/A,I10)') 'Correction to solution'
          ELSE
            WRITE(MP,'(/A,I10)') 'Residual'
          ENDIF
          K=MIN(N,10)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(1P,5D13.3)') (RESID(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
        END IF
        K=MIN(N,10)
        IF (LDIAG.GE.4) K = N
        WRITE(MP,'(/A,I10)') 'Solution'
        WRITE (MP,'(1P,5D13.3)') (X(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
      END IF
 500  RETURN
      END
      SUBROUTINE MA57ED(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,
     *                  IFACT,LIFACT,NEWIFC,LINEW,INFO)
      INTEGER N,IC,KEEP(*),LFACT,LNEW,LIFACT,LINEW,INFO(40)
      DOUBLE PRECISION FACT(LFACT),NEWFAC(LNEW)
      INTEGER IFACT(LIFACT),NEWIFC(LINEW)
      INTEGER APOSBB,ASTK,HOLD,I,ISTK,IWPOS,MOVE,NFRONT
      HOLD = N + 3
      INFO(1) = 0
      INFO(2) = 0
      IF (IC.GE.1) THEN
        IF (LINEW.LE.LIFACT) THEN
          INFO(1) = -7
          INFO(2) = LINEW
          RETURN
        ENDIF
        IWPOS = KEEP(HOLD+7)
        ISTK  = KEEP(HOLD+14)
        NFRONT = KEEP(HOLD+23)
        DO 10 I = 1,IWPOS+NFRONT-1
          NEWIFC(I) = IFACT(I)
   10   CONTINUE
        MOVE = LINEW - LIFACT
        DO 20 I = ISTK+1,LIFACT
          NEWIFC(I+MOVE) = IFACT(I)
   20   CONTINUE
          KEEP(HOLD+13) = KEEP(HOLD+13) + MOVE
          KEEP(HOLD+14) = ISTK + MOVE
          KEEP(HOLD+18) = KEEP(HOLD+18) + MOVE
      ENDIF
      IF (IC.NE.1) THEN
        IF (LNEW.LE.LFACT) THEN
          INFO(1) = -7
          INFO(2) = LNEW
          RETURN
        ENDIF
        APOSBB = KEEP(HOLD+9)
        ASTK   = KEEP(HOLD+15)
        DO 60 I = 1, APOSBB-1
          NEWFAC(I) = FACT(I)
   60   CONTINUE
        MOVE = LNEW - LFACT
        DO 70 I = ASTK+1,LFACT
          NEWFAC(I+MOVE) = FACT(I)
   70   CONTINUE
        KEEP(HOLD+12) = KEEP(HOLD+12) + MOVE
        KEEP(HOLD+15) = ASTK + MOVE
        KEEP(HOLD+19) = KEEP(HOLD+19) + MOVE
      ENDIF
      RETURN
      END
      SUBROUTINE MA57GD(N,NE,IRN,JCN,IW,IPE,COUNT,FLAG,IWFR,ICNTL,INFO)
      INTEGER N,NE,IRN(NE),JCN(NE),IW(NE*2+N),IPE(N),COUNT(N),
     +        FLAG(N),IWFR,ICNTL(20),INFO(40)
      INTRINSIC MAX,MIN
      INTEGER I,J,K,L,LDIAG,WP
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      IF (WP.LT.0) LDIAG = 0
      INFO(1) = 0
      INFO(3) = 0
      DO 10 I = 1,N
        FLAG(I) = 0
        COUNT(I) = 0
   10 CONTINUE
      DO 20 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          INFO(3) = INFO(3) + 1
          INFO(1) = 1
          IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
     +        '*** Warning message from subroutine MA57AD ***',
     +        ' INFO(1) =',INFO(1)
          IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
     +         K,'th entry (in row',I,' and column',J,') ignored'
        ELSE IF (I.NE.J) THEN
          COUNT(I) = COUNT(I) + 1
          COUNT(J) = COUNT(J) + 1
        END IF
   20 CONTINUE
      IPE(1) = COUNT(1)+1
      DO 30 I = 2,N
        IPE(I) = IPE(I-1) + COUNT(I)
   30 CONTINUE
      DO 40 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N .OR. I.EQ.J) GO TO 40
        IPE(I) = IPE(I) - 1
        IW(IPE(I)) = J
        IPE(J) = IPE(J) - 1
        IW(IPE(J)) = I
   40 CONTINUE
      INFO(4) = 0
      IWFR = 1
      DO 60 I = 1,N
        L = IPE(I)
        IPE(I) = IWFR
        DO 50 K = L,L+COUNT(I)-1
          J = IW(K)
          IF (FLAG(J).NE.I) THEN
            FLAG(J) = I
            IW(IWFR) = J
            IWFR = IWFR + 1
          ELSE
            IF (I.LT.J) INFO(4) = INFO(4) + 1
          END IF
   50   CONTINUE
        COUNT(I) = IWFR - IPE(I)
   60 CONTINUE
      IF (INFO(4).GT.0) THEN
        INFO(1) = INFO(1) + 2
        IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
     +      '*** Warning message from subroutine MA57AD ***',INFO(4),
     +      ' off-diagonal duplicate entries found'
      END IF
      END
      SUBROUTINE MA57JD(N,NE,IRN,JCN,PERM,IW,IPE,COUNT,FLAG,IWFR,
     +                  ICNTL,INFO)
      INTEGER N,NE,IRN(NE),JCN(NE),IW(NE+N),IPE(N),COUNT(N),
     +        PERM(N),FLAG(N),IWFR,ICNTL(20),INFO(40)
      INTRINSIC MAX,MIN
      INTEGER I,J,K,L,LDIAG,WP
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      IF (WP.LT.0) LDIAG = 0
      INFO(1) = 0
      DO 10 I = 1,N
        FLAG(I) = 0
        COUNT(I) = 0
   10 CONTINUE
      INFO(3) = 0
      DO 30 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          IRN(K) = 0
          JCN(K) = 0
          INFO(3) = INFO(3) + 1
          INFO(1) = 1
          IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
     +        '*** Warning message from subroutine MA57AD ***',
     +        ' INFO(1) =',INFO(1)
          IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
     +        K,'th entry (in row',I,' and column',J,') ignored'
        ELSE IF (PERM(I).LE.PERM(J)) THEN
          COUNT(I) = COUNT(I) + 1
        ELSE
          COUNT(J) = COUNT(J) + 1
        END IF
   30 CONTINUE
      IPE(1) = COUNT(1) + 1
      DO 40 I = 2,N
        IPE(I) = IPE(I-1) + COUNT(I) + 1
   40 CONTINUE
      DO 50 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 50
        IF (PERM(I).LE.PERM(J)) THEN
          IW(IPE(I)) = J
          IPE(I) = IPE(I) - 1
        ELSE
          IW(IPE(J)) = I
          IPE(J) = IPE(J) - 1
        END IF
   50 CONTINUE
      IWFR = 1
      INFO(4) = 0
      DO 70 I = 1,N
        L = IPE(I)
        IPE(I) = IWFR
        DO 60 K = L + 1,L + COUNT(I)
          J = IW(K)
          IF (FLAG(J).NE.I) THEN
            FLAG(J) = I
            IWFR = IWFR + 1
            IW(IWFR) = J
          ELSE
            IF (I.LT.J) INFO(4) = INFO(4) + 1
          END IF
   60   CONTINUE
        IF (IWFR.GT.IPE(I)) THEN
          IW(IPE(I)) = IWFR - IPE(I)
          IWFR = IWFR + 1
        ELSE
          IPE(I) = 0
        END IF
   70 CONTINUE
      IF (INFO(4).GT.0) THEN
        INFO(1) = INFO(1) + 2
        IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
     +      '*** Warning message from subroutine MA57AD ***',
     +      INFO(4),' off-diagonal duplicate entries found'
      END IF
      END
      SUBROUTINE MA57KD(N, IPE, IW, LW, IWFR, PERM, IPS, NV, FLAG,
     *                  NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER IW(LW), PERM(N), IPS(N), NV(N), FLAG(N)
      INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
      INTEGER LN,JP1,JS,LWFR,JP2,JE
      EXTERNAL MA57FD
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = PERM(I)
        IPS(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML=1,N
        MS = IPS(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
          DO 50 JP1=1,LN
            JP = JP + 1
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA57FD(N, IPE, IW, IP-1, LWFR, NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = MIN0(MINJS,PERM(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
   60     IPE(IE) = -ME
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPS(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      RETURN
      END
C** end of MA57KD**
      SUBROUTINE MA57FD(N, IPE, IW, LW, IWFR, NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER   IW(LW)
      INTEGER I,K1,LWFR,IR,K,K2
      NCMPA = NCMPA + 1
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
C--------------------------------------------------------------------
C-             Copyright CCLRC Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA57LD(N, IPE, NV, IPS, NE, NA, NODE, PERM, NSTEPS,
     *                  FILS, FRERE, ND, NEMIN, SUBORD)
      INTEGER N, NSTEPS
      INTEGER ND(N)
      INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
      INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N), PERM(N)
      INTEGER NEMIN
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,IFSON,INO
      INTEGER INOS,IB,IL,INT
      INTEGER IPERM
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
        NODE(I) = 0
        SUBORD(I) = 0
   10 CONTINUE
      NR = N + 1
      DO 50 I=1,N
        IF = -IPE(I)
        IF (NV(I).EQ.0) THEN
          IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
          SUBORD(IF) = I
          NODE(IF) = NODE(IF)+1
        ELSE
          IF (IF.NE.0) THEN
            IS = -IPS(IF)
            IF (IS.GT.0) IPE(I) = IS
            IPS(IF) = -I
          ELSE
            NR = NR - 1
            NE(NR) = I
          ENDIF
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (IPE(INS).LT.0) THEN
       INS       = -IPE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (IPE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = IPE(INS)
C?? I think this test is the wrong way round
      IF (NV(INB).GE.NV(INS)) THEN
C?? So reversed
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
        FILS(INF) = -INB
        IPS(INF)  = -INB
        IPE(INS)  = IPE(INB)
        IPE(INB)  = INS
      ELSE
        INSW = INFS
 1100   INFS = IPE(INSW)
        IF (INFS.NE.INS) THEN
          INSW = INFS
          GO TO 1100
        ENDIF
        IPE(INS) = IPE(INB)
        IPE(INB) = INS
        IPE(INSW)= INB
      ENDIF
        INS      = INB
        GO TO 1070
 1151 DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I = 0
      IPERM = 1
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
        IF (NR.GT.N) GO TO 161
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   CONTINUE
        DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   CONTINUE
C?? Do we want to expand for subordinate variables
        IPS(I) = K
        NE(IS) = NE(IS) + NODE(I) + 1
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        NODE(I) = IS
        PERM(I) = IPERM
        IPERM = IPERM + 1
        IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
          IN = SUBORD(IN)
          NODE(IN) = IS
          PERM(IN) = IPERM
          IPERM = IPERM + 1
          GO TO 777
  778   IF (NA(IS).NE.1) GO TO 90
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
   90   IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NA(IS).EQ.0) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        NODE(I) = IS-1
        IFSON = -FILS(I)
        IN = IFSON
 102    INO = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 102
        NV(INO) = 0
        IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
        IN = SUBORD(IN)
        NODE(IN) = IS-1
        GO TO 888
  889   SUBORD(IN) = INO
        IN = INO
        IF (SUBORD(IN).EQ.0) GO TO 887
        IN = SUBORD(IN)
        IPE(IN) = -I
  887   CONTINUE
      INOS = -FILS(INO)
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
        IF (INOS.EQ.0) THEN
          FRERE(INS) = -I
          GO TO 120
        ELSE
          FRERE(INS) =  INOS
        ENDIF
 107    IN = INOS
        IF (IN.EQ.0) GO TO 120
 108    INT = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 108
        FRERE(INT) = -I
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
          GO TO 160
        ELSE
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
  161 NSTEPS = IS - 1
      RETURN
      END
      SUBROUTINE MA57MD(N,NE,IRN,JCN,MAP,IRNPRM,
     +                  LROW,PERM,COUNT,IDIAG)
      INTEGER N,NE
      INTEGER IRN(NE),JCN(NE),MAP(NE),IRNPRM(N+NE),LROW(N),PERM(N),
     +        COUNT(N),
     +        IDIAG(N)
      INTEGER EXPNE,I,J,K
      DO 10 I = 1,N
        COUNT(I) = 1
        IDIAG(I) = 0
   10 CONTINUE
      EXPNE = NE + N
      DO 20 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MAX(I,J).GT.N .OR. MIN(I,J).LT.1) THEN
          EXPNE = EXPNE - 1
          GO TO 20
        ENDIF
        IF (I.EQ.J) THEN
          I = PERM(I)
          IF (IDIAG(I).GE.1) THEN
            COUNT(I) = COUNT(I) + 1
            IDIAG(I) = IDIAG(I) + 1
          ELSE
            IDIAG(I) = 1
            EXPNE = EXPNE - 1
          ENDIF
          GO TO 20
        ENDIF
        IF (PERM(I).LT.PERM(J)) THEN
          I = PERM(I)
          COUNT(I) = COUNT(I) + 1
        ELSE
          J = PERM(J)
          COUNT(J) = COUNT(J) + 1
        END IF
   20 CONTINUE
      LROW(1) = COUNT(1)
      IDIAG(1) = MAX(IDIAG(1),1)
      DO 30 I = 2,N
        LROW(I) = COUNT(I)
        COUNT(I) = COUNT(I-1) + LROW(I)
        IDIAG(I) = COUNT(I-1) + MAX(IDIAG(I),1)
   30 CONTINUE
      DO 35 I = 1,N
        K = PERM(I)
        IRNPRM(IDIAG(K)) = I
   35 CONTINUE
      DO 40 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          MAP(K) = 0
          GO TO 40
        ENDIF
        I = PERM(IRN(K))
        J = PERM(JCN(K))
        IF (I.EQ.J) THEN
          MAP(K) = IDIAG(I)
          IRNPRM(IDIAG(I)) = IRN(K)
          IDIAG(I) = IDIAG(I) - 1
        ELSE
          IF (I.GT.J) THEN
            MAP(K) = COUNT(J)
            IRNPRM(COUNT(J)) = IRN(K)
            COUNT(J) = COUNT(J) - 1
          ELSE
            MAP(K) = COUNT(I)
            IRNPRM(COUNT(I)) = JCN(K)
            COUNT(I) = COUNT(I) - 1
          ENDIF
        ENDIF
   40 CONTINUE
      IDIAG(1) = EXPNE
      RETURN
      END
      SUBROUTINE MA57ND(N,LENR,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
     *                  INFO,RINFO)
      INTEGER N,NSTEPS
      INTEGER LENR(N),LSTKI(N),LSTKR(N),NA(NSTEPS),
     +        ND(NSTEPS),NE(NSTEPS),INFO(40)
      DOUBLE PRECISION RINFO(20)
      INTEGER I,IORG,ISTKI,ISTKR,ITOP,ITREE,JORG,K,
     +        LSTK,NASSR,NELIM,NFR,NSTK,NTOTPV,NZ1,NZ2
      DOUBLE PRECISION DELIM
      INTRINSIC MAX
      DOUBLE PRECISION OPS,OPSASS
      INTEGER NIRADU,NIRNEC,NIRTOT,NRLADU,NRLNEC,NRLTOT,MAXFRT
      NZ1 = 0
      DO 40 I = 1,N
        NZ1 = NZ1 + LENR(I)
   40 CONTINUE
      NZ2 = NZ1
      ISTKI = 0
      ISTKR = 0
      OPS = 0.0D0
      OPSASS = 0.0D0
      NRLADU = 0
      NIRADU = 3
      NIRTOT = NZ1+N+5
      NRLTOT = NZ1
      NIRNEC = NZ2+N+5
      NRLNEC = NZ2
      NTOTPV = 0
      ITOP = 0
      MAXFRT = 0
      DO 100 ITREE = 1,NSTEPS
        NELIM = NE(ITREE)
        DELIM = NELIM
        NFR = ND(ITREE)
        MAXFRT = MAX(MAXFRT,NFR)
        NSTK = NA(ITREE)
        NASSR = NELIM*(NELIM+1)/2 + NFR*NFR
        NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
        NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)
        DO 70 IORG = 1,NELIM
          JORG = NTOTPV + IORG
          OPSASS = OPSASS + LENR(JORG)
          NZ2 = NZ2 - LENR(JORG)
   70   CONTINUE
        NTOTPV = NTOTPV + NELIM
        DO 80 K = 1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          OPSASS = OPSASS + LSTK
          LSTK = LSTKI(ITOP)
          ISTKI = ISTKI - LSTK
          ITOP = ITOP - 1
   80   CONTINUE
        NRLADU = NRLADU + (NELIM*(NELIM+1))/2 + (NFR-NELIM)*NELIM
        NIRADU = NIRADU + 2 + NFR
        OPS = OPS + (DELIM* (12*NFR+6*NFR*NFR - (DELIM+1)*
     +        (6*NFR+6-(2*DELIM+1))))/6 + DELIM
        IF (NFR.GT.NELIM) THEN
          ITOP = ITOP + 1
          LSTKR(ITOP) = ((NFR-NELIM)*(NFR-NELIM+1))/2
          LSTKI(ITOP) = NFR - NELIM + 1
          ISTKI = ISTKI + LSTKI(ITOP)
          ISTKR = ISTKR + LSTKR(ITOP)
        ENDIF
        IF (ITREE.EQ.NSTEPS) THEN
          NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
          NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
        ELSE
          NIRTOT = MAX(NIRTOT,NIRADU+(N-NTOTPV+2)+ISTKI+NZ1)
          NIRNEC = MAX(NIRNEC,NIRADU+(N-NTOTPV+2)+ISTKI+NZ2)
        ENDIF
  100 CONTINUE
      INFO(5)   = NRLADU
      INFO(6)   = NIRADU
      INFO(7)   = MAXFRT
      INFO(8)   = NSTEPS
      INFO(9)   = NRLTOT
      INFO(10)  = NIRTOT
      INFO(11)  = NRLNEC
      INFO(12)  = NIRNEC
      RINFO(1)  = OPSASS
      RINFO(2)  = OPS
      RETURN
      END
      SUBROUTINE MA57OD(N,NE,A,LA,IW,LIW,LROW,PERM,NSTEPS,NSTK,NODE,
     +                  DIAG,SCHNAB,PPOS,CNTL,ICNTL,INFO,RINFO,HOLD,
     +                  BIGA)
      INTEGER N,NE,LA
      DOUBLE PRECISION A(LA),DIAG(N),SCHNAB(*),CNTL(5),RINFO(20),BIGA
      INTEGER LIW,IW(LIW),LROW(N),PERM(N),NSTEPS,NSTK(NSTEPS),
     +        NODE(N),PPOS(N),ICNTL(20),INFO(40),HOLD(40)
      INTEGER ZCOL,RPOS
      DOUBLE PRECISION ZERO,HALF,ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
      INTEGER AINPUT
      DOUBLE PRECISION AMAX,AMULT1,AMULT2
      INTEGER APOS,APOSA,APOSB,APOSBB,APOSBK,APOSC,APOSI,APOSJ,APOSM,
     +        APOS1,APOS2,APOS3,APOS4,ASTK,ATRASH,BLK
      DOUBLE PRECISION DELTA,DETPIV
      INTEGER ELT
      DOUBLE PRECISION FLOPSA,FLOPSB,FLOPSX
      INTEGER I,I1,IASS,IBEG,IELL,IEND,IEXCH,IINPUT,INTSPA,
     +        IORG,IPIV,IPOS,IROW,ISNPIV,ISTK,ISWOP,IWNFS,IWPOS,
     +        J,JAY,JA1,JCOL,JJ,JJJ,JMAX,J1,J2,K,
     +        KB,KBLK,KCT,KR,KROW,K1,K2,L,LASPIV,LDIAG,LIELL,
     +        LP,LPIV, NBSTATIC
      LOGICAL LASTBK,LTWO
      INTEGER MAXFRT
      DOUBLE PRECISION MAXPIV
      INTEGER NASS,NBLK,NBLOC,NCMPBI,NCMPBR,NEIG,NELL,NFRONT,NIRBDU
      DOUBLE PRECISION NORMJ
      INTEGER NTWO
      LOGICAL SCHUR,LSTAT
      INTEGER MPIV,NPIV,NPOTPV,NRLBDU,NSC1,NST,
     +        NSTACK(2),NSTKAC(2),NTOTPV,
     +        NUMORG,OFFDAG,PHASE,PIVBLK
      DOUBLE PRECISION PIVOT
      INTEGER PIVSIZ,POSELT,POSPV1,POSPV2,PTRA,PTRIRN,RLSPA,
     +        SIZBLK,SIZC,SIZF,TRLSPA,TINSPA,TOTSTA(2),WP,ZCOUNT
      DOUBLE PRECISION RMAX,SWOP,TMAX,TOL,UU,ULOC,UTARG,STCTOL
      DOUBLE PRECISION FD15AD
C?? To identify bug
      INTRINSIC MIN,MAX,ABS
      EXTERNAL DGEMM,FD15AD,MA57PD,MA57WD
      NBLOC = ICNTL(11)
      TOL = CNTL(2)
      LP = ICNTL(1)
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      UU = MIN(CNTL(1),HALF)
      UU = MAX(UU,ZERO)
      LSTAT = .FALSE.
      IF (CNTL(4).GT.ZERO) THEN
        IF (CNTL(5).EQ.ZERO) LSTAT = .TRUE.
        UTARG = SQRT(UU/CNTL(4))*CNTL(4)
        STCTOL = BIGA*CNTL(4)
      ENDIF
      IF (HOLD(1).GT.0) THEN
        INFO(1) = 0
        NBLK = HOLD(2)
        NTWO = HOLD(3)
        INFO(23) = HOLD(4)
        NCMPBR = 0
        NCMPBI = 0
        NEIG   = HOLD(6)
        MAXFRT = HOLD(7)
        IWPOS  = HOLD(8)
        APOS   = HOLD(9)
        APOSBB = HOLD(10)
        NSTKAC(1) = HOLD(11)
        NSTKAC(2) = HOLD(12)
        AINPUT  = HOLD(13)
        IINPUT  = HOLD(14)
        ISTK    = HOLD(15)
        ASTK    = HOLD(16)
        INTSPA  = HOLD(17)
        RLSPA   = HOLD(18)
        PTRIRN  = HOLD(19)
        PTRA    = HOLD(20)
        NTOTPV  = HOLD(21)
        NPOTPV  = HOLD(22)
        NUMORG  = HOLD(23)
        NFRONT  = HOLD(24)
        NASS    = HOLD(25)
        IF (HOLD(1).EQ.1) NELL    = HOLD(27)
        IF (HOLD(1).EQ.2) NPIV    = HOLD(27)
        IASS    = HOLD(28)
        TINSPA  = HOLD(29)
        TRLSPA  = HOLD(30)
        TOTSTA(1) = HOLD(31)
        TOTSTA(2) = HOLD(32)
        NSTACK(1) = HOLD(33)
        NSTACK(2) = HOLD(34)
        INFO(32)  = HOLD(37)
        INFO(33)  = HOLD(38)
        INFO(34)  = HOLD(39)
        NBSTATIC  = HOLD(40)
        IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) ISNPIV = HOLD(35)
        IF (ICNTL(7).EQ.4) PHASE = HOLD(36)
        IF (HOLD(1).EQ.2) NSC1    = NFRONT-NPIV
        FLOPSA = RINFO(3)
        FLOPSB = RINFO(4)
        FLOPSX = RINFO(5)
        IF (HOLD(1).EQ.1) THEN
          HOLD(1) = 0
          GO TO 333
        ELSE
          IF (HOLD(1).EQ.3) THEN
            HOLD(1) = 0
            GO TO 555
          ELSE
            HOLD(1) = 0
            GO TO 444
          ENDIF
        ENDIF
      ENDIF
      NBSTATIC = 0
      NBLK = 0
      NTWO = 0
      NCMPBR = 0
      NCMPBI = 0
      FLOPSA = ZERO
      FLOPSB = ZERO
      FLOPSX = ZERO
      NEIG = 0
      MAXFRT  = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(17) = 0
      INFO(19) = 0
      INFO(40) = 0
      INFO(26) = 0
      INFO(27) = 0
      INFO(32) = 0
      INFO(33) = 0
      INFO(34) = 0
      INFO(23) = 0
      RINFO(3) = ZERO
      RINFO(4) = ZERO
      RINFO(5) = ZERO
      RINFO(15) = ZERO
      DO 10 I = 1,N
        PPOS(I) = N + 1
   10 CONTINUE
      IWPOS = 6
      IW(1) = 0
      IW(2) = 0
      IW(3) = 0
      IW(4) = 0
      IW(5) = 0
      APOSBB = 1
      NSTACK(1) = 0
      NSTACK(2) = 0
      NSTKAC(1) = NE
      NSTKAC(2) = NE
      TOTSTA(1) = NE
      TOTSTA(2) = NE
      INTSPA = NE+5+N
      RLSPA = NE
      TINSPA = NE+5+N
      TRLSPA = NE
      PTRIRN = LIW - NE + 1
      PTRA = LA - NE + 1
      ISTK = PTRIRN - 1
      ASTK = PTRA - 1
      AINPUT = PTRA
      IINPUT = PTRIRN
      NTOTPV = 0
      NPOTPV = 0
      IF (ICNTL(7).EQ.2 .OR. ICNTL(7).EQ.3) ISNPIV = 0
      IF (ICNTL(7).EQ.4) THEN
        PHASE = 1
        DO 19 I = 1,N
          DIAG(I) = ZERO
   19   CONTINUE
        APOS1 = PTRA-1
        J1 = PTRIRN
        DO 20 I = 1,N
          J2 = J1 + LROW(I) - 1
          DO 25 JJ = J1,J2
            J = IW(JJ)
            APOS1 = APOS1 + 1
            IF (J.EQ.PERM(I)) DIAG(J) = DIAG(J) + A(APOS1)
   25     CONTINUE
          J1 = J2 + 1
   20   CONTINUE
        SCHNAB(1) = ONE
        SCHNAB(5) = ZERO
        DO 21 I = 1,N
          SCHNAB(1) = MAX(SCHNAB(1),ABS(DIAG(I)))
          SCHNAB(5) = MIN(SCHNAB(5),DIAG(I))
   21   CONTINUE
        SCHNAB(4) = SCHNAB(1)
        SCHNAB(2) = FD15AD('E')**(1.0/3.0)
        SCHNAB(3) = 0.1
        RINFO(15) = FD15AD('H')
        DELTA     = ZERO
      ENDIF
      IASS = 1
 2160 CONTINUE
        NUMORG = 0
        DO 30 I = NPOTPV + 1,N
          J = PERM(I)
          IF (ABS(NODE(J)).GT.IASS) GO TO 40
          IW(IWPOS+NUMORG) = J
          NUMORG = NUMORG + 1
          PPOS(J) = NUMORG
   30   CONTINUE
   40   NASS = NUMORG
        NELL = NSTK(IASS)
        IELL = ISTK + 1
        DO 70 ELT = 1,NELL
          DO 50 JJ = IELL + 1,IELL + IW(IELL)
            J = IW(JJ)
            IF (NODE(J).GT.IASS) GO TO 50
            IF (PPOS(J).LE.N) GO TO 50
            IW(IWPOS+NASS) = J
            NASS = NASS + 1
            PPOS(J) = NASS
   50     CONTINUE
          IELL = IELL + IW(IELL) + 1
   70   CONTINUE
        IWNFS = IWPOS + NASS
        J1 = PTRIRN
        DO 90 IORG = 1,NUMORG
          J2 = J1 + LROW(NPOTPV+IORG) - 1
          DO 80 JJ = J1,J2
            J = IW(JJ)
            IF (PPOS(J).LE.N) GO TO 80
            IW(IWNFS) = J
            IWNFS = IWNFS + 1
            PPOS(J) = IWNFS - IWPOS
   80     CONTINUE
          J1 = J2 + 1
   90   CONTINUE
        IELL = ISTK + 1
        DO 170 ELT = 1,NELL
          J1 = IELL+1
          J2 = IELL+IW(IELL)
          DO 150 JJ = J1,J2
            J = IW(JJ)
            IF (PPOS(J).LE.N) GO TO 150
            IW(IWNFS) = J
            IWNFS = IWNFS + 1
            PPOS(J) = IWNFS - IWPOS
  150     CONTINUE
          IELL = J2 + 1
  170   CONTINUE
        NFRONT = IWNFS - IWPOS
        MAXFRT = MAX(MAXFRT,NFRONT)
        IF (INFO(1).NE.-3) THEN
          APOS = APOSBB + (NASS*(NASS+1))/2
        ELSE
          APOS = 1
        END IF
        RLSPA  = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+NSTKAC(1))
        TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+TOTSTA(1))
  333   IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
          CALL MA57PD(A,IW,ASTK,AINPUT,PTRA,.TRUE.)
          NCMPBR = NCMPBR + 1
          IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
            IF (ICNTL(8).NE.0) THEN
              DO 334 I = APOSBB,ASTK
                A(I) = ZERO
  334         CONTINUE
              HOLD(1) = 1
              HOLD(2) = NBLK
              HOLD(3) = NTWO
              HOLD(4) = INFO(23)
              HOLD(5) = NCMPBI
              HOLD(6) = NEIG
              HOLD(7) = MAXFRT
              HOLD(8) = IWPOS
              HOLD(9) = APOS
              HOLD(10) = APOSBB
              HOLD(11) = NSTKAC(1)
              HOLD(12) = NSTKAC(2)
              HOLD(13) = AINPUT
              HOLD(14) = IINPUT
              HOLD(15) = ISTK
              HOLD(16) = ASTK
              HOLD(17) = INTSPA
              HOLD(18) = RLSPA
              HOLD(19) = PTRIRN
              HOLD(20) = PTRA
              HOLD(21) = NTOTPV
              HOLD(22) = NPOTPV
              HOLD(23) = NUMORG
              HOLD(24) = NFRONT
              HOLD(25) = NASS
              HOLD(27) = NELL
              HOLD(28) = IASS
              HOLD(29) = TINSPA
              HOLD(30) = TRLSPA
              HOLD(31) = TOTSTA(1)
              HOLD(32) = TOTSTA(2)
              HOLD(33) = NSTACK(1)
              HOLD(34) = NSTACK(2)
              IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
              IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
              HOLD(37) = INFO(32)
              HOLD(38) = INFO(33)
              HOLD(39) = INFO(34)
              RINFO(3) =FLOPSA
              RINFO(4) =FLOPSB
              RINFO(5) =FLOPSX
              HOLD(40) = NBSTATIC
              INFO(35) = HOLD(40)
              INFO(1) = 10
              RETURN
            ELSE
              INFO(40) = INFO(40) + APOS - 1
              APOS = 1
              APOSBB = 1
              INFO(1) = -3
              IF (NFRONT*NFRONT.GT.ASTK) THEN
                INFO(17) = MAX(INFO(17),RLSPA)
                IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
                INFO(2) = LA
                RETURN
              ENDIF
            ENDIF
          ENDIF
        END IF
        ATRASH = APOS + NFRONT*NFRONT - 1
        DO 210 JJ = APOS,ATRASH
          A(JJ) = ZERO
  210   CONTINUE
        J1 = PTRIRN
        DO 230 IORG = 1,NUMORG
          J = PERM(NPOTPV+IORG)
          APOSI = APOS + (PPOS(J)-1)*NFRONT - 1
          J2 = J1 + LROW(NPOTPV+IORG) - 1
          FLOPSA = FLOPSA + J2 - J1 + 1
          DO 220 JJ = J1,J2
            JAY = IW(JJ)
CCC
              APOS2 = APOSI + PPOS(JAY)
            A(APOS2) = A(APOS2) + A(PTRA)
            PTRA = PTRA + 1
  220     CONTINUE
          NSTKAC(1) = NSTKAC(1) - J2 + J1 - 1
          J1 = J2 + 1
  230   CONTINUE
        NSTKAC(2) = NSTKAC(2) - J1 + PTRIRN
        PTRIRN = J1
        NPOTPV = NPOTPV + NUMORG
C???
C?? Depends if we need lower triangle and whether all entries are
        DO 380 ELT = 1,NELL
          POSELT = ASTK + 1
          LIELL = IW(ISTK+1)
          J1 = ISTK + 2
          J2 = ISTK+1 + LIELL
          FLOPSA = FLOPSA + (LIELL*(LIELL+1))/2
          DO 250 JJ = J1,J2
            J = IW(JJ)
            APOS2 = APOS + (PPOS(J)-1)*NFRONT
            APOS1 = POSELT
            DO 240 JJJ=JJ,J2
              JAY = IW(JJJ)
C???          APOS3 = APOS2 + PPOS(JAY) - 1
C???          A(APOS3) = A(APOS3) + A(APOS1)
C???          APOS5 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
C???          IF (APOS3.NE.APOS5) A(APOS5) = A(APOS5) + A(APOS1)
              IF (PPOS(JAY) .GE. PPOS(J)) THEN
                APOS3 = APOS2 + PPOS(JAY) - 1
              ELSE
                APOS3 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
              ENDIF
              A(APOS3) = A(APOS3) + A(APOS1)
              APOS1 = APOS1 + 1
  240       CONTINUE
            POSELT = POSELT + LIELL - (JJ-J1)
  250     CONTINUE
          NSTKAC(2) = NSTKAC(2) - (J2-ISTK)
          NSTACK(2) = NSTACK(2) - (J2-ISTK)
          TOTSTA(2) = TOTSTA(2) - (J2-ISTK)
          ISTK = J2
          ASTK = ASTK + (LIELL*(LIELL+1))/2
          NSTKAC(1) = NSTKAC(1) - (LIELL*(LIELL+1))/2
          NSTACK(1) = NSTACK(1) - (LIELL*(LIELL+1))/2
          TOTSTA(1) = TOTSTA(1) - (LIELL*(LIELL+1))/2
  380   CONTINUE
C1122     CONTINUE
        PIVBLK = MIN(NBLOC,NASS)
        APOSBK = APOS
        NPIV = 0
        ULOC = UU
        DO 918 BLK = 1,NASS
        IF (NPIV+PIVBLK .GE. NASS) THEN
          LASTBK = .TRUE.
          SIZBLK = NASS - NPIV
        ELSE
          LASTBK = .FALSE.
          SIZBLK = PIVBLK
        ENDIF
        LASPIV = NPIV
        MPIV = 0
        KR = 0
CCC Set to following to force 2 by 2 pivots in Nocedal examples
        KCT = SIZBLK + 1
  920   CONTINUE
          KR = KR + 1
          KCT = KCT - 1
          IF (KCT.EQ.0) GO TO 930
          IF (KR.GT.SIZBLK) KR = MPIV + 1
          IPIV = LASPIV + KR
            APOSI = APOS + (IPIV-1)*NFRONT
            POSPV1 = APOSI + IPIV - 1
            PIVOT = A(POSPV1)
   29       IF (ICNTL(7).EQ.4 .AND. PHASE.EQ.2) THEN
              IF (INFO(27).EQ.0) INFO(27) = NTOTPV + 1
              NORMJ = ZERO
              DO 28 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
                NORMJ = NORMJ + ABS(A(I))
   28         CONTINUE
              DELTA = MAX(ZERO,
     *                    - A(POSPV1) + MAX(NORMJ,SCHNAB(2)*SCHNAB(1)))
              A(POSPV1) = A(POSPV1) + DELTA
              IF (A(POSPV1).EQ.ZERO) GO TO 970
              RINFO(15) = MIN(RINFO(15),A(POSPV1))
              DIAG(PERM(NTOTPV+1)) = DELTA
              PIVSIZ = 1
              GO TO 811
            ENDIF
            IF (ICNTL(7).GT.1) THEN
              IF (ABS(PIVOT).LE.CNTL(2)) THEN
                IF (ICNTL(7).LT.4) GO TO 970
                PHASE = 2
                GO TO 29
              ENDIF
              IF (NTOTPV.EQ.0) THEN
                IF (PIVOT.GT.ZERO) ISNPIV = 1
                IF (PIVOT.LT.ZERO) ISNPIV = -1
              ELSE
                IF (ICNTL(7).EQ.2 .AND. ISNPIV*PIVOT.LT.ZERO) GO TO 980
                IF (ICNTL(7).EQ.3 .AND. ISNPIV*PIVOT.LT.ZERO) THEN
                    INFO(26) = INFO(26) + 1
                    ISNPIV = -ISNPIV
                ENDIF
              ENDIF
              IF (ICNTL(7).EQ.4) THEN
                IF (PIVOT.GE.SCHNAB(1)*SCHNAB(2) .AND.
     *              SCHNAB(5).GE.-SCHNAB(3)*SCHNAB(4)) THEN
                  SCHNAB(5) = ZERO
                  SCHNAB(4) = ZERO
                  DO 22 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
                    J = IW(IWPOS+NPIV+I-POSPV1)
                    DIAG(J) = DIAG(J) - A(I)*A(I)/PIVOT
                    SCHNAB(5) = MIN(DIAG(J),SCHNAB(5))
                    SCHNAB(4) = MAX(DIAG(J),SCHNAB(4))
                    IF (DIAG(J).LT.-SCHNAB(3)*SCHNAB(1)) THEN
                      PHASE = 2
                      GO TO 29
                    ENDIF
   22             CONTINUE
                  DIAG(PERM(NTOTPV+1)) = ZERO
                  RINFO(15) = MIN(RINFO(15),PIVOT)
                ELSE
                  PHASE = 2
                  GO TO 29
                ENDIF
              ENDIF
              PIVSIZ = 1
              GO TO 811
            ENDIF
            AMAX = ZERO
            JMAX = 0
            DO 110 K = 1, IPIV - NPIV - 1
              IF (ABS(A(POSPV1-K*NFRONT)).GT.AMAX) THEN
                AMAX = ABS(A(POSPV1-K*NFRONT))
                JMAX = IPIV - K
              ENDIF
  110       CONTINUE
            DO 111 K =  1, MIN(NASS,LASPIV+PIVBLK) - IPIV
              IF (ABS(A(POSPV1+K)).GT.AMAX) THEN
                AMAX = ABS(A(POSPV1+K))
                JMAX = IPIV + K
              ENDIF
  111       CONTINUE
            RMAX = ZERO
            DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NFRONT-IPIV
               RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
 112        CONTINUE
            IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
              GO TO 920
            END IF
            IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 920
            PIVSIZ = 0
            IF (ABS(PIVOT).GT.ULOC*MAX(RMAX,AMAX)) THEN
              PIVSIZ = 1
              A(POSPV1) = PIVOT
              GO TO 810
            END IF
            IF (NPIV+1.EQ.NASS) THEN
              A(POSPV1) = PIVOT
              GO TO 920
            END IF
            IF (AMAX.LE.TOL) GO TO 920
            IF (RMAX.LT.AMAX) THEN
              RMAX = ZERO
              DO 113 K = 1, IPIV - NPIV - 1
                IF (IPIV-K.EQ.JMAX) GO TO 113
                RMAX=MAX(RMAX,ABS(A(POSPV1-K*NFRONT)))
  113         CONTINUE
              DO 114 K =  1, NFRONT - IPIV
                IF (IPIV+K.EQ.JMAX) GO TO 114
                RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
  114         CONTINUE
            ENDIF
            APOSJ = APOS + (JMAX-1)*NFRONT
            POSPV2 = APOSJ + JMAX - 1
            IF (IPIV.GT.JMAX) THEN
              OFFDAG = APOSJ + IPIV - 1
            ELSE
              OFFDAG = APOSI + JMAX - 1
            END IF
            TMAX = ZERO
            DO 115 K = 1, JMAX - NPIV - 1
              IF (JMAX-K.EQ.IPIV) GO TO 115
              TMAX=MAX(TMAX,ABS(A(POSPV2-K*NFRONT)))
  115       CONTINUE
            DO 116 K =  1, NFRONT - JMAX
              IF (JMAX+K.EQ.IPIV) GO TO 116
              TMAX = MAX(TMAX,ABS(A(POSPV2+K)))
  116       CONTINUE
            DETPIV = A(POSPV1)*A(POSPV2) - AMAX*AMAX
            MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
            IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
            IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 920
            PIVSIZ = 2
            IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*ULOC.GT.
     +          ABS(DETPIV)) GO TO 920
            IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*ULOC.GT.
     +          ABS(DETPIV)) GO TO 920
  810       LPIV = IPIV
            IF (PIVSIZ.EQ.2) LPIV = MIN(IPIV,JMAX)
CCC         KR = MAX(KR,NPIV+PIVSIZ)
            KR = MAX(KR,MPIV+PIVSIZ)
            KCT = SIZBLK - MPIV - PIVSIZ + 1
            DO 860 KROW = NPIV,NPIV + PIVSIZ - 1
              IF (LPIV.EQ.KROW+1) GO TO 850
              JA1 = APOS + (LPIV-1)
              J1 = APOS + KROW
              DO 820 JJ = 1,KROW
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
                JA1 = JA1 + NFRONT
                J1 = J1 + NFRONT
  820         CONTINUE
              JA1 = JA1 + NFRONT
              J1 = J1 + 1
              DO 830 JJ = 1,LPIV - KROW - 2
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
                JA1 = JA1 + NFRONT
                J1 = J1 + 1
  830         CONTINUE
              SWOP = A(APOS+KROW* (NFRONT+1))
              A(APOS+KROW* (NFRONT+1)) = A(JA1)
              A(JA1) = SWOP
              DO 840 JJ = 1,NFRONT - LPIV
                JA1 = JA1 + 1
                J1 = J1 + 1
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
  840         CONTINUE
              IPOS = IWPOS + KROW
              IEXCH = IWPOS + LPIV - 1
              ISWOP = IW(IPOS)
              IW(IPOS) = IW(IEXCH)
              IW(IEXCH) = ISWOP
  850         LPIV = MAX(IPIV,JMAX)
  860       CONTINUE
  811       POSPV1 = APOS + NPIV* (NFRONT+1)
            POSPV2 = POSPV1 + NFRONT + 1
            IF (PIVSIZ.EQ.1) THEN
              FLOPSB = FLOPSB + ONE
              A(POSPV1) = ONE/A(POSPV1)
              IF (A(POSPV1).LT.ZERO) NEIG = NEIG + 1
              J1 = POSPV1 + 1
              J2 = POSPV1 + NASS - (NPIV+1)
              IBEG = POSPV1 + NFRONT + 1
              IEND = APOS + (NPIV+1)*NFRONT + NFRONT - 1
              DO 880 JJ = J1,J2
                AMULT1 = -A(JJ)*A(POSPV1)
                IF (.NOT.LASTBK) A(POSPV1+(JJ-J1+1)*NFRONT) = A(JJ)
                JCOL = JJ
                FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
                IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 871
CDIR$            IVDEP
                DO 870 IROW = IBEG,IEND
                  A(IROW) = A(IROW) + AMULT1*A(JCOL)
                  JCOL = JCOL + 1
  870           CONTINUE
  871           A(JJ) = AMULT1
                IBEG = IBEG + NFRONT + 1
                IEND = IEND + NFRONT
  880         CONTINUE
              NPIV = NPIV + 1
              MPIV = MPIV + 1
              NTOTPV = NTOTPV + 1
              IF (MPIV.EQ.SIZBLK) GO TO 930
            ELSE
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6.0
              SWOP = A(POSPV2)
              IF (DETPIV.LT.ZERO) THEN
                NEIG = NEIG + 1
              ELSE
                IF (SWOP.LT.ZERO) NEIG = NEIG + 2
              END IF
              A(POSPV2) = A(POSPV1)/DETPIV
              A(POSPV1) = SWOP/DETPIV
              A(OFFDAG) = -A(OFFDAG)/DETPIV
              J1 = POSPV1 + 2
              J2 = POSPV1 + NASS - (NPIV+1)
              IBEG = POSPV2 + NFRONT + 1
              IEND = APOS + (NPIV+2)*NFRONT + NFRONT - 1
              DO 900 JJ = J1,J2
                K1 = JJ
                K2 = JJ + NFRONT
                AMULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
                AMULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
                IF (.NOT.LASTBK) THEN
                  A(POSPV1 + (JJ-J1+2)*NFRONT) = A(K1)
                  A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = A(K2)
                ENDIF
                FLOPSB = FLOPSB + (IEND-IBEG+1)*4 + 6
                IF (MPIV+JJ-J1+3.GT.PIVBLK) GO TO 891
CDIR$            IVDEP
                DO 890 IROW = IBEG,IEND
                  A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
                  K1 = K1 + 1
                  K2 = K2 + 1
  890           CONTINUE
  891           A(JJ) = AMULT1
                A(JJ+NFRONT) = AMULT2
                IBEG = IBEG + NFRONT + 1
                IEND = IEND + NFRONT
  900         CONTINUE
              IPOS = IWPOS + NPIV
              IW(IPOS) = -IW(IPOS)
              IW(IPOS+1) = -IW(IPOS+1)
              NPIV = NPIV + 2
              MPIV = MPIV + 2
              NTOTPV = NTOTPV + 2
              NTWO = NTWO + 1
              IF (MPIV.EQ.SIZBLK) GO TO 930
            END IF
        GO TO 920
 930    IF (LASTBK) THEN
          IF (NPIV.EQ.NASS) GO TO 935
          IF (.NOT. LSTAT)  GO TO 935
          ULOC = ULOC/10.0D0
          IF (ULOC.LT.UTARG) THEN
            ULOC = ULOC * 10.0D0
            GO TO 9919
          ENDIF
          KCT = SIZBLK + 1 - MPIV
          GO TO 920
        ENDIF
        IF (MPIV.EQ.0) THEN
          PIVBLK = 2*PIVBLK
          GO TO 918
        ENDIF
        KBLK = (NASS-(LASPIV+PIVBLK))/PIVBLK
        L = NASS - (LASPIV+PIVBLK)
        APOS4 = APOS+(LASPIV+PIVBLK)*(NFRONT+1)
        DO 931 KB = 1,KBLK
          FLOPSX = FLOPSX + PIVBLK*(PIVBLK-1)*MPIV
          CALL DGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,ONE,
     +               A(APOSBK+PIVBLK*KB),NFRONT,
     +               A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,ONE,
     +               A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
          IF (NFRONT.GT.NASS)
     +    CALL DGEMM('N','T',NFRONT-NASS,PIVBLK,MPIV,ONE,
     +               A(APOSBK+NASS-LASPIV),NFRONT,
     +               A(APOSBK+PIVBLK*KB),NFRONT,ONE,
     +               A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)
  931   CONTINUE
       SIZC = NASS - (KBLK+1)*PIVBLK - LASPIV
       SIZF = NFRONT - (KBLK+1)*PIVBLK - LASPIV
       APOSA = APOSBK + (KBLK+1)*PIVBLK
       DO 934 K = 1,MPIV
         APOSB = APOSBK + NFRONT*PIVBLK*(KBLK+1) + K - 1
         APOSM = APOSBK + PIVBLK*(KBLK+1) + (K-1)*NFRONT
         APOSC = APOSBK + PIVBLK*(KBLK+1)*(NFRONT+1)
         DO 933 JJ = 1,SIZC
            DO 932 J = JJ,SIZC
              A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSB)
  932       CONTINUE
            DO 936 J = SIZC+1,SIZF
              A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSM)
  936       CONTINUE
            APOSC = APOSC + NFRONT
            APOSB = APOSB + NFRONT
            APOSM = APOSM + 1
  933     CONTINUE
          APOSA = APOSA + NFRONT
  934   CONTINUE
        APOSBK = APOSBK + MPIV*(NFRONT+1)
        LASPIV = NPIV
  918   CONTINUE
CCC
 9919      IPIV = LASPIV+MPIV
 9920      IPIV = IPIV + 1
CADD Probably not needed .. use only IPIV
           APOSI = APOS + (IPIV-1)*NFRONT
           POSPV1 = APOSI + IPIV - 1
           PIVOT = A(POSPV1)
CADD
CCC        PIVSIZ = 1
CCC        LPIV = IPIV
CCC        AMAX = ZERO
CCC        DO 9876 K = 1, IPIV - NPIV - 1
CCC          AMAX = MAX(AMAX,ABS(A(POSPV1-K*NFRONT)))
CCC76      CONTINUE
CCC        DO 9878 K =  1, NFRONT - IPIV
CCC          AMAX = MAX(AMAX,ABS(A(POSPV1+K)))
CCC78      CONTINUE
           IF (ABS(A(POSPV1)).LT.STCTOL) THEN
               PIVOT = STCTOL
              IF (A(POSPV1) .LT. ZERO) THEN
                 A(POSPV1) = -PIVOT
                 PIVOT     = -PIVOT
              ELSE
                 A(POSPV1) = PIVOT
              ENDIF
              NBSTATIC = NBSTATIC + 1
           ENDIF
           FLOPSB = FLOPSB + ONE
           A(POSPV1) = ONE/A(POSPV1)
           IF (A(POSPV1).LT.ZERO) NEIG = NEIG + 1
           J1 = POSPV1 + 1
           J2 = POSPV1 + NASS - (NPIV+1)
           IBEG = POSPV1 + NFRONT + 1
           IEND = APOSI + 2*NFRONT - 1
           DO 9880 JJ = J1,J2
              AMULT1 = -A(JJ)*A(POSPV1)
              JCOL = JJ
              FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
CDIR$            IVDEP
              DO 9870 IROW = IBEG,IEND
                 A(IROW) = A(IROW) + AMULT1*A(JCOL)
                 JCOL = JCOL + 1
 9870         CONTINUE
 9871         A(JJ) = AMULT1
              IBEG = IBEG + NFRONT + 1
              IEND = IEND + NFRONT
 9880      CONTINUE
           NPIV = NPIV + 1
           MPIV = MPIV + 1
           NTOTPV = NTOTPV + 1
           IF (MPIV.LT.SIZBLK) GO TO 9920
C********************************
C********************************
  935   SCHUR = (NBLOC.LT.(NFRONT-NASS) .AND. NPIV.GE.NBLOC)
        IF (ICNTL(16).EQ.1) THEN
          ZCOUNT = 0
          APOS4 = APOS + NPIV*NFRONT + NPIV
          APOSB = APOS4 + NFRONT
          APOSC = APOS4 + 1
          DO 4444 I = 2,NASS-NPIV
            DO 4443 J = 1,I-1
              A(APOSB) = A(APOSC)
              APOSB = APOSB + 1
              APOSC = APOSC + NFRONT
 4443       CONTINUE
            APOSB = APOS4 + NFRONT*I
            APOSC = APOS4 + I
 4444     CONTINUE
          I = NASS - NPIV
 4445     CONTINUE
          IF (ZCOUNT.EQ.I) GO TO 4450
          APOSB = APOS4 + (I-1)*NFRONT
          DO 4446 J = 1,NFRONT-NPIV
            IF (ABS(A(APOSB+J-1)).GT.TOL) GO TO 4449
 4446     CONTINUE
          ZCOUNT = ZCOUNT + 1
          DO 4447 J = 1,NFRONT-NPIV
            A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4447     CONTINUE
          DO 4448 J = 1,NFRONT-NPIV
            A(APOS4+NFRONT*(ZCOUNT-1)+J-1) = ZERO
 4448     CONTINUE
          ISWOP = IW(IWPOS+NPIV+ZCOUNT-1)
          IW(IWPOS+NPIV+ZCOUNT-1) = IW(IWPOS+NPIV+I-1)
          IW(IWPOS+NPIV+I-1) = ISWOP
          GO TO 4445
 4449     I = I - 1
          GO TO 4445
 4450     CONTINUE
        ELSE
          ZCOUNT = 0
        ENDIF
        NSC1 = NFRONT - NPIV - ZCOUNT
        IF (IASS.NE.NSTEPS) INFO(23) = INFO(23) + NASS - NPIV
        IF (CNTL(4).GT.ZERO .AND. INFO(23).GT.CNTL(5)*N) LSTAT = .TRUE.
        IF (NSC1.EQ.0) GO TO 1830
        IF (.NOT.SCHUR) THEN
          RLSPA = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
     +                      NSTKAC(1))
          TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
     +                      TOTSTA(1))
          NSTKAC(1) = NSTKAC(1) + ((NSC1+1)*NSC1)/2
          NSTACK(1) = NSTACK(1) + ((NSC1+1)*NSC1)/2
          TOTSTA(1) = TOTSTA(1) + ((NSC1+1)*NSC1)/2
          APOSI = APOS + NFRONT*NFRONT - 1
          DO 1370 JJ = 1,NFRONT-NPIV
            J = APOSI
            DO 1360 JJJ = 1,JJ
                A(ASTK) = A(J)
                ASTK = ASTK - 1
                J = J - 1
 1360       CONTINUE
            APOSI = APOSI - NFRONT
 1370     CONTINUE
          APOS4 = ASTK + 1
          J1 = IWPOS
          LTWO = .FALSE.
          POSPV1 = APOS
          DO 1450 I1 = 1,NPIV
            IF (LTWO) GO TO 1440
            APOSI = APOS + (I1-1)*NFRONT + NASS
            J2 = APOS + NFRONT* (I1-1) + NFRONT - 1
CCC What happens here ??
            APOSC = APOS4 +
     *       ((NASS-NPIV-ZCOUNT)*(2*NFRONT-NPIV-ZCOUNT-NASS+1))/2
            IF (IW(J1).GT.0) THEN
              FLOPSB = FLOPSB + (NFRONT-NASS) +
     *                          (NFRONT-NASS)* (NFRONT-NASS+1)
              DO 1410 JJ = APOSI,J2
                AMULT1 = -A(JJ)*A(POSPV1)
                DO 1400 JJJ = JJ,J2
                  A(APOSC) = A(APOSC) + AMULT1*A(JJJ)
                  APOSC = APOSC + 1
 1400           CONTINUE
                A(JJ) = AMULT1
 1410         CONTINUE
              J1 = J1 + 1
            ELSE
              POSPV2 = POSPV1 + NFRONT + 1
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6* (NFRONT-NASS) +
     +                 2* (NFRONT-NASS)* (NFRONT-NASS+1)
              DO 1430 JJ = APOSI,J2
                AMULT1 = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
                AMULT2 = -A(POSPV2)*A(JJ+NFRONT) - A(OFFDAG)*A(JJ)
                DO 1420 JJJ = JJ,J2
                  A(APOSC) = A(APOSC) + AMULT1*A(JJJ) +
     +                       AMULT2*A(JJJ+NFRONT)
                  APOSC = APOSC + 1
 1420           CONTINUE
                A(JJ) = AMULT1
                A(JJ+NFRONT) = AMULT2
 1430         CONTINUE
              J1 = J1 + 2
              POSPV1 = POSPV2
              LTWO = .TRUE.
              GO TO 1450
            END IF
 1440       LTWO = .FALSE.
            POSPV1 = POSPV1 + NFRONT + 1
 1450     CONTINUE
        ELSE
          APOS4 = APOS+NASS*(NFRONT+1)
        APOS3 = APOS+NASS*NFRONT
        J1 = IWPOS
        LTWO = .FALSE.
        POSPV1 = APOS
          DO 1490 I = 1,NPIV
            IF (LTWO) GO TO 1480
            APOSI = APOS + (I-1)*NFRONT + NASS
            POSELT = APOS3 + I - 1
            IF (IW(J1).GT.0) THEN
              FLOPSB = FLOPSB + (NFRONT-NASS)
              DO 1460 JJ = APOSI,APOS + NFRONT*I - 1
                A(POSELT) = A(JJ)
                A(JJ) = -A(JJ)*A(POSPV1)
                POSELT = POSELT + NFRONT
 1460         CONTINUE
              J1 = J1 + 1
            ELSE
              POSPV2 = POSPV1 + NFRONT + 1
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6* (NFRONT-NASS)
              DO 1470 JJ = APOSI,APOS + NFRONT*I - 1
                A(POSELT) = A(JJ)
                A(POSELT+1) = A(JJ+NFRONT)
                A(JJ) = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
                A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
     +                         A(OFFDAG)*A(POSELT)
                POSELT = POSELT + NFRONT
 1470         CONTINUE
              J1 = J1 + 2
              POSPV1 = POSPV2
              LTWO = .TRUE.
              GO TO 1490
            END IF
 1480       LTWO = .FALSE.
            POSPV1 = POSPV1 + NFRONT + 1
 1490     CONTINUE
          FLOPSB = FLOPSB + NPIV* (NFRONT-NASS)**2 +
     *                      NPIV* (NFRONT-NASS)
          KBLK = ( NFRONT-NASS)/NBLOC
          L =  NFRONT - NASS
          DO 1500 KB = 1,KBLK
            FLOPSX = FLOPSX + NBLOC* (NBLOC-1)* (NPIV)
            CALL DGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,ONE,
     +                 A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
     +                 A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,ONE,
     +                 A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
 1500     CONTINUE
          DO 1550 I = 1 + KBLK*NBLOC,L
            APOSA = APOS + NASS
            APOSB = APOS3 +(I-1)*NFRONT
            APOSC = APOS4 + (I-1)*NFRONT - 1
            DO 1540 K = 1,NPIV
              DO 1530 J = I,L
                A(APOSC+J) = A(APOSC+J) + A(APOSA+J-1)*A(APOSB)
 1530         CONTINUE
              APOSA = APOSA + NFRONT
              APOSB = APOSB + 1
 1540       CONTINUE
 1550     CONTINUE
          JA1 = APOS+NFRONT*NFRONT-1
          NSTKAC(1) = NSTKAC(1) + ((NSC1+1)* (NSC1))/2
          NSTACK(1) = NSTACK(1) + ((NSC1+1)* (NSC1))/2
          TOTSTA(1) = TOTSTA(1) + ((NSC1+1)* (NSC1))/2
          DO 1710 I = NSC1,1,-1
            DO 1700 JJ = JA1,JA1-(NSC1-I),-1
              A(ASTK) = A(JJ)
              ASTK = ASTK - 1
 1700       CONTINUE
            JA1 = JA1 - NFRONT
 1710     CONTINUE
        END IF
        NSTKAC(2) = NSTKAC(2) + NSC1 + 1
        NSTACK(2) = NSTACK(2) + NSC1 + 1
        TOTSTA(2) = TOTSTA(2) + NSC1 + 1
 1830   IF (IASS.EQ.NSTEPS) THEN
          INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+NSTKAC(2))
          TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+TOTSTA(2))
          GO TO 2158
        ELSE
          INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+NSTKAC(2))
          TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+TOTSTA(2))
        ENDIF
  444   NST = 0
        IF (NSC1.GT.0) NST = NSC1 + 1
        IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
          CALL MA57PD(A,IW,ISTK,IINPUT,PTRIRN,.FALSE.)
          NCMPBI = NCMPBI + 1
          IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
            IF (ICNTL(8).NE.0) THEN
              HOLD(1) = 2
              HOLD(2) = NBLK
              HOLD(3) = NTWO
              HOLD(4) = INFO(23)
              HOLD(5) = NCMPBI
              HOLD(6) = NEIG
              HOLD(7) = MAXFRT
              HOLD(8) = IWPOS
              HOLD(9) = APOS
              HOLD(10) = APOSBB
              HOLD(11) = NSTKAC(1)
              HOLD(12) = NSTKAC(2)
              HOLD(13) = AINPUT
              HOLD(14) = IINPUT
              HOLD(15) = ISTK
              HOLD(16) = ASTK
              HOLD(17) = INTSPA
              HOLD(18) = RLSPA
              HOLD(19) = PTRIRN
              HOLD(20) = PTRA
              HOLD(21) = NTOTPV
              HOLD(22) = NPOTPV
              HOLD(23) = NUMORG
              HOLD(24) = NFRONT
              HOLD(25) = NASS
              HOLD(27) = NPIV
              HOLD(28) = IASS
              HOLD(29) = TINSPA
              HOLD(30) = TRLSPA
              HOLD(31) = TOTSTA(1)
              HOLD(32) = TOTSTA(2)
              HOLD(33) = NSTACK(1)
              HOLD(34) = NSTACK(2)
              IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
              IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
              HOLD(37) = INFO(32)
              HOLD(38) = INFO(33)
              HOLD(39) = INFO(34)
              NSC1    = NFRONT-NPIV
              RINFO(3) =FLOPSA
              RINFO(4) =FLOPSB
              RINFO(5) =FLOPSX
              INFO(1) = 11
              HOLD(40) = NBSTATIC
              INFO(35) = HOLD(40)
            ELSE
              INFO(1)  = -4
              INFO(2)  = LIW
              INFO(18) = INTSPA
            ENDIF
            RETURN
          END IF
        END IF
        IF (NSC1.GT.0) THEN
          DO 1720 I = 1,NSC1
            IW(ISTK) = IW(IWPOS+NFRONT-I)
            ISTK = ISTK - 1
 1720     CONTINUE
          IW(ISTK) = NSC1
          ISTK = ISTK - 1
        ENDIF
        DO 1840 JJ = IWPOS + NPIV,IWPOS + NFRONT - 1
          J = ABS(IW(JJ))
          PPOS(J) = N + 1
 1840   CONTINUE
C********************************
C********************************
 2158   IF (NPIV.EQ.0) GO TO 2159
        NBLK = NBLK + 1
        IW(IWPOS-2) = NFRONT
        IW(IWPOS-1) = NPIV
        IWPOS = IWPOS + NFRONT + 2
        IF (INFO(1).EQ.-3) THEN
          INFO(40) = INFO(40) + (NPIV * (2*NFRONT-NPIV+1))/2
          GO TO 2159
        END IF
        APOS2 = APOSBB
        DO 2130 I = 1,NPIV
          JA1 = APOS + (I-1)* (NFRONT+1)
          DO 2120 J = I,NPIV
            A(APOS2) = A(JA1)
            IF (A(APOS2).EQ.ZERO) INFO(32) = INFO(32) + 1
            APOS2 = APOS2 + 1
            JA1 = JA1 + 1
 2120     CONTINUE
 2130   CONTINUE
        RPOS = APOS2
        DO 2150 I = 1,NPIV
          JA1 = APOS + (I-1)*NFRONT + NPIV
          DO 2140 J = 1,NFRONT - NPIV
            A(APOS2) = A(JA1)
            APOS2 = APOS2 + 1
            JA1 = JA1 + 1
 2140     CONTINUE
 2150   CONTINUE
        APOSBB = APOS2
        DO 2152 J = 1,NFRONT-NPIV
        APOS2 = RPOS+J-1
        ZCOL = 1
          DO 2151 I = 1,NPIV
            IF (A(APOS2).EQ.ZERO) INFO(33) = INFO(33)+1
            IF (A(APOS2).NE.ZERO) ZCOL = 0
            APOS2 = APOS2 + NFRONT - NPIV
 2151     CONTINUE
        IF (ZCOL.EQ.1) INFO(34) = INFO(34)+1
 2152   CONTINUE
 2159   IASS = IASS + 1
      IF (IASS.LE.NSTEPS) THEN
        IW(IWPOS-2) = 0
        IW(IWPOS-1) = 0
        GO TO 2160
      ENDIF
C2160 CONTINUE
      INFO(35) = NBSTATIC
      IF (INFO(1).EQ.-3) THEN
        INFO(2)  = LA
        INFO(17) = MAX(INFO(17),RLSPA)
        IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
        RETURN
      END IF
      GO TO 1000
 970  INFO(1) = -5
      INFO(2) = NTOTPV + 1
      RINFO(20) = PIVOT
      IF (LDIAG.GT.0 .AND. LP.GE.0)
     *    WRITE(LP,99992) INFO(1),PIVOT,CNTL(2),INFO(2),ICNTL(7)
99992 FORMAT (/'*** Error message from routine MA57BD **',
     *       '   INFO(1) = ',I3/'Pivot has value ',D16.8,' when ',
     *       'CNTL(2) has value ',D16.8/
     *       'at stage',I11,2X,'when ICNTL(7) =',I3)
      RETURN
 980  INFO(1) = -6
      INFO(2) = NTOTPV + 1
      RINFO(20) = PIVOT
      IF (LDIAG.GT.0 .AND. LP.GE.0)
     *    WRITE(LP,99993) INFO(1),INFO(2),ICNTL(7)
99993 FORMAT (/'*** Error message from routine MA57BD **',
     *       '   INFO(1) = ',I3/'Change in sign of pivot at stage',
     *       I10,2X,'when ICNTL(7) = ',I3)
      RETURN
 1000 NRLBDU = APOSBB - 1
      NIRBDU = IWPOS - 3
      IF (NTOTPV.NE.N) THEN
        INFO(1) = 4
        RINFO(20) = PIVOT
        IF (LDIAG.GT.0 .AND. WP.GE.0)
     *      WRITE(WP,99994) INFO(1),NTOTPV
99994 FORMAT (/'*** Warning message from routine MA57BD **',
     *         '   INFO(1) =',I2/5X, 'Matrix is singular, rank =', I5)
      ENDIF
  555 IF (NTOTPV.NE.N .AND. ICNTL(16).EQ.1) THEN
        IF (NIRBDU+3*(N-NTOTPV) .GT. LIW
     +    .OR. NRLBDU+(N-NTOTPV)+NTWO .GT. LA) THEN
          IF (ICNTL(8).NE.0) THEN
            HOLD(1) = 3
            HOLD(2) = NBLK
            HOLD(3) = NTWO
            HOLD(4) = INFO(23)
            HOLD(5) = NCMPBI
            HOLD(6) = NEIG
            HOLD(7) = MAXFRT
            HOLD(8) = IWPOS
            HOLD(9) = APOS
            HOLD(10) = APOSBB
            HOLD(11) = NSTKAC(1)
            HOLD(12) = NSTKAC(2)
            HOLD(13) = AINPUT
            HOLD(14) = IINPUT
            HOLD(15) = ISTK
            HOLD(16) = ASTK
            HOLD(17) = INTSPA
            HOLD(18) = RLSPA
            HOLD(19) = PTRIRN
            HOLD(20) = PTRA
            HOLD(21) = NTOTPV
            HOLD(22) = NPOTPV
            HOLD(23) = NUMORG
            HOLD(24) = NFRONT
            HOLD(25) = NASS
            HOLD(27) = NPIV
            HOLD(28) = IASS
            HOLD(29) = TINSPA
            HOLD(30) = TRLSPA
            HOLD(31) = TOTSTA(1)
            HOLD(32) = TOTSTA(2)
            HOLD(33) = NSTACK(1)
            HOLD(34) = NSTACK(2)
            IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
            IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
            HOLD(37) = INFO(32)
            HOLD(38) = INFO(33)
            HOLD(39) = INFO(34)
            NSC1    = NFRONT-NPIV
            RINFO(3) =FLOPSA
            RINFO(4) =FLOPSB
            RINFO(5) =FLOPSX
            INFO(1) = 11
            HOLD(40) = NBSTATIC
            INFO(35) = HOLD(40)
          ELSE
            IF (NIRBDU+3*(N-NTOTPV) .GT. LIW) THEN
              INFO(1)  = -4
              INFO(2)  = LIW
              INFO(18) = MAX(INTSPA,NIRBDU+3*(N-NTOTPV))
            ELSE
              INFO(1)  = -3
              INFO(2) = LA
              INFO(17) = MAX(INFO(17),RLSPA,NRLBDU+(N-NTOTPV)+NTWO)
              IF (ICNTL(7).EQ.4) INFO(17) =
     +          MAX(INFO(17),RLSPA + N,NRLBDU+(N-NTOTPV)+NTWO)
            ENDIF
          ENDIF
          RETURN
        ENDIF
      ENDIF
      IF (N.NE.NTOTPV) THEN
        DO 3331 I = 1,N
          PPOS(I) = 0
 3331   CONTINUE
        IWPOS = 4
        DO 3332 I = 1,NBLK
          NFRONT = IW(IWPOS)
          NPIV = IW(IWPOS+1)
          DO 3330 J = IWPOS+2,IWPOS+NPIV+1
            PPOS(ABS(IW(J))) = 1
 3330     CONTINUE
          IWPOS = IWPOS + NFRONT + 2
 3332   CONTINUE
        K= 0
        DO 3333 I=1,N
          IF (PPOS(I).EQ.0) THEN
            K=K+1
            NBLK = NBLK + 1
            NRLBDU = NRLBDU+1
            A(NRLBDU) = ONE
            IW(NIRBDU+1) = 1
            IW(NIRBDU+2) = 1
            IW(NIRBDU+3) = I
            NIRBDU = NIRBDU+3
          ENDIF
 3333   CONTINUE
      ENDIF
      INFO(14) = NRLBDU
      IW(1) = NRLBDU + 1
      IW(2) = NRLBDU + NTWO
      INFO(15) = IW(2)
      IW(3) = NBLK
      INFO(31) = NBLK
      CALL MA57WD(A,LA,IW,LIW,NRLBDU)
      INFO(16) = NIRBDU
      INFO(18) = INTSPA
      INFO(20) = TINSPA
      INFO(17) = RLSPA
      INFO(19) = TRLSPA
      INFO(21) = MAXFRT
      INFO(22) = NTWO
      INFO(24) = NEIG
      INFO(25) = NTOTPV
      INFO(28) = NCMPBR
      INFO(29) = NCMPBI
      RINFO(3) = FLOPSA
      RINFO(4) = FLOPSB
      RINFO(5) = FLOPSX
      IF (INFO(27).GT.0) THEN
        RINFO(14) = ZERO
        DO 332 I = 1,N
          RINFO(14) = MAX(RINFO(14),DIAG(I))
 332    CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE MA57PD(A,IW,J1,J2,ITOP,REAL)
      INTEGER ITOP,J1,J2
      LOGICAL REAL
      DOUBLE PRECISION A(*)
      INTEGER IW(*)
      INTEGER IPOS,JJ
      IF (J2.EQ.ITOP) GO TO 50
      IPOS = ITOP - 1
      IF (REAL) THEN
        DO 10 JJ = J2-1,J1+1,-1
          A(IPOS) = A(JJ)
          IPOS = IPOS - 1
   10   CONTINUE
      ELSE
        DO 20 JJ = J2-1,J1+1,-1
          IW(IPOS) = IW(JJ)
          IPOS = IPOS - 1
   20   CONTINUE
      ENDIF
      J2 = ITOP
      J1 = IPOS
   50 RETURN
      END
      SUBROUTINE MA57WD(A,LA,IW,LIW,NRLBDU)
      INTEGER LA,LIW
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)
      INTEGER NRLBDU
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER APOS,IBLK,IROW,IWPOS,J,JPIV,NCOLS,NROWS
      APOS = 1
      IWPOS = 6
      DO 40 IBLK = 1,IW(3)
        NCOLS = IW(IWPOS-2)
        NROWS = IW(IWPOS-1)
        JPIV = 1
        DO 30 IROW = 1,NROWS
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 10
          IF (IW(IWPOS+IROW-1).LT.0) THEN
            JPIV = 2
            NRLBDU = NRLBDU + 1
            A(NRLBDU) = A(APOS+1)
            A(APOS+1) = ZERO
          END IF
   10     DO 20 J = APOS + 1,APOS + NROWS - IROW
            A(J) = -A(J)
   20     CONTINUE
          APOS = APOS + NROWS - IROW + 1
   30   CONTINUE
        APOS = APOS + NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + NCOLS + 2
   40 CONTINUE
      END
      SUBROUTINE MA57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,K1,K2,
     +        NCOLS,NROWS
      DOUBLE PRECISION W1,W2
      APOS = 1
      IWPOS = 4
      DO 270 IBLK = 1,IFACT(3)
        IW1(IBLK) = IWPOS
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = RHS(II)
   10     CONTINUE
          CALL DTPSV('L','N','U',NROWS,FACT(APOS),W,1)
          APOS = APOS + (NROWS* (NROWS+1))/2
          IF (NCOLS.GT.NROWS) CALL DGEMV('N',NCOLS-NROWS,NROWS,
     +                                  ONE,FACT(APOS),NCOLS-NROWS,
     +                                  W,1,ONE,W(NROWS+1),1)
          APOS = APOS + NROWS* (NCOLS-NROWS)
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   35     CONTINUE
        ELSE
        J1 = IWPOS
        J2 = IWPOS + NROWS - 1
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          W1 = RHS(ABS(IFACT(J1)))
          K = APOS
          DO 100 J = J1+1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) - FACT(K)*W1
            K = K + 1
  100     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS-1,2
          K1 = APOS
          K2 = APOS+NCOLS-NROWS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          W2 = RHS(ABS(IFACT(IWPOS+IPIV)))
          DO 133 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K1) + W2*FACT(K2)
            K1 = K1 + 1
            K2 = K2 + 1
  133     CONTINUE
          APOS = K2
  136   CONTINUE
        IF (MOD(NROWS,2).EQ.1) THEN
          K = APOS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          DO 137 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
            K = K + 1
  137     CONTINUE
          APOS = K
        ENDIF
      END IF
      IWPOS = IWPOS + NCOLS
  270 CONTINUE
      END
      SUBROUTINE MA57YD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,J1,J2,K,K2,LROW,NCOLS,NROWS
      DOUBLE PRECISION W1,W2
      APOS = IFACT(1)
      APOS2 = IFACT(2)
      DO 380 IBLK = IFACT(3),1,-1
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 5 I = NROWS + 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = RHS(II)
    5     CONTINUE
          DO 10 IPIV = NROWS,1,-1
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            APOS = APOS - (NROWS+1-IPIV)
            W(IPIV) = RHS(IRHS)*FACT(APOS)
   10     CONTINUE
          JPIV = -1
          DO 20 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              W(IPIV) = RHS(IRHS1)*FACT(APOS2) + W(IPIV)
              IF (JPIV.EQ.1) APOS2 = APOS2 - 1
              JPIV = -JPIV
            END IF
   20     CONTINUE
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMV('T',K,NROWS,ONE,
     +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
     +                           W(NROWS+1),1,ONE,W,1)
          CALL DTPSV('L','T','U',NROWS,FACT(APOS),W,1)
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   60     CONTINUE
        ELSE
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1
          JPIV = -1
          DO 210 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            LROW = NROWS + 1 - IPIV
            IF (IRHS.GT.0) THEN
              APOS = APOS - LROW
              RHS(IRHS) = RHS(IRHS)*FACT(APOS)
            ELSE
              IF (JPIV.EQ.-1) THEN
                IRHS1 = -IFACT(IWPOS+IPIV-2)
                IRHS2 = -IRHS
                APOS = APOS - LROW - LROW - 1
                W1 = RHS(IRHS1)*FACT(APOS) +
     +               RHS(IRHS2)*FACT(APOS2)
                RHS(IRHS2) = RHS(IRHS1)*FACT(APOS2) +
     +                         RHS(IRHS2)*FACT(APOS+LROW+1)
                RHS(IRHS1) = W1
                APOS2 = APOS2 - 1
              END IF
              JPIV = -JPIV
            END IF
  210     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2
          K = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS-1,2
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W1 = RHS(IRHS)
            IRHS1 = ABS(IFACT(IWPOS+IPIV))
            W2 = RHS(IRHS1)
            K2 = K+(NCOLS-NROWS)
            DO 215 J = J1,J2
              II = ABS(IFACT(J))
              W1 = W1 + FACT(K)*RHS(II)
              W2 = W2 + FACT(K2)*RHS(II)
              K = K + 1
              K2 = K2 + 1
  215       CONTINUE
            RHS(IRHS) = W1
            RHS(IRHS1) = W2
            K = K2
  220     CONTINUE
          IF (MOD(NROWS,2).EQ.1) THEN
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W1 = RHS(IRHS)
            DO 216 J = J1,J2
              W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
              K = K + 1
  216       CONTINUE
            RHS(IRHS) = W1
          ENDIF
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            W1 = RHS(IRHS)
            K = APOS + 1
            DO 230 J = J1,J2
              W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)))
              K = K + 1
  230       CONTINUE
            RHS(IRHS) = W1
            J1 = J1 - 1
  260     CONTINUE
        END IF
  380 CONTINUE
      END
      SUBROUTINE MA57VD(N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
      INTEGER IWFR,LW,N,NZ
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
      INTEGER ICNTL(*),INFO(*)
      INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LR,N1,NDUP
      INFO(2) = 0
      DO 10 I = 1,N
        IPE(I) = 0
   10 CONTINUE
      LR = NZ
      IF (NZ.EQ.0) GO TO 120
      DO 110 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 90
        ELSE IF (I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 90
        ELSE
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF
   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA57AD',
     +          '  *** INFO(1) =',I2)
        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF
   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')
   80   I = 0
        J = 0
        GO TO 100
   90   IPE(I) = IPE(I) + 1
        IPE(J) = IPE(J) + 1
  100   IW(K) = J
        LR = LR + 1
        IW(LR) = I
  110 CONTINUE
  120 IQ(1) = 1
      N1 = N - 1
      IF (N1.LE.0) GO TO 140
      DO 130 I = 1,N1
        FLAG(I) = 0
        IF (IPE(I).EQ.0) IPE(I) = -1
        IQ(I+1) = IPE(I) + IQ(I) + 1
        IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
      FLAG(N) = 0
      IF (LR.GE.LAST) GO TO 160
      K1 = LR + 1
      DO 150 K = K1,LAST
        IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
      IWFR = LAST + 1
      IF (NZ.EQ.0) GO TO 230
      DO 220 K = 1,NZ
        J = IW(K)
        IF (J.LE.0) GO TO 220
        L = K
        IW(K) = 0
        DO 210 ID = 1,NZ
          IF (L.GT.NZ) GO TO 170
          L = L + NZ
          GO TO 180
  170     L = L - NZ
  180     I = IW(L)
          IW(L) = 0
          IF (I.LT.J) GO TO 190
          L = IQ(J) + 1
          IQ(J) = L
          JN = IW(L)
          IW(L) = -I
          GO TO 200
  190     L = IQ(I) + 1
          IQ(I) = L
          JN = IW(L)
          IW(L) = -J
  200     J = JN
          IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
  230 NDUP = 0
      DO 280 I = 1,N
        K1 = IPE(I) + 1
        K2 = IQ(I)
        IF (K1.LE.K2) GO TO 240
        IPE(I) = 0
        IQ(I) = 0
        GO TO 280
  240   DO 260 K = K1,K2
          J = -IW(K)
          IF (J.LE.0) GO TO 270
          L = IQ(J) + 1
          IQ(J) = L
          IW(L) = I
          IW(K) = J
          IF (FLAG(J).NE.I) GO TO 250
          NDUP = NDUP + 1
          IW(L) = 0
          IW(K) = 0
  250     FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
        IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
      IF (NDUP.EQ.0) GO TO 310
      IWFR = 1
      DO 300 I = 1,N
        K1 = IPE(I) + 1
        IF (K1.EQ.1) GO TO 300
        K2 = IQ(I) + IPE(I)
        L = IWFR
        IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 290 K = K1,K2
          IF (IW(K).EQ.0) GO TO 290
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
  290   CONTINUE
        IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN
      END
      SUBROUTINE MA57HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
     +                 NCMPA,FRATIO)
      DOUBLE PRECISION FRATIO
      INTEGER IWFR,LW,N,IOVFLO,NCMPA
      INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
      INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
     +        KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
     +        NP0,NS,NVPIV,NVROOT,ROOT
      EXTERNAL MA57ZD
      INTRINSIC ABS,MIN
      DO 10 I = 1,N
        IPD(I) = 0
        NV(I) = 1
        FLAG(I) = IOVFLO
   10 CONTINUE
      MD = 1
      NCMPA = 0
      NFLG = IOVFLO
      NEL = 0
      ROOT = N+1
      NVROOT = 0
      DO 30 IS = 1,N
        K = IPE(IS)
        IF (K.LE.0) GO TO 20
        ID = IW(K) + 1
        NS = IPD(ID)
        IF (NS.GT.0) LST(NS) = IS
        NXT(IS) = NS
        IPD(ID) = IS
        LST(IS) = -ID
        GO TO 30
   20   NEL = NEL + 1
        FLAG(IS) = -1
        NXT(IS) = 0
        LST(IS) = 0
   30 CONTINUE
      DO 340 ML = 1,N
        IF (NEL+NVROOT+1.GE.N) GO TO 350
        DO 40 ID = MD,N
          MS = IPD(ID)
          IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
        NVPIV = NV(MS)
        NS = NXT(MS)
        NXT(MS) = 0
        LST(MS) = 0
        IF (NS.GT.0) LST(NS) = -ID
        IPD(ID) = NS
        ME = MS
        NEL = NEL + NVPIV
        IDN = 0
        KP = IPE(ME)
        FLAG(MS) = -1
        IP = IWFR
        LEN = IW(KP)
        DO 140 KP1 = 1,LEN
          KP = KP + 1
          KE = IW(KP)
          IF (FLAG(KE).LE.-2) GO TO 60
          IF (FLAG(KE).LE.0) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 140
             KE = ROOT
             IF (FLAG(KE).LE.0) GO TO 140
          END IF
          JP = KP - 1
          LN = LEN - KP1 + 1
          IE = MS
          GO TO 70
   60     IE = KE
          JP = IPE(IE)
          LN = IW(JP)
   70     DO 130 JP1 = 1,LN
            JP = JP + 1
            IS = IW(JP)
            IF (FLAG(IS).LE.0) THEN
               IF (IPE(IS).EQ.-ROOT) THEN
                  IS = ROOT
                  IW(JP) = ROOT
                  IF (FLAG(IS).LE.0) GO TO 130
               ELSE
                  GO TO 130
               END IF
            END IF
            FLAG(IS) = 0
            IF (IWFR.LT.LW) GO TO 100
            IPE(MS) = KP
            IW(KP) = LEN - KP1
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA57ZD(N,IPE,IW,IP-1,LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 90
            DO 80 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   80       CONTINUE
   90       IP = LWFR
            JP = IPE(IE)
            KP = IPE(ME)
  100       IW(IWFR) = IS
            IDN = IDN + NV(IS)
            IWFR = IWFR + 1
            LS = LST(IS)
            LST(IS) = 0
            NS = NXT(IS)
            NXT(IS) = 0
            IF (NS.GT.0) LST(NS) = LS
            IF (LS.LT.0) THEN
              LS = -LS
              IPD(LS) = NS
            ELSE IF (LS.GT.0) THEN
              NXT(LS) = NS
            END IF
  130     CONTINUE
          IF (IE.EQ.MS) GO TO 150
          IPE(IE) = -ME
          FLAG(IE) = -1
  140   CONTINUE
  150   NV(MS) = IDN + NVPIV
        IF (IWFR.EQ.IP) GO TO 330
        K1 = IP
        K2 = IWFR - 1
        LIMIT = NINT(FRATIO*(N-NEL))
        DO 310 K = K1,K2
          IS = IW(K)
          IF (IS.EQ.ROOT) GO TO 310
          IF (NFLG.GT.2) GO TO 170
          DO 160 I = 1,N
            IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
            IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160     CONTINUE
          NFLG = IOVFLO
  170     NFLG = NFLG - 1
          ID = IDN
          KP1 = IPE(IS) + 1
          NP = KP1
          KP2 = IW(KP1-1) + KP1 - 1
          DO 220 KP = KP1,KP2
            KE = IW(KP)
          IF (FLAG(KE).EQ.-1) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 220
             KE = ROOT
             IW(KP) = ROOT
             IF (FLAG(KE).EQ.-1) GO TO 220
          END IF
          IF (FLAG(KE).GE.0) GO TO 230
            JP1 = IPE(KE) + 1
            JP2 = IW(JP1-1) + JP1 - 1
            IDL = ID
            DO 190 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).LE.NFLG) GO TO 190
              ID = ID + NV(JS)
              FLAG(JS) = NFLG
  190       CONTINUE
            IF (ID.GT.IDL) GO TO 210
            DO 200 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).NE.0) GO TO 210
  200       CONTINUE
            IPE(KE) = -ME
            FLAG(KE) = -1
            GO TO 220
  210       IW(NP) = KE
            FLAG(KE) = -NFLG
            NP = NP + 1
  220     CONTINUE
          NP0 = NP
          GO TO 250
  230     KP0 = KP
          NP0 = NP
          DO 240 KP = KP0,KP2
            KS = IW(KP)
            IF (FLAG(KS).LE.NFLG) THEN
               IF (IPE(KS).EQ.-ROOT) THEN
                  KS = ROOT
                  IW(KP) = ROOT
                  IF (FLAG(KS).LE.NFLG) GO TO 240
               ELSE
                  GO TO 240
               END IF
            END IF
            ID = ID + NV(KS)
            FLAG(KS) = NFLG
            IW(NP) = KS
            NP = NP + 1
  240     CONTINUE
  250     IF (ID.GE.LIMIT) GO TO 295
          IW(NP) = IW(NP0)
          IW(NP0) = IW(KP1)
          IW(KP1) = ME
          IW(KP1-1) = NP - KP1 + 1
          JS = IPD(ID)
          DO 280 L = 1,N
            IF (JS.LE.0) GO TO 300
            KP1 = IPE(JS) + 1
            IF (IW(KP1).NE.ME) GO TO 300
            KP2 = KP1 - 1 + IW(KP1-1)
            DO 260 KP = KP1,KP2
              IE = IW(KP)
              IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260       CONTINUE
            GO TO 290
  270       JS = NXT(JS)
  280     CONTINUE
  290     IPE(JS) = -IS
          NV(IS) = NV(IS) + NV(JS)
          NV(JS) = 0
          FLAG(JS) = -1
          NS = NXT(JS)
          LS = LST(JS)
          IF (NS.GT.0) LST(NS) = IS
          IF (LS.GT.0) NXT(LS) = IS
          LST(IS) = LS
          NXT(IS) = NS
          LST(JS) = 0
          NXT(JS) = 0
          IF (IPD(ID).EQ.JS) IPD(ID) = IS
          GO TO 310
  295     IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IW(K) = ROOT
            IPE(IS) = -ROOT
            NV(ROOT) = NV(ROOT) + NV(IS)
            NV(IS) = 0
            FLAG(IS) = -1
          END IF
          NVROOT = NV(ROOT)
          GO TO 310
  300     NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
          MD = MIN(MD,ID)
  310   CONTINUE
        DO 320 K = K1,K2
          IS = IW(K)
          IF (NV(IS).EQ.0) GO TO 320
          FLAG(IS) = NFLG
          IW(IP) = IS
          IP = IP + 1
  320   CONTINUE
        IWFR = K1
        FLAG(ME) = -NFLG
        IW(IP) = IW(K1)
        IW(K1) = IP - K1
        IPE(ME) = K1
        IWFR = IP + 1
        GO TO 335
  330   IPE(ME) = 0
  335   CONTINUE
  340 CONTINUE
  350 DO 360 IS = 1,N
        IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
          IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IPE(IS) = -ROOT
          END IF
          NVROOT = NVROOT + NV(IS)
          NV(IS) = 0
         END IF
  360 CONTINUE
      DO 370 IE = 1,N
        IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
      IF(NVROOT.GT.0)NV(ROOT)=NVROOT
      END
      SUBROUTINE MA57ZD(N,IPE,IW,LW,IWFR,NCMPA)
      INTEGER IWFR,LW,N,NCMPA
      INTEGER IPE(N),IW(LW)
      INTEGER I,IR,K,K1,K2,LWFR
      NCMPA = NCMPA + 1
      DO 10 I = 1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR = 1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K = LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K = K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
C *******************************************************************
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
C           Council for the Central Laboratory of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.
C 23 May 2007 Version 1.1.0. Absolute value of hash taken to cover the
C            case of integer overflow.
C            Comments with character in column 2 corrected.
C 2 August 2007 Version 2.0.0 Dense row handling added, error & warning
C            messages added, iovflo added, interface changed. MC47I/ID
C            added.
C 31 October 2007 Version 2.1.0 Corrected tree formation when handling
C            full variables
C

      SUBROUTINE MC47ID(ICNTL)
      INTEGER ICNTL(10)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 1
      ICNTL(5) = 2139062143
      DO 100 I=6,10
        ICNTL(I) = 0
 100  CONTINUE
      RETURN
      END
      SUBROUTINE MC47AD(N, NE, PE, IW, IWLEN,
     *      ICNTL,INFO, RINFO)
      INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), INFO(10)
      INTEGER ICNTL(10)
      DOUBLE PRECISION RINFO(10)
      INTEGER DEGREE
      DOUBLE PRECISION DUMMY(1)
      INTEGER ELEN,HEAD,I,II,I1,I2,J,LAST,LEN,LENIW,LP,MP,
     *        NEXT,NV,PFREE,W,WP
      INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,JNFO(10)
      EXTERNAL MC59AD,MC34AD,MC47BD
      DO 5 J = 1,10
         INFO(J) = 0
 5    CONTINUE
      LP = ICNTL(1)
      WP = ICNTL(2)
      MP = ICNTL(3)
      IF (N.LT.1) THEN
        INFO(1) = -1
        IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +       '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +       'N has value ',N
        GO TO 1000
      ENDIF
      IF (PE(1).LT.1) THEN
        IF (2*NE+N.GT.IWLEN) THEN
          INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN
          GO TO 1000
        ENDIF
      ELSE
        IF (NE+N.GT.IWLEN) THEN
          INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN
          GO TO 1000
        ENDIF
      ENDIF
      IF (MP.GE.0) THEN
        WRITE(MP,'(/A)') 'Entry to MC47A/AD'
        WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
     *                            ' entries'
        IF (PE(1).LT.0)  THEN
          WRITE(MP,'(A)') 'Matrix input in coordinate form'
          WRITE(MP,'(A/(4(I8,I8)))') 'Row and column indices',
     *          (IW(I),IW(NE+I),I=1,NE)
        ELSE
          WRITE(MP,'(A)') 'Matrix input by columns'
          DO 10 J=1,N
            WRITE(MP,'(A,I4/(10I8))') 'Column',J,
     *                                (IW(I),I=PE(J),PE(J+1)-1)
   10     CONTINUE
        ENDIF
      ENDIF
      LAST   = IWLEN  - N + 1
      ELEN   = LAST   - N
      NV     = ELEN   - N
      W      = NV     - N
      DEGREE = W      - N
      HEAD   = DEGREE - N
      NEXT   = HEAD   - N
      LEN    = NEXT   - N
      LENIW = LEN-1
      INFO(6) = 0
      INFO(7) = 0
      IF (PE(1).LT.0) THEN
        DO 20 I=1,NE
          IF (IW(I).LE.IW(NE+I)) THEN
            IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
              INFO(7) = INFO(7) + 1
            ELSE
              IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
            ENDIF
            IW(I)=0
          ENDIF
   20   CONTINUE
        ICT59(1) = 0
        ICT59(2) = 1
        ICT59(3) = 1
        ICT59(4) = LP
        ICT59(5) = -1
        ICT59(6) = 0
        CALL MC59AD(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
     *              N+1,PE,N+1,IW(2*NE+1),INFO59)
        IDUP  = INFO59(3)
        IOUT  = INFO59(4)
        JOUT  = INFO59(5)
      ELSE
        IDUP = 0
        IOUT = 0
        JOUT = 0
        DO 30 I = 1,N
          IW(NE+I) = 0
   30   CONTINUE
        DO 50 J=1,N
          I1 = PE(J)
          PE(J) = I1-(IOUT+IDUP)
          I2 = PE(J+1)-1
          IF (I2.LT.I1-1) THEN
            INFO(1) = -3
            GO TO 1000
          ENDIF
          DO 40 II = I1,I2
            I = IW(II)
            IF (I.LE.J .OR. I.GT.N) THEN
              IF (I.EQ.J) INFO(7) = INFO(7) + 1
              IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
              IOUT = IOUT + 1
            ELSE
              IF (IW(NE+I).EQ.J) THEN
                IDUP = IDUP + 1
              ELSE
                IW(NE+I)=J
                IW(II-(IOUT+IDUP)) = I
              ENDIF
            ENDIF
   40     CONTINUE
   50   CONTINUE
        PE(N+1) = NE - (IOUT+IDUP) + 1
      ENDIF
      IF (IDUP.GT.0) THEN
        INFO(1) = 1
        INFO(4) = IDUP
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of duplicates found: ',INFO(4)
      ELSE
        INFO(4) = 0
      ENDIF
      IF (IOUT+ JOUT - INFO(7) .GT.0 ) THEN
        INFO(1) = 1
        INFO(5) = IOUT + JOUT - INFO(7)
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of out of range entries found and ignored: ',
     +       INFO(5)
      ELSE
        INFO(5) = 0
      ENDIF
      IF (INFO(6).GT.0) THEN
         INFO(1) = 1
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of entries in upper triangle found and ignored: ',
     +        INFO(6)
      ENDIF
      IF (INFO(7).GT.0) THEN
         INFO(1) = 1
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of entries in diagonals found and ignored: ',
     +        INFO(7)
      ENDIF
      IF (NE-(IOUT+IDUP).EQ.0) THEN
        INFO(1) = -4
        IF (LP.GE.0) WRITE(LP,'(/A,I3/A)')
     +       '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +       'Matrix is null'
        GO TO 1000
      ENDIF
      IF (LENIW.LT.2*(PE(N+1)-1)) THEN
        INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN,
     +        'Should be at least', 2*(PE(N+1)-1)+8*N
        GO TO 1000
      ENDIF
      CALL MC34AD(N,IW,PE,.FALSE.,DUMMY,IW(W))
      PFREE = PE(N+1)
      DO 60 I=1,N
        IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE
      CALL MC47BD(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),
     *            IW(ELEN),IW(LAST),IW(DEGREE),
     *            IW(HEAD),IW(NEXT),IW(W), ICNTL,JNFO, RINFO)
      INFO(2) = JNFO(1)
      INFO(3) = PFREE+8*N
      INFO(8) = JNFO(2)
      IF (MP.GE.0) THEN
        WRITE(MP,'(/A)') 'Exit from MC47A/AD'
        WRITE(MP,'(A/(7I10))') 'INFO(1-10):',(INFO(I),I=1,10)
        WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Inverse permutation',
     *                         (IW(LAST+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
      ENDIF
 1000 RETURN
      END
      SUBROUTINE MC47BD (N, IWLEN, PE, PFREE, LEN, IW, NV,
     $                   ELEN, LAST, DEGREE,
     $                   HEAD, DENXT, W, ICNTL, JNFO, RJNFO)
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
     $        ELEN(N), LAST(N),  DEGREE(N),
     $         HEAD(N), DENXT(N), W(N), ICNTL(10), JNFO(10)
      DOUBLE PRECISION RJNFO(10)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $     IDUMMY, ILAST, INEXT, IOVFLO,J, JDUMMY, JLAST, JNEXT, K,
     $     KNT1, KNT2, KNT3, LASTD,  LENJ, LN, MAXMEM, ME,
     $     MEM, MINDEG, NBD, NCMPA, NDME, NEL, NELME, NEWMEM,
     $     NFULL, NLEFT, NRLADU, NVI, NVJ, NVPIV, P, P1, P2, P3, PDST,
     $     PEE, PEE1, PEND, PJ, PME, PME1, PME2, PN, PSRC, RSTRT,
     $     SLENME, THRESH, THRESM, WE, WFLG, WNVI,X
     $
      DOUBLE PRECISION RELDEN, SM, STD, OPS
      LOGICAL IDENSE
      INTRINSIC MAX, MIN, MOD
      DO 2 I = 1,10
         RJNFO(I) = 0.0
         JNFO(I) = 0
 2    CONTINUE
      DMAX = 0
      HMOD = MAX (1, N-1)
      IOVFLO = ICNTL(5)
      LASTD = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 1
      NBD   = 0
      NCMPA = 0
      NEL = 0
      NFULL  = 0
      NRLADU = 0
      RSTRT = 0
      OPS = 0.00
      THRESH = ICNTL(4)
      WFLG = 2
      IF (THRESH.GT.0) THEN
         THRESM  = 0
         RELDEN = 0.0
         SM = 0
         DO 5 I=1,N
            THRESM = MAX(THRESM, LEN(I))
            IF (LEN(I).GT.0) THEN
               RELDEN = RELDEN + LEN(I)
               SM = SM + (LEN(I) * LEN(I))
            END IF
            LAST (I) = 0
            HEAD (I) = 0
            NV (I) = 1
            DEGREE (I) = LEN (I)
            IF (DEGREE(I) .EQ. 0) THEN
               NEL = NEL + 1
               ELEN (I) = -NEL
               PE (I) = 0
               W (I) = 0
               NRLADU = NRLADU + 1
               OPS = OPS + 1
            ELSE
               W (I) = 1
               ELEN (I) = 0
            ENDIF
 5       CONTINUE
         IF (N .EQ. NEL) GOTO 265
         RELDEN = RELDEN/(N-NEL)
         SM = SM/(N-NEL-NFULL) - RELDEN*RELDEN
         STD = SQRT(ABS(SM))
         IF (STD .LE. RELDEN) THEN
            THRESM = -1
         ELSE
            THRESM = INT(9*RELDEN + 0.5*STD*((STD/(RELDEN+0.01))**1.5)+
     *           2*RELDEN*RELDEN/(STD+0.01) +1)
         END IF
      ELSE
         THRESM = THRESH
         DO 10 I = 1, N
            LAST (I) = 0
            HEAD (I) = 0
            NV (I) = 1
            DEGREE (I) = LEN (I)
            IF (DEGREE(I) .EQ. 0) THEN
               NEL = NEL + 1
               ELEN (I) = -NEL
               PE (I) = 0
               W (I) = 0
               NRLADU = NRLADU + 1
               OPS = OPS + 1
            ELSE
               W (I) = 1
               ELEN (I) = 0
            ENDIF
 10      CONTINUE
      ENDIF
      IF (THRESM.GE.0) THEN
         IF (THRESM.GE.N) THEN
            THRESM = -1
         ELSE IF (THRESM.EQ.0) THEN
            THRESM = N
         ENDIF
      ENDIF
      DO 20 I = 1, N
         DEG = DEGREE (I)
         IF (DEG .GT. 0) THEN
            IF ( (THRESM.GE.0) .AND.
     &           (DEG+1.GE.THRESM.OR.DEG+1.GE.N-NEL )) THEN
               NBD = NBD+1
               IF (DEG+1.NE.N-NEL) THEN
                  DEGREE(I) = DEGREE(I)+N+1
                  DEG = N
                  INEXT = HEAD (DEG)
                  IF (INEXT .NE. 0) LAST (INEXT) = I
                  DENXT (I) = INEXT
                  HEAD (DEG) = I
                  LAST(I)  = 0
                  IF (LASTD.EQ.0) THEN
                     LASTD=I
                  END IF
               ELSE
                  NFULL = NFULL+1
                  DEGREE(I) = N+1
                  DEG = N
                  IF (LASTD.EQ.0) THEN
                     LASTD     = I
                     HEAD(DEG) = I
                     DENXT(I)   = 0
                     LAST(I)   = 0
                  ELSE
                        DENXT(LASTD) = I
                        LAST(I)     = LASTD
                        LASTD       = I
                        DENXT(I)     = 0
                  ENDIF
               ENDIF
            ELSE
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               DENXT (I) = INEXT
               HEAD (DEG) = I
            ENDIF
         ENDIF
 20   CONTINUE
      IF (NBD.EQ.0 .AND. THRESH.GT.0) THEN
         THRESM = -1
      END IF
 30   IF (NEL .LT. N) THEN
         DO 40 DEG = MINDEG, N
            ME = HEAD (DEG)
            IF (ME .GT. 0) GO TO 50
 40      CONTINUE
 50      MINDEG = DEG
         IF (DEG.LT.N)  THEN
            INEXT = DENXT (ME)
            IF (INEXT .NE. 0) LAST (INEXT) = 0
            HEAD (DEG) = INEXT
         ELSE
            IF (DEGREE(ME).EQ.N+1) GO TO 263
            RSTRT = RSTRT + 1
            RELDEN = 0.0
            SM = 0
            IF (WFLG .GT. IOVFLO-NBD-1) THEN
               DO  51 X = 1, N
                  IF (W (X) .NE. 0) W (X) = 1
 51            CONTINUE
               WFLG = 2
            END IF
            WFLG = WFLG + 1
            DO 57 IDUMMY = 1,N
               INEXT = DENXT (ME)
               IF (INEXT .NE. 0) THEN
                  LAST (INEXT) = 0
               ELSE
                  LASTD = 0
               ENDIF
               DENXT(ME) = 0
               W(ME)      = WFLG
               P1 = PE(ME)
               P2 = P1 + LEN(ME) -1
               LN       = P1
               ELN      = P1
               DO 55 P=P1,P2
                  E= IW(P)
                  IF (W(E).EQ.WFLG) GO TO 55
                  W(E) = WFLG
                  DO 52 JDUMMY = 1,N
                     IF ( PE(E) .GE. 0 ) GOTO 53
                     E = -PE(E)
                     IF (W(E) .EQ.WFLG) GOTO 55
                     W(E) = WFLG
 52               CONTINUE
 53               IF (ELEN(E).LT.0) THEN
                     DENXT(E) = DENXT(E) - NV(ME)
                     IW(LN) = IW(ELN)
                     IW(ELN) = E
                     LN  = LN+1
                     ELN = ELN + 1
                     PEE1 = PE(E)
                     DO 54 PEE = PEE1, PEE1+LEN(E)-1
                        X = IW(PEE)
                        IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                           DENXT(ME) = DENXT(ME) + NV(X)
                           W(X) = WFLG
                        ENDIF
 54                  CONTINUE
                  ELSE
                     DENXT(ME) = DENXT(ME) + NV(E)
                     IW(LN)=E
                     LN = LN+1
                  ENDIF
 55            CONTINUE
               WFLG     = WFLG + 1
               LEN(ME)  = LN-P1
               ELEN(ME) = ELN- P1
               NDME = DENXT(ME)+NV(ME)
               IF (DENXT(ME).EQ.0) DENXT(ME) =1
               IF (NDME .LT. NBD) THEN
                  RELDEN = RELDEN + NV(ME)*NDME
                  SM = SM + NV(ME)*NDME*NDME
                  DEGREE(ME) = DENXT(ME)
                  DEG = DEGREE(ME)
                  MINDEG = MIN(DEG,MINDEG)
                  JNEXT = HEAD(DEG)
                  IF (JNEXT.NE. 0) LAST (JNEXT) = ME
                  DENXT(ME) = JNEXT
                  HEAD(DEG) = ME
               ELSE
                  DEGREE(ME) = N+1
                  DEG = DENXT(ME)
                  MINDEG = MIN(DEG,MINDEG)
                  DEG = N
                  P1 = PE(ME)
                  P2 = P1 + ELEN(ME) - 1
                  DO 56 PJ=P1,P2
                     E= IW(PJ)
                     DENXT (E) = DENXT(E) + NV(ME)
 56               CONTINUE
                  DEG = N
                  NFULL = NFULL +NV(ME)
                  IF (LASTD.EQ.0) THEN
                     LASTD     = ME
                     HEAD(N) = ME
                     DENXT(ME)   = 0
                     LAST(ME)   = 0
                     IF (INEXT.EQ.0) INEXT = LASTD
                  ELSE
                        DENXT(LASTD) = ME
                        LAST(ME)     = LASTD
                        LASTD        = ME
                        DENXT(ME)     = 0
                        IF (INEXT.EQ.0) INEXT = LASTD
                  ENDIF
               END IF
               ME    = INEXT
               IF (ME.EQ.0) GO TO 58
               IF (DEGREE(ME).LE.(N+1) ) GOTO 58
 57         CONTINUE
 58         HEAD (N) = ME
            IF (NBD.EQ.NFULL) THEN
               RELDEN = 0
               SM = 0
            ELSE
               RELDEN = (RELDEN + NFULL*NBD)/(NBD)
               SM = (SM + NFULL*NBD*NBD)/(NBD) - RELDEN*RELDEN
            END IF
            STD = SQRT(ABS(SM))
            THRESM = INT(9*RELDEN+0.5*STD*((STD/(RELDEN + 0.01))**1.5)
     *           + 2*RELDEN*RELDEN/(STD+0.01) +1)
            THRESM = MIN(THRESM,NBD)
            IF (THRESM.GE.NBD) THEN
               THRESM = N
            END IF
            NBD = NFULL
            GOTO 30
         ENDIF
         ELENME = ELEN (ME)
         ELEN (ME) = - (NEL + 1)
         NVPIV = NV (ME)
         NEL = NEL + NVPIV
         DENXT(ME) = 0
         NV (ME) = -NVPIV
         DEGME = 0
         IF (ELENME .EQ. 0) THEN
            PME1 = PE (ME)
            PME2 = PME1 - 1
            DO 60 P = PME1, PME1 + LEN (ME) - 1
               I = IW (P)
               NVI = NV (I)
               IF (NVI .GT. 0) THEN
                  DEGME = DEGME + NVI
                  NV (I) = -NVI
                  PME2 = PME2 + 1
                  IW (PME2) = I
                  IF (DEGREE(I).LE.N) THEN
                     ILAST = LAST (I)
                     INEXT = DENXT (I)
                     IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                     IF (ILAST .NE. 0) THEN
                        DENXT (ILAST) = INEXT
                     ELSE
                        HEAD (DEGREE (I)) = INEXT
                     ENDIF
                  ELSE
                     DENXT(ME) = DENXT(ME) + NVI
                  ENDIF
               ENDIF
 60         CONTINUE
            NEWMEM = 0
         ELSE
            P  = PE (ME)
            PME1 = PFREE
            SLENME = LEN (ME) - ELENME
            DO 120 KNT1 = 1, ELENME
               E = IW (P)
               P = P + 1
               PJ = PE (E)
               LN = LEN (E)
               DO 110 KNT2 = 1, LN
                  I = IW (PJ)
                  PJ = PJ + 1
                  NVI = NV (I)
                  IF (NVI .GT. 0) THEN
                     IF (PFREE .GT. IWLEN) THEN
                        PE (ME) = P
                        LEN (ME) = LEN (ME) - KNT1
                        IF (LEN (ME) .EQ. 0) PE (ME) = 0
                        PE (E) = PJ
                        LEN (E) = LN - KNT2
                        IF (LEN (E) .EQ. 0) PE (E) = 0
                        NCMPA = NCMPA + 1
                        DO 70 J = 1, N
                           PN = PE (J)
                           IF (PN .GT. 0) THEN
                              PE (J) = IW (PN)
                              IW (PN) = -J
                           ENDIF
 70                     CONTINUE
                        PDST = 1
                        PSRC = 1
                        PEND = PME1 - 1
                        DO 91 IDUMMY = 1, IWLEN
                           IF (PSRC .GT. PEND) GO TO 95
                           J = -IW (PSRC)
                           PSRC = PSRC + 1
                           IF (J .GT. 0) THEN
                              IW (PDST) = PE (J)
                              PE (J) = PDST
                              PDST = PDST + 1
                              LENJ = LEN (J)
                              DO 90 KNT3 = 0, LENJ - 2
                                 IW (PDST + KNT3) = IW (PSRC + KNT3)
 90                           CONTINUE
                              PDST = PDST + LENJ - 1
                              PSRC = PSRC + LENJ - 1
                           ENDIF
 91                     END DO
 95                     P1 = PDST
                        DO 100 PSRC = PME1, PFREE - 1
                           IW (PDST) = IW (PSRC)
                           PDST = PDST + 1
 100                    CONTINUE
                        PME1 = P1
                        PFREE = PDST
                        PJ = PE (E)
                        P = PE (ME)
                     ENDIF
                     DEGME = DEGME + NVI
                     NV (I) = -NVI
                     IW (PFREE) = I
                     PFREE = PFREE + 1
                     IF (DEGREE(I).LE.N) THEN
                        ILAST = LAST (I)
                        INEXT = DENXT (I)
                        IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                        IF (ILAST .NE. 0) THEN
                           DENXT (ILAST) = INEXT
                        ELSE
                           HEAD (DEGREE (I)) = INEXT
                        ENDIF
                     ELSE
                        DENXT(ME) = DENXT(ME) + NVI
                     ENDIF
                  ENDIF
 110           CONTINUE
                  PE (E) = -ME
                  W (E) = 0
 120        CONTINUE
            KNT1 = ELENME + 1
            E = ME
            PJ = P
            LN = SLENME
            DO 126 KNT2 = 1, LN
               I = IW (PJ)
               PJ = PJ + 1
               NVI = NV (I)
               IF (NVI .GT. 0) THEN
                  IF (PFREE .GT. IWLEN) THEN
                     PE (ME) = P
                     LEN (ME) = LEN (ME) - KNT1
                     IF (LEN (ME) .EQ. 0) PE (ME) = 0
                     PE (E) = PJ
                     LEN (E) = LN - KNT2
                     IF (LEN (E) .EQ. 0) PE (E) = 0
                     NCMPA = NCMPA + 1
                     DO 121 J = 1, N
                        PN = PE (J)
                        IF (PN .GT. 0) THEN
                           PE (J) = IW (PN)
                           IW (PN) = -J
                        ENDIF
 121                 CONTINUE
                     PDST = 1
                     PSRC = 1
                     PEND = PME1 - 1
                     DO 123 IDUMMY = 1,IWLEN
                        IF (PSRC .GT. PEND) GO TO 124
                        J = -IW (PSRC)
                        PSRC = PSRC + 1
                        IF (J .GT. 0) THEN
                           IW (PDST) = PE (J)
                           PE (J) = PDST
                           PDST = PDST + 1
                           LENJ = LEN (J)
                           DO 122 KNT3 = 0, LENJ - 2
                              IW (PDST + KNT3) = IW (PSRC + KNT3)
 122                       CONTINUE
                           PDST = PDST + LENJ - 1
                           PSRC = PSRC + LENJ - 1
                        ENDIF
 123                 END DO
 124                 P1 = PDST
                     DO 125 PSRC = PME1, PFREE - 1
                        IW (PDST) = IW (PSRC)
                        PDST = PDST + 1
 125                 CONTINUE
                     PME1 = P1
                     PFREE = PDST
                     PJ = PE (E)
                     P = PE (ME)
                  END IF
                  DEGME = DEGME + NVI
                  NV (I) = -NVI
                  IW (PFREE) = I
                  PFREE = PFREE + 1
                  IF (DEGREE(I).LE.N) THEN
                     ILAST = LAST (I)
                     INEXT = DENXT (I)
                     IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                     IF (ILAST .NE. 0) THEN
                        DENXT (ILAST) = INEXT
                     ELSE
                        HEAD (DEGREE (I)) = INEXT
                     ENDIF
                  ELSE
                     DENXT(ME) = DENXT(ME) + NVI
                  ENDIF
               ENDIF
 126           CONTINUE
            PME2 = PFREE - 1
            NEWMEM = PFREE - PME1
            MEM = MEM + NEWMEM
            MAXMEM = MAX (MAXMEM, MEM)
         ENDIF
         DEGREE (ME) = DEGME
         PE (ME) = PME1
         LEN (ME) = PME2 - PME1 + 1
         IF (WFLG .GT. IOVFLO-N) THEN
            DO 130 X = 1, N
               IF (W (X) .NE. 0) W (X) = 1
 130        CONTINUE
            WFLG = 2
         ENDIF
         IF (NBD.GT.0) THEN
            DO 150 PME = PME1, PME2
               I = IW (PME)
               IF (DEGREE(I).GT.N) GOTO 150
               ELN = ELEN (I)
               IF (ELN .GT. 0) THEN
                  NVI = -NV (I)
                  WNVI = WFLG - NVI
                  DO 140 P = PE (I), PE (I) + ELN - 1
                     E = IW (P)
                     WE = W (E)
                     IF (WE .GE. WFLG) THEN
                        WE = WE - NVI
                     ELSE IF (WE .NE. 0) THEN
                        WE = DEGREE (E) + WNVI - DENXT(E)
                     ENDIF
                     W (E) = WE
 140              CONTINUE
               ENDIF
 150        CONTINUE
         ELSE
            DO 152 PME = PME1, PME2
               I = IW (PME)
               ELN = ELEN (I)
               IF (ELN .GT. 0) THEN
                  NVI = -NV (I)
                  WNVI = WFLG - NVI
                  DO 151 P = PE (I), PE (I) + ELN - 1
                     E = IW (P)
                     WE = W (E)
                     IF (WE .GE. WFLG) THEN
                        WE = WE - NVI
                     ELSE IF (WE .NE. 0) THEN
                        WE = DEGREE (E) + WNVI
                     ENDIF
                     W (E) = WE
 151              CONTINUE
               ENDIF
 152        CONTINUE
         END IF
         IF (NBD.GT.0) THEN
            DO 180 PME = PME1, PME2
               I = IW (PME)
               IF (DEGREE(I).GT.N) GOTO 180
               P1 = PE (I)
               P2 = P1 + ELEN (I) - 1
               PN = P1
               HASH = 0
               DEG = 0
               DO 160 P = P1, P2
                  E = IW (P)
                  DEXT = W (E) - WFLG
                  IF (DEXT .GT. 0) THEN
                     DEG = DEG + DEXT
                     IW (PN) = E
                     PN = PN + 1
                     HASH = HASH+E
                  ELSE IF ((DEXT .EQ. 0) .AND.
     &                    (DENXT(ME).EQ.NBD)) THEN
                     PE (E) = -ME
                     W (E)  = 0
                  ELSE IF (DEXT.EQ.0) THEN
                     IW(PN) = E
                     PN     = PN+1
                     HASH = HASH + E
                  ENDIF
 160           CONTINUE
               ELEN (I) = PN - P1 + 1
               P3 = PN
               DO 170 P = P2 + 1, P1 + LEN (I) - 1
                  J = IW (P)
                  NVJ = NV (J)
                  IF (NVJ .GT. 0) THEN
                     IF (DEGREE(J).LE.N) DEG=DEG+NVJ
                     IW (PN) = J
                     PN = PN + 1
                     HASH = HASH + J
                  ENDIF
 170           CONTINUE
               IF ((DEG .EQ. 0).AND.(DENXT(ME).EQ.NBD)) THEN
                  PE (I) = -ME
                  NVI = -NV (I)
                  DEGME = DEGME - NVI
                  NVPIV = NVPIV + NVI
                  NEL = NEL + NVI
                  NV (I) = 0
                  ELEN (I) = 0
               ELSE
                  DEGREE(I) = MIN (DEG+NBD-DENXT(ME), DEGREE(I))
                  IW (PN) = IW (P3)
                  IW (P3) = IW (P1)
                  IW (P1) = ME
                  LEN (I) = PN - P1 + 1
                  HASH = ABS(MOD (HASH, HMOD)) + 1
                  J = HEAD (HASH)
                  IF (J .LE. 0) THEN
                     DENXT (I) = -J
                     HEAD (HASH) = -I
                  ELSE
                     DENXT (I) = LAST (J)
                     LAST (J) = I
                  ENDIF
                  LAST (I) = HASH
               ENDIF
 180        CONTINUE
         ELSE
            DO 183 PME = PME1, PME2
               I = IW (PME)
               P1 = PE (I)
               P2 = P1 + ELEN (I) - 1
               PN = P1
               HASH = 0
               DEG = 0
               DO 181 P = P1, P2
                  E = IW (P)
                  DEXT = W (E) - WFLG
                  IF (DEXT .GT. 0) THEN
                     DEG = DEG + DEXT
                     IW (PN) = E
                     PN = PN + 1
                     HASH = HASH + E
                  ELSE IF (DEXT .EQ. 0) THEN
                     PE (E) = -ME
                     W (E)  = 0
                  ENDIF
 181           CONTINUE
               ELEN (I) = PN - P1 + 1
               P3 = PN
               DO 182 P = P2 + 1, P1 + LEN (I) - 1
                  J = IW (P)
                  NVJ = NV (J)
                  IF (NVJ .GT. 0) THEN
                     DEG=DEG+NVJ
                     IW (PN) = J
                     PN = PN + 1
                     HASH = HASH + J
                  ENDIF
 182           CONTINUE
               IF (DEG .EQ. 0) THEN
                  PE (I) = -ME
                  NVI = -NV (I)
                  DEGME = DEGME - NVI
                  NVPIV = NVPIV + NVI
                  NEL = NEL + NVI
                  NV (I) = 0
                  ELEN (I) = 0
               ELSE
                  DEGREE(I) = MIN (DEG,  DEGREE(I))
                  IW (PN) = IW (P3)
                  IW (P3) = IW (P1)
                  IW (P1) = ME
                  LEN (I) = PN - P1 + 1
                  HASH = ABS(MOD (HASH, HMOD)) + 1
                  J = HEAD (HASH)
                  IF (J .LE. 0) THEN
                     DENXT (I) = -J
                     HEAD (HASH) = -I
                  ELSE
                     DENXT (I) = LAST (J)
                     LAST (J) = I
                  ENDIF
                  LAST (I) = HASH
               ENDIF
 183        CONTINUE
         END IF
         DEGREE (ME) = DEGME
         DMAX = MAX (DMAX, DEGME)
         WFLG = WFLG + DMAX
         IF (WFLG .GE. IOVFLO - N) THEN
            DO 190 X = 1, N
               IF (W (X) .NE. 0) W (X) = 1
 190        CONTINUE
            WFLG = 2
         ENDIF
         DO 250 PME = PME1, PME2
            I = IW (PME)
            IF ( (NV(I).GE.0) .OR. (DEGREE(I).GT.N) ) GO TO 250
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
               I = -J
               HEAD (HASH) = 0
            ELSE
               I = LAST (J)
               LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
            DO 247 JDUMMY = 1,N
               IF (DENXT (I) .EQ. 0) GO TO 250
               LN = LEN (I)
               ELN = ELEN (I)
               DO 210 P = PE (I) + 1, PE (I) + LN - 1
                  W (IW (P)) = WFLG
 210           CONTINUE
               JLAST = I
               J = DENXT (I)
               DO 245 IDUMMY=1,N
                  IF (J .EQ. 0) GO TO 246
                  IF (LEN (J) .NE. LN) GO TO 240
                  IF (ELEN (J) .NE. ELN) GO TO 240
                  DO 230 P = PE (J) + 1, PE (J) + LN - 1
                     IF (W (IW (P)) .NE. WFLG) GO TO 240
 230              CONTINUE
                  PE (J) = -I
                  NV (I) = NV (I) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                  J = DENXT (J)
                  DENXT (JLAST) = J
                  GO TO 245
 240              CONTINUE
                  JLAST = J
                  J = DENXT (J)
 245           CONTINUE
 246           WFLG = WFLG + 1
               I = DENXT (I)
               IF (I .EQ. 0) GO TO 250
 247         CONTINUE
 250      CONTINUE
          P = PME1
          NLEFT = N - NEL
          DO 260 PME = PME1, PME2
             I = IW (PME)
             NVI = -NV (I)
             IF (NVI .LE. 0) GO TO 260
             NV (I) = NVI
             IF (DEGREE(I).GT.N) GO TO 258
             DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
             DEGREE (I) = DEG
             IDENSE = .FALSE.
             IF (THRESM.GE.0) THEN
                IF ((DEG+NVI .GE. THRESM).OR.
     &               (DEG+NVI .GE. NLEFT)) THEN
                   IF (THRESM.EQ.N) THEN
                      IF ((ELEN(I).LE.2) .AND.((DEG+NVI).EQ.NLEFT)
     &                     .AND. NBD.EQ.NFULL ) THEN
                         DEGREE(I) = N+1
                         IDENSE = .TRUE.
                      ENDIF
                   ELSE
                      IDENSE = .TRUE.
                      IF ((ELEN(I).LE.2).AND. ((DEG+NVI).EQ.NLEFT)
     &                     .AND. NBD.EQ.NFULL ) THEN
                         DEGREE(I) = N+1
                      ELSE
                         DEGREE(I) = N+1+DEGREE(I)
                      ENDIF
                   ENDIF
                ENDIF
                IF (IDENSE) THEN
                   P1 = PE(I)
                   P2 = P1 + ELEN(I) - 1
                   DO 255 PJ=P1,P2
                      E= IW(PJ)
                      DENXT (E) = DENXT(E) + NVI
 255               CONTINUE
                   NBD = NBD+NVI
                   DEG = N
                   IF (DEGREE(I).EQ.N+1) THEN
                      NFULL = NFULL +NVI
                      IF (LASTD.EQ.0) THEN
                         LASTD     = I
                         HEAD(DEG) = I
                         DENXT(I)   = 0
                         LAST(I)   = 0
                      ELSE
                            DENXT(LASTD) = I
                            LAST(I)     = LASTD
                            LASTD       = I
                            DENXT(I)     = 0
                      ENDIF
                   ELSE
                      INEXT = HEAD(DEG)
                      IF (INEXT .NE. 0) LAST (INEXT) = I
                      DENXT (I) = INEXT
                      HEAD (DEG) = I
                      LAST(I)    = 0
                      IF (LASTD.EQ.0) LASTD=I
                   ENDIF
                ENDIF
             ENDIF
             IF (.NOT.IDENSE) THEN
                INEXT = HEAD (DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
                DENXT (I) = INEXT
                LAST (I) = 0
                HEAD (DEG) = I
             ENDIF
             MINDEG = MIN (MINDEG, DEG)
 258         CONTINUE
             IW (P) = I
             P = P + 1
 260      CONTINUE
          OPS = OPS + DEGME*NVPIV + DEGME * NVPIV*NVPIV +
     *         DEGME*DEGME*NVPIV + NVPIV*NVPIV*NVPIV/3 +
     *         NVPIV*NVPIV/2 + NVPIV/6 + NVPIV
          NRLADU = NRLADU + (NVPIV*(NVPIV+1))/2 + (DEGME*NVPIV)
          NV (ME) = NVPIV + DEGME
          LEN (ME) = P - PME1
          IF (LEN (ME) .EQ. 0) THEN
             PE (ME) = 0
             W (ME) = 0
          ENDIF
          IF (NEWMEM .NE. 0) THEN
             PFREE = P
             MEM = MEM - NEWMEM + LEN (ME)
          ENDIF
          GO TO 30
       ENDIF
       GO TO 265
 263   NELME    = -(NEL+1)
       DO 264 X=1,N
          IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
          ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
          ENDIF
 264   CONTINUE
       ELEN(ME) = NELME
       NV(ME)   = NBD
       NRLADU = NRLADU + (NBD*(NBD+1))/2
       OPS = OPS + NBD*NBD*NBD/3 + NBD*NBD/2 + NBD/6 + NBD
       PE(ME)   = 0
 265   CONTINUE
       DO 290 I = 1, N
          IF (ELEN (I) .EQ. 0) THEN
             J = -PE (I)
             DO 270 JDUMMY = 1,N
                IF (ELEN (J) .LT. 0) GO TO 275
                J = -PE (J)
 270         CONTINUE
 275         E = J
             K = -ELEN (E)
             J = I
             DO 280 IDUMMY = 1,N
                IF (ELEN (J) .LT. 0) GO TO 285
                JNEXT = -PE (J)
                PE (J) = -E
                IF (ELEN (J) .EQ. 0) THEN
                   ELEN (J) = K
                   K = K + 1
                ENDIF
                J = JNEXT
 280         CONTINUE
 285         ELEN (E) = -K
          ENDIF
 290   CONTINUE
       DO 300 I = 1, N
          K = ABS (ELEN (I))
          LAST (K) = I
          ELEN (I) = K
 300   CONTINUE
       RJNFO(1) = OPS
       RJNFO(2) = NRLADU
       JNFO(1) = NCMPA
       JNFO(2) = RSTRT
       PFREE = MAXMEM
       RETURN
       END
* *******************************************************************
* COPYRIGHT (c) 1988 Hyprotech UK and
* Council for the Central Laboratory of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION EST
      INTEGER KASE,N
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTRINSIC ABS,SIGN,NINT,DBLE
      IF (N.LE.0) THEN
        KASE = -1
        RETURN
      END IF
      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN
      END IF
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
      GO TO (100,200,300,400,500) JUMP
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
        GO TO 510
      END IF
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
  300 CONTINUE
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
      GO TO 410
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220
      END IF
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
  510 KASE = 0
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1987 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 10 Feb 1993
C       Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
      INTEGER N
      LOGICAL YESA
      DOUBLE PRECISION A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
      OLDTAU = JCOLST(N+1) - 1
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   10   CONTINUE
   20 CONTINUE
      NEWTAU = 2*OLDTAU - NDIAG
      IPKP1 = OLDTAU + 1
      CKP1 = NEWTAU + 1
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
        LENK = I2 - I1
        JSTART = CKP1
        IPKP1 = I1
        I2 = I2 - 1
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
        JCOLST(J) = JSTART
        CKP1 = CKP1 - IW(J)
        IW(J) = LENK
   40 CONTINUE
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
      INTRINSIC MAX
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
      LP = ICNTL(4)
      MP = ICNTL(5)
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF
      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF
      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF
      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF
      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF
      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0
      PART = 0
      IF (ICNTL6.NE.0) PART = 1
      IF (ICNTL2.EQ.0) THEN
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)
      ELSE IF (ICNTL2.EQ.1) THEN
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)
      ELSE IF (ICNTL2.EQ.2) THEN
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
          IF (INFO(1).EQ.-9) GO TO 40
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)
      END IF
      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70
   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70
   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN
 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE
      KNE = 0
      IOUT = 0
      JOUT = 0
      IF (LCHECK) THEN
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
        IF (KNE.EQ.0) GO TO 130
      ELSE
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF
      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE
      IF (LA.EQ.1) THEN
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF
  130 CONTINUE
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      INTEGER LA,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE
      IF (LA.GT.1) THEN
        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE
        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE
        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
      INTEGER LA,NC,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC ABS
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************
      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
      INTEGER I,J,K,KSTART,KSTOP,NZJ
      IDUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER
      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C *******************************************************************
C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date July 1999
CCCCC PACKAGE MC64A/AD
CCCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and
CCCCC         Jacko Koster (jacko.koster@uninett.no)

C 12th July 2004 Version 1.0.0. Version numbering added.

C 30/07/04  Version 1.1.0. Permutation array flagged negative to
C           indicate dependent columns in singular case.  Calls to
C           MC64F/FD changed to avoid unsafe reference to array L.
C 21st February 2005 Version 1.2.0. FD05 dependence changed to FD15.
C  7th March 2005 Version 1.3.0. Scan of dense columns avoided in
C           the cheap assignment phase in MC64W/WD.
C 15th May 2007 Version 1.4.0. Minor change made in MC64W/WD to avoid
C           crash when the input array contains NaNs.
C 28th June 2007  Version 1.5.0. Permutation array flagged negative to
C           indicate dependent columns in singular case for all values
C           of the JOB parameter (previously was only for JOB 4 and 5).
C           Also some very minor changes concerning tests in DO loops
C           being moved from beginning to end of DO loop.


      SUBROUTINE MC64ID(ICNTL)
      IMPLICIT NONE
      INTEGER ICNTL(10)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      DOUBLE PRECISION A(NE),DW(LDW)
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
      EXTERNAL FD15AD,MC21AD,MC64BD,MC64RD,MC64SD,MC64WD
      DOUBLE PRECISION FD15AD
      INTRINSIC ABS,LOG
      RINF = FD15AD('H')
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE
      IF (JOB.EQ.1) THEN
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
        CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
        CALL MC64RD(N,NE,IP,IW,DW)
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
      ENDIF
   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF
   99 RETURN
 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END
C**********************************************************************
      SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),D(N)
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
      INTRINSIC ABS,MIN
      EXTERNAL FD15AD,MC64DD,MC64ED,MC64FD
      DOUBLE PRECISION FD15AD
      RINF = FD15AD('H')
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
        CSP = MINONE
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
          ENDIF
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.MINONE) GO TO 190
        BV = MIN(BV,CSP)
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = - I
  320 CONTINUE
 1000 RETURN
      END
C**********************************************************************
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI
      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
      ENDIF
   40 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
      INTEGER TODO(TDLEN)
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
        TD = TD - 2
        GO TO 425
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2
  425   CONTINUE
        IF (TD.EQ.0) GO TO 400
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
        TD = TD - 2
        GO TO 425
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      DOUBLE PRECISION A(NE)
      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
      EXTERNAL FD15AD,MC64QD,MC64UD
      DOUBLE PRECISION FD15AD
      RINF = FD15AD('H')
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))
      NUM = NUMX
      IF (NUM.NE.N) THEN
        BMAX = RINF
      ELSE
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF
      BVAL = 0.0
      BMIN = 0.0
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
              LEN(J) = II - IP(J) + 1
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
   81     MOD = 1
        ELSE
          BMAX = BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))
   90 CONTINUE
   99 IF (NUMX.EQ.N) GO TO 1000
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) =  - J
  320 CONTINUE
 1000 RETURN
      END
C**********************************************************************
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA
      NVAL = 0
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1
      IF (ID.EQ.1) THEN
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM
      NFC = 0
      ID0 = (ID-1)*N
      DO 100 JORD = NUM0+1,N
        ID1 = ID0 + JORD
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = LENC(J)
   30     OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
   70   CONTINUE
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
   95   IF (NUM.EQ.NUM1) THEN
          LAST = JORD
          GO TO 101
        ENDIF
  100 CONTINUE
      LAST = N
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D)
      IMPLICIT NONE
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),U(N),D(N)
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
      DOUBLE PRECISION RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
      EXTERNAL FD15AD,MC64DD,MC64ED,MC64FD
      DOUBLE PRECISION FD15AD
      RINF = FD15AD('H')
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
        IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
        IF (K1.GT.K2) GO TO 95
        I0 = IRN(K1)
        VJ = A(K1) - U(I0)
        K0 = K1
        DO 50 K = K1+1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
        CSP = RINF
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64ED(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = VJ + A(K)-U(I)
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                LPOS = L(I)
                IF (LPOS.NE.0)
     *            CALL MC64FD(LPOS,QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.RINF) GO TO 190
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J)
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE
      IF (NUM.EQ.N) GO TO 1100
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1977 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      INTEGER LICN,N,NUMNZ
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
      EXTERNAL MC21BD
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER LICN,N,NUMNZ
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENR(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
  100 CONTINUE
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
      END
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

      DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C       
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C          EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C          HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15AD = EPSILON(1.0D0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15AD = TINY(1.0D0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15AD = HUGE(1.0D0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15AD = DBLE(RADIX(1.0D0))
      ELSE
         FD15AD = 0.0D0
      ENDIF
      RETURN
      END
      SUBROUTINE METIS_NODEND(N,IPTR,IRN,METFTN,METOPT,INVPRM,PERM)
C Dummy routine that is called if MeTiS is not linked.
      INTEGER N
      INTEGER IPTR(N+1),IRN(*),METFTN,METOPT(8),INVPRM(N),PERM(N)
      PERM(1) = -1
      RETURN
      END
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      SUBROUTINE DTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
************************************************************************
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTPSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix, supplied in packed form.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  AP     - DOUBLE PRECISION array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U' or 'u', the array AP must
*           contain the upper triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
*           respectively, and so on.
*           Before entry with UPLO = 'L' or 'l', the array AP must
*           contain the lower triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
*           respectively, and so on.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, K, KK, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTPSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of AP are
*     accessed sequentially with one pass through AP.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KK = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/AP( KK )
                     TEMP = X( J )
                     K    = KK     - 1
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*AP( K )
                        K      = K      - 1
   10                CONTINUE
                  END IF
                  KK = KK - J
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/AP( KK )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, K = KK - 1, KK - J + 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*AP( K )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
                  KK = KK - J
   40          CONTINUE
            END IF
         ELSE
            KK = 1
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/AP( KK )
                     TEMP = X( J )
                     K    = KK     + 1
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*AP( K )
                        K      = K      + 1
   50                CONTINUE
                  END IF
                  KK = KK + ( N - J + 1 )
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/AP( KK )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, K = KK + 1, KK + N - J
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*AP( K )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
                  KK = KK + ( N - J + 1 )
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A' )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KK = 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  K    = KK
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - AP( K )*X( I )
                     K    = K    + 1
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/AP( KK + J - 1 )
                  X( J ) = TEMP
                  KK     = KK   + J
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, K = KK, KK + J - 2
                     TEMP = TEMP - AP( K )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/AP( KK + J - 1 )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  KK      = KK   + J
  120          CONTINUE
            END IF
         ELSE
            KK = ( N*( N + 1 ) )/2
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  K = KK
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - AP( K )*X( I )
                     K    = K    - 1
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/AP( KK - N + J )
                  X( J ) = TEMP
                  KK     = KK   - ( N - J + 1 )
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, K = KK, KK - ( N - ( J + 1 ) ), -1
                     TEMP = TEMP - AP( K )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/AP( KK - N + J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  KK      = KK   - (N - J + 1 )
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTPSV .
*
      END
C       Toolpack tool decs employed.
C       Arg dimension set to *.
C
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
 
c      FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
c      JACK DONGARRA, LINPACK, 3/11/78.
 
c      .. Scalar Arguments ..
      INTEGER INCX,N
c      ..
c      .. Array Arguments ..
      DOUBLE PRECISION DX(*)
c      ..
c      .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
c      ..
c      .. Intrinsic Functions ..
      INTRINSIC DABS
c      ..
c      .. Executable Statements ..
 
      IDAMAX = 0
      IF (N.LT.1) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
 
c         CODE FOR INCREMENT NOT EQUAL TO 1
 
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
        IF (DABS(DX(IX)).LE.DMAX) GO TO 5
        IDAMAX = I
        DMAX = DABS(DX(IX))
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
 
c         CODE FOR INCREMENT EQUAL TO 1
      
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
        IF (DABS(DX(I)).LE.DMAX) GO TO 30
        IDAMAX = I
        DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
       
      END
C       A BLAS routine modified for use with HSL
C       Toolpack tool decs employed.
C       Contained comment lines which caused Decs to fail.
C
      SUBROUTINE XERBLA(SRNAME,INFO)
C     .. Scalar Arguments ..
      INTEGER INFO
      CHARACTER SRNAME*6
C     ..
C
C  Purpose
C  =======
C
C  XERBLA  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Installers should consider modifying the STOP statement in order to
C  call system-specific exception-handling facilities.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*6.
C           On entry, SRNAME specifies the name of the routine which
C           called XERBLA.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Executable Statements ..
C
      WRITE (*,FMT=99999) SRNAME,INFO
C
      STOP
C
99999 FORMAT (' ** On entry to ',A6,' parameter number ',I2,
     +       ' had an illegal value')
C
C     End of XERBLA.
C
      END
      LOGICAL FUNCTION LSAME ( CA, CB )
*     .. Scalar Arguments ..
      CHARACTER*1            CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME  tests if CA is the same letter as CB regardless of case.
*
*  N.B. This version of the routine is only correct for ASCII code.
*       Installers must modify the routine for other character-codes.
*
*       For EBCDIC systems the constant IOFF must be changed to -64.
*       For CDC systems using 6-12 bit representations, the system-
*       specific code in comments must be activated.
*
*  Parameters
*  ==========
*
*  CA     - CHARACTER*1
*  CB     - CHARACTER*1
*           On entry, CA and CB specify characters to be compared.
*           Unchanged on exit.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  -- Written on 11-October-1988.
*     Richard Hanson, Sandia National Labs.
*     Jeremy Du Croz, Nag Central Office.
*
*     .. Parameters ..
      INTEGER                IOFF
      PARAMETER            ( IOFF=32 )
*     .. Intrinsic Functions ..
      INTRINSIC              ICHAR
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
*
*     Now test for equivalence
*
      IF ( .NOT.LSAME ) THEN
         LSAME = ICHAR(CA) - IOFF .EQ. ICHAR(CB)
      END IF
      IF ( .NOT.LSAME ) THEN
         LSAME = ICHAR(CA) .EQ. ICHAR(CB) - IOFF
      END IF
*
      RETURN
*
*  The following comments contain code for CDC systems using 6-12 bit
*  representations.
*
*     .. Parameters ..
*     INTEGER                ICIRFX
*     PARAMETER            ( ICIRFX=62 )
*     .. Scalar Arguments ..
*     CHARACTER*1            CB
*     .. Array Arguments ..
*     CHARACTER*1            CA(*)
*     .. Local Scalars ..
*     INTEGER                IVAL
*     .. Intrinsic Functions ..
*     INTRINSIC              ICHAR, CHAR
*     .. Executable Statements ..
*
*     See if the first character in string CA equals string CB.
*
*     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
*
*     IF (LSAME) RETURN
*
*     The characters are not identical. Now check them for equivalence.
*     Look for the 'escape' character, circumflex, followed by the
*     letter.
*
*     IVAL = ICHAR(CA(2))
*     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
*        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
*     END IF
*
*     RETURN
*
*     End of LSAME.
*
      END
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
*
************************************************************************
*
*     File of the DOUBLE PRECISION  Level-2 BLAS.
*     ===========================================
*
*     SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
*    $                   BETA, Y, INCY )
*
*     SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
*    $                   BETA, Y, INCY )
*
*     SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
*    $                   BETA, Y, INCY )
*
*     SUBROUTINE DSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
*    $                   BETA, Y, INCY )
*
*     SUBROUTINE DSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
*
*     SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
*     SUBROUTINE DTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
*     SUBROUTINE DTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
*     SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*
*     SUBROUTINE DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*
*     SUBROUTINE DTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
*
*     SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
*     SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
*
*     SUBROUTINE DSPR  ( UPLO, N, ALPHA, X, INCX, AP )
*
*     SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*
*     SUBROUTINE DSPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
*
*     See:
*
*        Dongarra J. J., Du Croz J. J., Hammarling S.  and Hanson R. J..
*        An  extended  set of Fortran  Basic Linear Algebra Subprograms.
*
*        Technical  Memoranda  Nos. 41 (revision 3) and 81,  Mathematics
*        and  Computer Science  Division,  Argonne  National Laboratory,
*        9700 South Cass Avenue, Argonne, Illinois 60439, US.
*
*        Or
*
*        NAG  Technical Reports TR3/87 and TR4/87,  Numerical Algorithms
*        Group  Ltd.,  NAG  Central  Office,  256  Banbury  Road, Oxford
*        OX2 7DE, UK,  and  Numerical Algorithms Group Inc.,  1101  31st
*        Street,  Suite 100,  Downers Grove,  Illinois 60515-1263,  USA.
*
************************************************************************
*
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END

C     lssinfo:
C
C     0: Success.
C     1: Matrix not positive definite.
C     2: Rank deficient matrix.
C     6: Insufficient space to store the linear system.
C     7: Insufficient double precision working space.
C     8: Insufficient integer working space.

C     ******************************************************************
C     ******************************************************************

      logical function lss(lsssub)

C     SCALAR ARGUMENTS
      character * 4 lsssub

      lsssub = 'MA57'
      lss    = .true.

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssini(sclsys,acneig,usefac)

      implicit none

C     SCALAR ARGUMENTS
      logical acneig,sclsys,usefac

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/
C     COMMON SCALARS
      double precision bignum,macheps,macheps12,macheps13,macheps23

C     COMMON BLOCKS
      common /machcon/ bignum,macheps,macheps12,macheps13,macheps23
      save   /machcon/

      lacneig = acneig
      lusefac = usefac
      lsclsys = sclsys

      call ma57id(cntl,icntl)

C     Suppress monitoring, warning and error messages
      icntl(5) = 0

C     Chooses AMD pivot ordering using MC47
C     icntl(6) = 0

      if ( lsclsys ) then
          icntl(15) = 1
      else
          icntl(15) = 0
      end if
      
      if ( .not. lacneig ) then
          icntl(7)  = 2
          icntl(8)  = 1
          
          cntl(2)   = macheps23
      end if
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssana(nsys,hnnz,hlin,hcol,hval,hdiag,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,lssinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(10,9000) nsysmax,nsys
          end if
          
          lssinfo = 6
          return
      end if

      call ma57ad(nsys,hnnz,hlin,hcol,wintmax,keep,iwork,icntl,info,
     +rinfo)

      if ( info(1) .eq. 0 ) then
        ! SUCCESS

          lssinfo = 0
          return
      end if

C     UNHANDLED ERROR

      if ( iprintctl(3) ) then
          write(* ,9030) info(1)
          write(10,9030) info(1)
      end if
      
      stop

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSANA-MA57 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSANA-MA57 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssfac(nsys,hnnz,hlin,hcol,hval,hdiag,d,pind,pval,
     +nneigv,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,pind,nneigv,lssinfo
      double precision pval

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz),d(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer i,idiag
      double precision d1

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(10,9000) nsysmax,nsys
          end if
          
          lssinfo = 6
          return
      end if

      do i = 1,nsys
          idiag       = hdiag(i)
          w(i)        = hval(idiag)
          hval(idiag) = hval(idiag) + d(i)
      end do

      call ma57bd(nsys,hnnz,hval,fact,fnnzmax,ifact,wintmax,wintmax,
     +keep,iwork,icntl,cntl,info,rinfo)

      if ( lsclsys ) then
          do i = 1,nsys
              s(i) = fact(fnnzmax-nsys-1+i)
          end do    
      end if

      do i = 1,nsys
          idiag       = hdiag(i)
          hval(idiag) = w(i)
      end do
      
      if ( info(1) .eq. 0 ) then

          d1 = hval(hdiag(1)) + d(1)

          if ( d1 .gt. 0.0d0 ) then
            ! SUCCESS
          
              lssinfo = 0
          else
            ! MATRIX IS NEGATIVE DEFINITE

              pind    = 1
              pval    = abs( d1 )
              
              lssinfo = 1
          end if

      else if ( info(1) .eq. 5 .or. info(1) .eq. -6 ) then
        ! MATRIX NOT POSITIVE DEFINITE

          pind    = info(2)
          pval    = abs( rinfo(20) )
          
          lssinfo = 1

      else if ( info(1) .eq. 4 .or. info(1) .eq. -5 ) then
        ! RANK DEFICIENT MATRIX

          pind    = info(2)
          pval    = abs( rinfo(20) )

          lssinfo = 2

      else if ( info(1) .eq. -3 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          if ( iprintctl(3) ) then
              write(* ,9010) fnnzmax,info(17)
              write(10,9010) fnnzmax,info(17)
          end if
          
          lssinfo = 7
          
      else if ( info(1) .eq. -4 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          if ( iprintctl(3) ) then
              write(* ,9020) wintmax,info(18)
              write(10,9020) wintmax,info(18)
          end if

          lssinfo = 8

      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) info(1)
              write(10,9030) info(1)
          end if

          stop

      end if

C     NUMBER OF NEGATIVE EIGENVALUES
      nneigv = info(24)

      if ( lusefac .and. lssinfo .eq. 0 ) then

C        Define matrix D^{-1} (stored in sdiag)
          
          do i = 1,nsys
              w(i) = 0.0d0
          end do
          
          do i = 1,nsys
              
              w(i) = 1.0d0

              call ma57cd(3,nsys,fact,fnnzmax,ifact,wintmax,1,w,nsys,
     +        work,nsysmax,iwork,icntl,info)

              sdiag(i) = w(i)
              w(i)     = 0.0d0
          end do
          
      end if

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSFAC-MA57 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSFAC-MA57 WARNING: Insufficient double precision ',
     +            'working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9020 format(/,1X,'LSSFAC-MA57 WARNING: Insufficient integer working ',
     +            'space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSFAC-MA57 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssol(nsys,sol)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      double precision sol(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

      call ma57cd(1,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +nsysmax,iwork,icntl,info)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssoltr(job,nsys,sol)

      implicit none

C     SCALAR ARGUMENTS
      character * 1 job
      integer nsys

C     ARRAY ARGUMENTS
      double precision sol(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer i

      if ( job .eq. 'T' .or. job .eq. 't' ) then

          call ma57cd(2,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +    nsysmax,iwork,icntl,info)

          do i = 1,nsys
              sol(i) = sol(i) * sqrt( sdiag(i) )
          end do
          
          do i = 1,nsys
              work(i) = sol(i)
          end do

          do i = 1,nsys
              sol(keep(i)) = work(i)
          end do

      else

          do i = 1,nsys
              work(i) = sol(i)
          end do

          do i = 1,nsys
              sol(keep(i)) = work(i)
          end do

          do i = 1,nsys
              sol(i) = sol(i) * sqrt( sdiag(i) )
          end do
          
          call ma57cd(4,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +    nsysmax,iwork,icntl,info)

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsspermvec(nsys,v)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      integer v(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer i

      do i = 1,nsys
          iwork(keep(i)) = v(i)
      end do
      
      do i = 1,nsys
          v(i) = iwork(i)
      end do
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssunpermvec(nsys,v)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     ARRAY ARGUMENTS
      integer v(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer i

      do i = 1,nsys
          iwork(i) = v(keep(i))
      end do
      
      do i = 1,nsys
          v(i) = iwork(i)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsspermind(hnnz,hlin,hcol)

      implicit none

C     SCALAR ARGUMENTS
      integer hnnz

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer col,i,lin

      do i = 1,hnnz
          col = hcol(i)
          lin = hlin(i)
          
          hcol(i) = keep(col)
          hlin(i) = keep(lin)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssunpermind(nsys,hnnz,hlin,hcol)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer col,i,lin

      do i = 1,nsys
          invp(i)=0
      end do
      
      do i = 1,nsys
          invp(keep(i)) = i
      end do
      
      do i = 1,hnnz
          col = hcol(i)
          lin = hlin(i)
          
          hcol(i) = invp(col)
          hlin(i) = invp(lin)
      end do

      end

C     ******************************************************************
C     ******************************************************************

      double precision function lssgetd(j)

      implicit none

C     SCALAR ARGUMENTS
      integer j

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

      lssgetd = sqrt( 1.0d0 / sdiag(invp(j)) )

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lsssetrow(nsys)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer apos,i,nblk,iblk,ipiv,ix,iwpos,j1,j2,k,ncols,nrows,totlen,
     +        nsys2

C     Count the number of elements of fact and ifact used to form each
C     row of L^T

      totlen  = 0
      nsys2   = 2*nsys
      
      do i = 1,nsys
          iwork(nsys+i) = 0
      end do

      do i = 1,nsys
          iwork(nsys2+i) = 0
      end do

      do i = 1,nsys
          iwork(i) = 0
      end do

      apos  = 1
      iwpos = 4
      nblk  = ifact(3)
      
      do iblk = 1,nblk
          ncols = ifact(iwpos)
          nrows = ifact(iwpos+1)
          iwpos = iwpos + 2
          j1    = iwpos
          j2    = iwpos + nrows - 1
          
          do ipiv = 1,nrows
              apos = apos + 1
              ix   = abs( ifact(j1) )
              k    = apos
              
              if ( j1+1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k    = max( k, k+j2-j1 )
              apos = k
              j1   = j1 + 1
          end do
          
          j2 = iwpos + ncols - 1
          
          do ipiv = 1,nrows-1,2
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k  = apos+ncols-nrows
              ix = abs( ifact(iwpos+ipiv) )
              
              if ( j1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end do
          
          if ( mod(nrows,2) .eq. 1 ) then
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  iwork(nsys2+ix) = iwork(nsys2+ix) + 3
                  totlen          = totlen + 3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end if
          
          iwpos = iwpos + ncols
          
      end do
      
      iwork(nsys+1) = 1
      
      do i = 2,nsys
          iwork(nsys+i) = iwork(nsys+i-1) + iwork(nsys2+i-1)
      end do
      
      do i = 1,nsys
          iwork(i) = iwork(nsys+i)
      end do

C     Store the indices of elements of fact and ifact used to form each
C     row of L^T

      apos  = 1
      iwpos = 4
      nblk  = ifact(3)
      
      do iblk = 1,nblk
          ncols = ifact(iwpos)
          nrows = ifact(iwpos+1)
          iwpos = iwpos + 2
          j1    = iwpos
          j2    = iwpos + nrows - 1
          
          do ipiv = 1,nrows
              apos = apos + 1
              ix   = abs( ifact(j1) )
              k    = apos
              
              if ( j1+1 .le. j2 ) then
                  posfac(iwork(ix))   = j1+1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k   = max( k, k+j2-j1 )
              apos = k
              j1   = j1 + 1
          end do
          
          j2 = iwpos + ncols - 1
          
          do ipiv = 1,nrows-1,2
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  posfac(iwork(ix))   = j1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = -k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k  = apos+ncols-nrows
              ix = abs( ifact(iwpos+ipiv) )
              
              if ( j1 .le. j2 ) then
                  posfac(iwork(ix))   = j1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = -k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end do
          
          if ( mod(nrows,2) .eq. 1 ) then
              k  = apos
              ix = abs( ifact(iwpos+ipiv-1) )
              
              if ( j1 .le. j2 ) then
                  posfac(iwork(ix))   = j1
                  posfac(iwork(ix)+1) = j2
                  posfac(iwork(ix)+2) = -k
                  iwork(ix)           = iwork(ix)+3
              end if
              
              k    = max( k, k+j2-j1+1 )
              apos = k
          end if
          
          iwpos = iwpos + ncols
          
      end do
      
C     Array of null elements

      do i = 1,nsys
          iwork(i) = 0
      end do

C     Compute inverse permutation

      do i = 1,nsys
          invp(keep(i)) = i
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssgetrow(nsys,idx,rownnz,rowind,rowval)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,idx,rownnz

C     ARRAY ARGUMENTS
      integer rowind(nsys)
      double precision rowval(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer i,j,j1,j2,k,ind,ista,iend,ix
      double precision sclf

      ind  = invp(idx)
      ista = iwork(nsys+ind)
      iend = ista + iwork(2*nsys+ind) - 1

      if ( lsclsys ) then
          sclf = s(ind)
      else
          sclf = 1.0d0
      end if

      rownnz         = 1
      rowind(rownnz) = ind
      rowval(rownnz) = sclf
      iwork(ind)     = rownnz

      do i = ista,iend,3
          j1 = posfac(i)
          j2 = posfac(i+1)
          k  = posfac(i+2)
          
          if ( k .gt. 0 ) then
              
              do j = j1,j2
                  ix = abs( ifact(j) )
                  
                  if ( iwork(ix) .eq. 0 ) then
                      if ( fact(k) .ne. 0.0d0 ) then
                          rownnz         = rownnz + 1
                          rowind(rownnz) = ix
                          rowval(rownnz) = sclf * fact(k)
                          iwork(ix)      = rownnz
                      end if
                  else
                      rowval(iwork(ix)) = 
     +                rowval(iwork(ix)) + sclf * fact(k)
                  end if
                  
                  k = k + 1
              end do
          else
              k = abs( k )
              
              do j = j1,j2
                  ix = abs( ifact(j) )
                  
                  if ( iwork(ix) .eq. 0 ) then
                      if ( fact(k) .ne. 0.0d0 ) then
                          rownnz         = rownnz + 1
                          rowind(rownnz) = ix
                          rowval(rownnz) = - sclf * fact(k)
                          iwork(ix)      = rownnz
                      end if
                  else
                      rowval(iwork(ix)) = 
     +                rowval(iwork(ix)) - sclf * fact(k)
                  end if
                  
                  k = k + 1
              end do
          end if
      end do

      if ( lsclsys ) then
          do i = 1,rownnz
              rowval(i) = rowval(i) / s(rowind(i))
          end do
      end if

C     Set used positions of iwork back to 0

      do i = 1,rownnz
          iwork(rowind(i)) = 0
      end do

      do i = 1,rownnz
          rowval(i) = rowval(i) / sqrt( sdiag(ind) )
      end do

      do i = 1,rownnz
          rowind(i) = keep(rowind(i))
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine lssafsol(nsys,hnnz,hlin,hcol,hval,hdiag,d,sol,lssinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,lssinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz),d(nsys),sol(nsys)

C     PARAMETERS

      integer mmax,nmax,nsmax,jcnnzmax,hnnzmax,fnnzmax,wintmax,nsysmax

      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( nsmax     =     1000 )
      parameter ( jcnnzmax  = 10000000 )
      parameter ( hnnzmax   = 10000000 )
      parameter ( fnnzmax   = 10000000 )
      parameter ( wintmax   = 10000000 )
      parameter ( nsysmax   =   100000 )
C     COMMON SCALARS
      logical iprintctl(6)
      integer iprintinn,iprintout,mprint,nprint

C     COMMON BLOCKS
      common /outdat/ iprintctl,iprintinn,iprintout,nprint,mprint
      save   /outdat/

C     COMMON SCALARS
      logical lacneig,lsclsys,lusefac

C     COMMON ARRAYS
      integer icntl(20),iwork(5*nsysmax),posfac(wintmax),invp(nsysmax),
     +        ifact(wintmax),keep(wintmax),info(40)
      double precision cntl(5),work(nsysmax),w(nsysmax),sdiag(nsysmax),
     +        fact(fnnzmax),s(nsysmax),rinfo(20)

C     COMMON BLOCKS
      common /lssdat/ cntl,work,w,sdiag,fact,s,rinfo,icntl,iwork,posfac,
     +                invp,ifact,keep,info,lacneig,lsclsys,lusefac
      save   /lssdat/

C     LOCAL SCALARS
      integer i,tmp

C     LOCAL ARRAYS
      integer keeptmp(wintmax),infotmp(40)
      double precision rinfotmp(20)

      tmp       = icntl(15)
      icntl(15) = 0

      if ( nsys .gt. nsysmax ) then
        ! INSUFFICIENT SPACE TO STORE THE LINEAR SYSTEM

          lssinfo = 6

          if ( iprintctl(3) ) then
              write(* ,9000) nsysmax,nsys
              write(10,9000) nsysmax,nsys
          end if
          
          go to 500
          
      end if

      call ma57ad(nsys,hnnz,hlin,hcol,wintmax,keeptmp,iwork,icntl,
     +infotmp,rinfotmp)

      if ( ( infotmp(1) .eq. 0 ) .or. ( infotmp(1) .eq. 1 ) ) then
        ! SUCCESS

          lssinfo = 0

      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) infotmp(1)
              write(10,9030) infotmp(1)
          end if
          
          stop
          
      end if
      
      do i = 1,nsys
          hval(hdiag(i)) = hval(hdiag(i)) + d(i)
      end do

      if ( lsclsys ) then
          do i = 1,hnnz
              hval(i) = hval(i) * s(hlin(i)) * s(hcol(i))
          end do
      end if

      call ma57bd(nsys,hnnz,hval,fact,fnnzmax,ifact,wintmax,wintmax,
     +keeptmp,iwork,icntl,cntl,infotmp,rinfotmp)

      if ( lsclsys ) then
          do i = 1,hnnz
              hval(i) = hval(i) * s(hlin(i)) * s(hcol(i))
          end do
      end if

      do i = 1,nsys
          hval(hdiag(i)) = hval(hdiag(i)) - d(i)
      end do

      if ( ( infotmp(1) .eq. 0 ) .or. ( infotmp(1) .eq. 1 ) ) then
        ! SUCCESS

          lssinfo = 0

      else if ( info(1) .eq. 5 .or. info(1) .eq. -6 ) then
        ! MATRIX NOT POSITIVE DEFINITE

          lssinfo = 1

      else if ( infotmp(1) .eq. 4 .or. infotmp(1) .eq. -5 ) then
        ! RANK DEFICIENT MATRIX

          lssinfo = 2
         
      else if ( infotmp(1) .eq. -3 ) then
        ! INSUFFICIENT DOUBLE PRECISION WORKING SPACE

          lssinfo = 7

          if ( iprintctl(3) ) then
              write(* ,9010) fnnzmax,infotmp(17)
              write(10,9010) fnnzmax,infotmp(17)
          end if

          go to 500
          
      else if ( infotmp(1) .eq. -4 ) then
        ! INSUFFICIENT INTEGER WORKING SPACE

          lssinfo = 8

          if ( iprintctl(3) ) then
              write(* ,9020) wintmax,infotmp(18)
              write(10,9020) wintmax,infotmp(18)
          end if
          
          go to 500
          
      else
        ! UNHANDLED ERROR

          if ( iprintctl(3) ) then
              write(* ,9030) infotmp(1)
              write(10,9030) infotmp(1)
          end if
          
          stop
          
      end if

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if
      
      call ma57cd(1,nsys,fact,fnnzmax,ifact,wintmax,1,sol,nsys,work,
     +nsysmax,iwork,icntl,infotmp)

      if ( lsclsys ) then
          do i = 1,nsys
              sol(i) = sol(i) * s(i)
          end do
      end if

 500  continue

      icntl(15) = tmp

C     NON-EXECUTABLE STATEMENTS

 9000 format(/,1X,'LSSAFSOL-MA57 WARNING: Insufficient space to store ',
     +            'linear system. Increase',
     +       /,1X,'parameter nsysmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9010 format(/,1X,'LSSAFSOL-MA57 WARNING: Insufficient double ',
     +            'precision working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9020 format(/,1X,'LSSAFSOL-MA57 WARNING: Insufficient integer ',
     +            'working space. Increase',
     +       /,1X,'parameter fnnzmax from ',I16,' to at least ',I16,
     +       /,1X,'if you would like to try a direct linear solver ',
     +            'again.')
 9030 format(/,1X,'LSSAFSOL-MA57 ERROR: Unhandled error ',I16,'.',
     +       /,1X,'See documentation for details.')

      end
C     sclinfo:
C
C     0: Success.
C     1: Invalid row or column index.
C     2: Duplicate entry.
C     6: Insufficient space to store the input matrix.
C     7: Insufficient double precision working space.

C     ******************************************************************
C     ******************************************************************

      logical function scl(sclsub)

C     SCALAR ARGUMENTS
      character * 4 sclsub

      scl = .false.
      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclini()

      implicit none

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclana(nsys,hnnz,hlin,hcol,hval,hdiag,sclinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,sclinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz),hdiag(nsys)
      double precision hval(hnnz)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine sclsys(nsys,hnnz,hlin,hcol,hval,s,sclinfo)

      implicit none

C     SCALAR ARGUMENTS
      integer nsys,hnnz,sclinfo

C     ARRAY ARGUMENTS
      integer hlin(hnnz),hcol(hnnz)
      double precision hval(hnnz),s(nsys)

      end
