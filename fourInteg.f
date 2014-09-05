
C 	DONE: Reproduce of Example 2 (a0) of
C     http://www.math24.net/definition-of-fourier-series.html
C   DONE: Reproduce of Rayeligh phase function
C
	implicit none

C	functions
	real*8 Rayleigh

C	variables
	integer MAX_FOUR_TERM
	parameter(MAX_FOUR_TERM=500)  ! 500: arbitrary
	real*8 as, bs
	dimension as(0:MAX_FOUR_TERM)
	dimension bs(MAX_FOUR_TERM)
	integer idx
	real*8 x0, y0
	real*8 mu, theta, pi
	real*8 radian
	real*8 res
	integer qmcTerms
C	parameter(qmcTerms=1000)
	integer ni
	real*8 xrange
	integer argc
	character*80 argv

	integer fourTerms  ! terms for Fourier expansion 
					   ! 2: most fit for Rayleigh

	argc = iargc()
	if (argc .le. 1) then
		print *, 'Command invalid: Execute as follows'
		print *, '   [cmd] [Num Fourier terms (1..499)] [Num QMC Integ]'
		print *, 'For 2 Fourier terms, and 100 QMC terms'
		print *, '   ./a.out 2 1000 > res.2_1000'
		stop
	end if

	call getarg(1, argv)
	read(argv,*) fourTerms
	call getarg(2, argv)
	read(argv,*) qmcTerms

C	fourTerms = 2  ! most fit for Rayleigh

	pi = acos(-1d0)
	xrange = 2d0 * pi

	as = 0d0
	bs = 0d0

	do idx = 1, qmcTerms
		call Halton_sequence(idx, x0, y0)
		theta = (x0 * 360d0)   ! [0., 360.]
		radian = theta * pi / 180d0

		as(0) = as(0) + Rayleigh(radian) * cos(radian * dble(0))
		do ni=1,fourTerms
			as(ni) = as(ni) + Rayleigh(radian) * cos(radian * dble(ni))
			bs(ni) = bs(ni) + Rayleigh(radian) * sin(radian * dble(ni))
		end do
	end do

	as = as * xrange
	bs = bs * xrange
	as = as / pi / dble(qmcTerms)
	bs = bs / pi / dble(qmcTerms)		

	theta = 0d0  ! [0., 180.]
	do while(theta <= 180d0)
		mu = cos(theta * pi / 180d0) 
		
		res = 0.5d0 * as(0)
		do ni=1,fourTerms
			res = res + as(ni) * cos(mu * dble(ni))
			res = res + bs(ni) * sin(mu * dble(ni))
		end do

		print *, theta, res, Rayleigh(mu)

		theta = theta + 1d-1
	end do

	end

!----------------------------------------------------------------------

	function Rayleigh(radian)
	implicit none
	real*8 radian
	real*8 x, Rayleigh

	x = cos(radian)
	Rayleigh = 0.75d0 * (1 + x * x)

	end

!----------------------------------------------------------------------
! 2007/09/17 16:00
!
!  base 2 for x 
!  base 3 for y
!  --------------------------------------------
!  When we pair them up, we get a sequence of points in a
!  unit square: (1/2, 1/3),(1/4, 2/3),(3/4, 1/9),(1/8, 2/9),
!  (3/8, 4/9),(5/8, 5/9),(7/8, 7/9),(1/16, 8/9),(3/16, 1/27)
! -----------------------------------------------------

      subroutine Halton_sequence(i0,x0,y0)

      integer xbase,ybase
      real*8 invxbase,invybase
      integer i0,j0,k0
      integer inp
      
      real*8 x0, y0, facx,facy


      xbase=2
      ybase=3


      inp=i0
      
      invxbase= 1.0 / dble(xbase)
      facx = 1.0 / dble(xbase)
      
      invybase = 1.0 / dble(ybase)
      facy = 1.0 / dble(ybase)
      
      x0 = 0.0
      y0 = 0.0
      
      do while (inp >0) 
         x0 = x0 + (mod(inp,xbase)) * invxbase
         inp = inp / xbase
         invxbase = invxbase * facx
      end do
      
      inp = i0
      do while (inp >0)
         y0 = y0 + (mod(inp,ybase)) * invybase
         inp = inp / ybase
         invybase = invybase * facy
      end do
      

      return
      end

!----------------------------------------------------------------------

