
	implicit none
	real*8 y
	real*8 x
	real*8 Rayleigh
C	real*8 a0, a1, b1
	real*8 as, bs
	dimension as(0:10)
	dimension bs(10)
	integer idx
	real*8 x0, y0
	real*8 mu
	real*8 theta
	real*8 pi
	real*8 res
	integer max_loop
	parameter(max_loop=5000)
	integer ni

	pi = acos(-1d0)
C	print *, 'fourier integration'

C   - [x] get y(:) as function of Rayleigh scattering
C     CCX = CX(I) * CX(I)	
C     PR(1,1) = 0.75d0 * (1.d0 + CCX)
C   - [x] get x(.) for given idx
C   - [-] get y(.) for x(.) by DINPOL()
C   - summation	

C	x = -1d0
C	do while(x <= 1d0) 
C		y = Rayleigh(x)
C		print *, x, y
C		x = x + 1d-1
C	end do

	as(0) = 0d0
	do idx = 1, 10
		as(idx) = 0d0
		bs(idx) = 0d0
	end do

	do idx = 1, max_loop
		call Halton_sequence(idx, x0, y0)
		mu = (x0 * 2d0) - 1d0            ! [-1.0,1.0]
		theta = acos(mu) * 180.0 / pi
C		print *, idx, theta

		as(0) = as(0) + Rayleigh(mu) * cos(mu * dble(0))
		do ni=1,10
			as(ni) = as(ni) + Rayleigh(mu) * cos(mu * dble(ni))
			bs(ni) = bs(ni) + Rayleigh(mu) * sin(mu * dble(ni))
		end do
C		b1 = b1 + Rayleigh(mu) * sin(mu * 1)
C		print *, idx, b1 / idx
	end do

	as(0) = as(0) / pi / max_loop
	do idx = 1, 10
		as(idx) = as(idx) / pi / max_loop
		bs(idx) = bs(idx) / pi / max_loop
	end do

	theta = 0d0
	do while(theta <= 180d0)
		mu = cos(theta * pi / 180d0) 
		
		res = 0.5d0 * as(0)
		do ni=1,10
			res = res + as(ni) * cos(mu * dble(ni))
			res = res + bs(ni) * sin(mu * dble(ni))
		end do

C		print *, theta, Rayleigh(mu)
		print *, theta, res
		theta = theta + 1d-1
	end do

	end

!----------------------------------------------------------------------

	function Rayleigh(x)
	implicit none
	real*8 x, Rayleigh

	Rayleigh = 0.75d0 * (1 + x * x)
	end

	subroutine testFunc


	return 
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

