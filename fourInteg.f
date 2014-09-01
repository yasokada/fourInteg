
	implicit none
	real*8 y
	real*8 x
	real*8 Rayleigh
	real*8 a0, a1, b1
	integer idx
	real*8 x0, y0
	real*8 mu
	real*8 theta
	real*8 pi
	integer max_loop
	parameter(max_loop=5000)

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

	a0 = 0d0
	a1 = 0d0
	b1 = 0d0

	do idx = 1, max_loop
		call Halton_sequence(idx, x0, y0)
		mu = (x0 * 2d0) - 1d0            ! -1.0 to 1.0
		theta = acos(mu) * 180.0 / pi
C		print *, idx, theta
		a0 = a0 + Rayleigh(mu) * cos(theta * 0)
		a1 = a1 + Rayleigh(mu) * cos(theta * 1)
		b1 = b1 + Rayleigh(mu) * sin(theta * 1)
C		print *, idx, b1 / idx
	end do

	a0 = a0 / pi / max_loop
	a1 = a1 / pi / max_loop
	b1 = b1 / pi / max_loop

C	print *, 'a0=', a0
C	print *, 'a1=', a1
C	print *, 'b1=', b1

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

