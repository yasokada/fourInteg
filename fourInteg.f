
	implicit none
	real*8 y
	real*8 x
	real*8 Rayleigh

	print *, 'fourier integration'

C   - get y(:) as function of Rayleigh scattering
C     CCX = CX(I) * CX(I)	
C     PR(1,1) = 0.75d0 * (1.d0 + CCX)
C   - get x(.) for given idx
C   - get y(.) for x(.) by DINPOL()
C   - summation	


	x = 0.d0
	do while(x <= 1.d0) 
		y = Rayleigh(x)
		print *, x, y
		x = x + 1d-1
	end do

C	stop "main"
	end


	function Rayleigh(x)
	implicit none
	real*8 x, Rayleigh

	Rayleigh = 0.75d0 * (1 + x * x)
	end

	subroutine testFunc


	return 
	end