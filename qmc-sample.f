! Usage : ./a.out

      implicit none

      character*20 ctmp
      integer i0,j0,k0     
      real*8 x0,y0
   
      real*8 sum, res

      sum = 0d0

      do k0 = 1, 1000
         call Halton_sequence(k0, x0, y0)
C         ! case 1: integration 0..1 for y=x. result = 1/2 [checked]
C         sum = sum + x0

C         ! case 2: integration 0..1 for y=x*x. result = 1/3 [checked]
C         sum = sum + (x0 * x0)

         ! case 3: integartion 0..1 for y=(2x+1). result = 1+1 = 2 [checked]     
         sum = sum + (2d0 * x0 + 1d0)

         res = sum / k0
         print *, k0, res
      end do


      !uncomment following to obtain azimuth and polar angle
      !
      !!!for azimuth angle
      !x0 = x0 * 360.0  

      !!!cos(theta) to theta[deg]
      !y0 = (y0 * 2.0) - 1.0
      !y0 = acos(y0) * 180.0 / acos(-1.d0) 

      
      stop
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

