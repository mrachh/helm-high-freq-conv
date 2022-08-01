      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x(:),y(:)
      complex *16 ztmp,z2,ima
      data ima/(0.0d0,1.0d0)/
      
      n = 100
      a = -1.0d0
      b = 0.0d0

      allocate(x(4*n),y(4*n))
      h = 0.1d0
      do i=1,n
        t = -1 + 2.0d0*(i-1.0d0)/(n-1.0d0)
        x(i) = -t
        call crn_fconvgauss(t,a,b,h,y(i),der,der2)
        write(33,*) x(i),y(i)
      enddo
      do i=1,n
        ztmp = (x(i)-1.0d0) + ima*(y(i)+1.0d0)
        z2 = ima*ztmp
        x(i+n) = real(z2) - 1.0d0
        y(i+n) = imag(z2) - 1.0d0
        write(33,*) x(i+n),y(i+n)
      enddo
      do i=1,n
        ztmp = (x(i)-1.0d0) + ima*(y(i)+1.0d0)
        z2 = ima*ima*ztmp
        x(i+2*n) = real(z2) - 1.0d0
        y(i+2*n) = imag(z2) - 3.0d0
        write(33,*) x(i+2*n),y(i+2*n)
      enddo
      do i=1,n
        ztmp = (x(i)-1.0d0) + ima*(y(i)+1.0d0)
        z2 = ima*ima*ima*ztmp
        x(i+3*n) = real(z2) + 1.0d0
        y(i+3*n) = imag(z2) - 3.0d0
        write(33,*) x(i+3*n),y(i+3*n)
      enddo
      
      
      
      
      stop
      end


c
c
c
c
c
        subroutine crn_fconvgauss(x, a, b, h, val, der, der2)
        implicit real *8 (a-h,o-z)
c
c       this routine computes the convolution
c
c         ( a*abs(x)+b ) \star 1/(sqrt(2*pi)*h^2) exp(-x^2/(2*h^2))
c
c       this effectively smoothes off the corner from the abs(x) function
c
c
        done=1
        two=2
        pi=4*atan(done)
c
c       . . . formulas are computed via maple
c
        x2=x/sqrt(two)/h
        call qerrfun(x2,verf)
        val=a*x*verf+b+sqrt(two/pi)*a*h*exp(-x*x/two/h/h)
c
        fnorm=1/sqrt(2*pi)*exp(-x*x/2)
        der=a*verf
        der2=a*sqrt(two/pi)/h*exp(-x*x/two/h/h)
c
        return
        end
c
c
c
c
c
