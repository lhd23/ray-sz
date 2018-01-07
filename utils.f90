MODULE utils
  USE precision1, only : dp,qp
  IMPLICIT NONE

  interface rj
      module procedure rj_c,rj_r
  end interface rj

  interface rd
      module procedure rd_c,rd_r
  end interface rd
  
  interface rf
      module procedure rf_c,rf_r
  end interface rf

  interface rc
      module procedure rc_c,rc_r
  end interface rc

CONTAINS

  FUNCTION gaussint(func,a,b)
    !Use this integrator when lower bound is a singularity
    USE abscissae
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: a,b
    REAL(DP) :: gaussint
    INTERFACE
            FUNCTION func(x)
              USE precision1
              IMPLICIT NONE
              REAL(DP), INTENT(IN) :: x
              REAL(DP) :: func
            END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: m = 33
    INTEGER :: i
    REAL(DP) :: xm,xr,s,dx,f1,f2
    REAL(DP), DIMENSION(m) :: x,w

    call kwadratury(m,x,w)
    xm = 0.5_dp*(b+a)
    xr = 0.5_dp*(b-a)
    s = 0._dp

    do i = 1,m
       dx = xr*x(i)
       f1 = func(dx+xm)
       f2 = func(-dx+xm)
       s = s + w(i)*(f1 + f2)

       if (isnan(f1).or.isnan(f2)) then
          ! write(*,*) 'Warning: check integrand is real in interval!'
          gaussint = xr*s
          RETURN
       end if

    end do

    gaussint = xr*s

  END FUNCTION gaussint

  ! Taken from CAMB (version MAY 2016)
  FUNCTION rombint2(f,a,b,tol, maxit, minsteps)
    !  Rombint returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dp) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.

    ! Modified by AL to specify max iterations and minimum number of steps
    ! (min steps useful to stop wrong results on periodic or sharp functions) 
    implicit none
    integer, parameter :: MAXITER=20,MAXJ=5
    dimension g(MAXJ+1)
    real(dp) f
    external f
    real(dp) :: rombint2
    real(dp), intent(in) :: a,b,tol
    integer, intent(in):: maxit,minsteps

    integer :: nint, i, k, jmax, j
    real(dp) :: h, gmax, error, g, g0, g1, fourj

    h=0.5_dp*(b-a)
    gmax=h*(f(a)+f(b))
    g(1)=gmax
    nint=1
    error=1.0e20_dp
    i=0
    do
       i=i+1
       if (i > maxit.or.(i > 5.and.abs(error) < tol) .and. nint > minsteps) exit
       !  Calculate next trapezoidal rule approximation to integral.
       g0=0._dp
       do k=1,nint
          g0=g0+f(a+(k+k-1)*h)
       end do
       g0=0.5_dp*g(1)+h*g0
       h=0.5_dp*h
       nint=nint+nint
       jmax=min(i,MAXJ)
       fourj=1._dp
       do j=1,jmax
          !  Use Richardson extrapolation.
          fourj=4._dp*fourj
          g1=g0+(g0-g(j))/(fourj-1._dp)
          g(j)=g0
          g0=g1
       end do
       if (abs(g0).gt.tol) then
          error=1._dp-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
    end do

    rombint2=g0
    if (i > maxit .and. abs(error) > tol)  then
       write(*,*) 'Warning: Rombint2 failed to converge; '
       write (*,*)'integral, error, tol:', rombint2,error, tol
    end if
    
  END FUNCTION rombint2

  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    IMPLICIT NONE
    INTEGER :: n,nmax,i,k
    parameter (nmax=2000)
    REAL(dp) :: yp1,ypn,x(n),y(n),y2(n)
    REAL(dp) :: p,qn,sig,un,u(nmax)

    if (yp1 > .99e30) then
       y2(1)=0._dp
        u(1)=0._dp
    else
        y2(1)=-0.5_dp
        u(1)=(3._dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2._dp
        y2(i)=(sig-1._dp)/p
        u(i)=(6._dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do

    if(ypn > .99e30) then
        qn=0._dp
        un=0._dp
    else
        qn=0.5_dp
        un=(3._dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._dp)

    do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    enddo

    RETURN
  END SUBROUTINE spline

  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    IMPLICIT NONE
    INTEGER :: n,k,khi,klo
    REAL(dp) :: x,y,xa(n),y2a(n),ya(n),a,b,h

    klo=1
    khi=n
    do while (khi-klo > 1)
    k=(khi+klo)/2
    if(xa(k) > x) then
        khi=k
    else
        klo=k
    endif
    end do

    h=xa(khi)-xa(klo)

    if(h == 0._dp) stop 'bad xa input in splint'

    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ &
      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6._dp

    return

  END SUBROUTINE splint

  FUNCTION rtbis(func,x1,x2,xacc)
    !Taken from  numerical recipes. Bisection method
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc
    REAL(DP) :: rtbis
    INTERFACE
        FUNCTION func(x)
            USE precision1
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
        END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: MAXIT=200
    INTEGER :: j
    LOGICAL :: succes
    REAL(DP) :: a,b,dx,f,fmid,xmid
    a=x1
    b=x2
    fmid=func(b)
    f=func(a)
    if (f*fmid >= 0.0) then
        call zbrac(func,a,b,succes)
        if (.not. succes) STOP 'rtbis: root not bracketed'
    end if
    if (f < 0.0) then
            rtbis=a
            dx=b-a
    else
            rtbis=b
            dx=a-b
    end if
    do j=1,MAXIT
            dx=dx*0.5_dp
            xmid=rtbis+dx
            fmid=func(xmid)
            if (fmid <= 0.0) rtbis=xmid
            if (abs(dx) < xacc .or. fmid == 0.0) RETURN
    end do
    print *, 'rtbis: too many bisections'
  END FUNCTION rtbis

  SUBROUTINE zbrac(func,x1,x2,succes)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: x1,x2
    LOGICAL, INTENT(OUT) :: succes
    INTERFACE
        FUNCTION func(x)
            USE precision1
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
        END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: NTRY=50
    REAL(DP), PARAMETER :: FACTOR=1.6_dp
    INTEGER :: j
    REAL(DP) :: f1,f2
    if (x1 == x2) STOP 'zbrac: you have to guess an initial range'
    f1=func(x1)
    f2=func(x2)
    succes=.true.
    do j=1,NTRY
        if ((f1 > 0.0 .and. f2 < 0.0) .or. &
                (f1 < 0.0 .and. f2 > 0.0)) RETURN
        if (abs(f1) < abs(f2)) then
            x1=x1+FACTOR*(x1-x2)
            f1=func(x1)
        else
            x2=x2+FACTOR*(x2-x1)
            f2=func(x2)
        end if
    end do
    succes=.false.
  END SUBROUTINE zbrac

  FUNCTION rtsafe(funcd,x1,x2,xacc)
    !double precision version of nr recipe
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc
    REAL(DP) :: rtsafe
    INTERFACE
            SUBROUTINE funcd(x,fval,fderiv)
              USE precision1
              IMPLICIT NONE
              REAL(DP), INTENT(IN) :: x
              REAL(DP), INTENT(OUT) :: fval,fderiv
            END SUBROUTINE funcd
    END INTERFACE
    INTEGER, PARAMETER :: MAXIT=100
    INTEGER :: j
    REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,fl,df)
    call funcd(x2,fh,df)
    if ((fl > 0.0_dp .and. fh > 0.0_dp) .or. &
            (fl < 0.0_dp .and. fh < 0.0_dp)) &
            print *, 'root must be bracketed in rtsafe'
    if (fl == 0.0_dp) then
            rtsafe=x1
            RETURN
    else if (fh == 0.0_dp) then
            rtsafe=x2
            RETURN
    else if (fl < 0.0_dp) then
            xl=x1
            xh=x2
    else
            xh=x1
            xl=x2
    end if
    rtsafe=0.5_dp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe,f,df)
    do j=1,MAXIT
            if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0.0_dp .or. &
                 abs(2.0_dp*f) > abs(dxold*df) ) then
            dxold=dx
            dx=0.5_dp*(xh-xl)
            rtsafe=xl+dx
            if (xl == rtsafe) RETURN
       else
            dxold=dx
            dx=f/df
            temp=rtsafe
            rtsafe=rtsafe-dx
            if (temp == rtsafe) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe,f,df)
       if (f < 0.0_dp) then
            xl=rtsafe
       else
            xh=rtsafe
       end if
    end do
    print *, 'rtsafe: too many bisections'
  END FUNCTION rtsafe

  FUNCTION zbrent(func,x1,x2,tol)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,tol
    REAL(DP) :: zbrent
    INTERFACE
        FUNCTION func(x)
            USE precision1
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: x
            REAL(DP) :: func
        END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: ITMAX=100
    REAL(DP), PARAMETER :: EPS=epsilon(x1)
    INTEGER :: iter
    LOGICAL :: succes
    REAL(DP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if ((fa > 0.0_dp .and. fb > 0.0_dp) .or. (fa < 0.0_dp .and. fb < 0.0_dp)) then
        call zbrac(func,a,b,succes)
        if (.not. succes) STOP 'zbrent: root not bracketed'
    end if
    c=b
    fc=fb
    do iter=1,ITMAX
        if ((fb > 0.0_dp .and. fc > 0.0_dp) .or. (fb < 0.0_dp .and. fc < 0.0_dp)) then
            c=a
            fc=fa
            d=b-a
            e=d
        end if
        if (abs(fc) < abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        end if
        tol1=2.0_dp*EPS*abs(b)+0.5_dp*tol
        xm=0.5_dp*(c-b)
        if (abs(xm) <= tol1 .or. fb == 0.0_dp) then
            zbrent=b
            RETURN
        end if
        if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s=fb/fa
            if (a == c) then
                p=2.0_dp*xm*s
                q=1.0_dp-s
            else
                q=fa/fc
                r=fb/fc
                p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
                q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
            end if
            if (p > 0.0_dp) q=-q
            p=abs(p)
            if (2.0_dp*p  <  min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
                e=d
                d=p/q
            else
                d=xm
                e=d
            end if
        else
            d=xm
            e=d
        end if
        a=b
        fa=fb
        b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
        fb=func(b)
    end do
    print *, 'zbrent: exceeded maximum iterations'
    zbrent=b
  END FUNCTION zbrent

  subroutine fcubic_roots(p_coeff,q_coeff,roots)
    ! Find roots of depressed cubic:
    ! y^3 + py + q = 0
    use constants, only : PI,PIO3,THIRD,SQRT3
    implicit none
    real(dp), intent(in) :: p_coeff,q_coeff
    complex(dp), dimension(3), intent(out) :: roots
    real(dp) :: Dt,tht,ry1,ry2,iy2,ry3,iy3
    real(dp) :: Q,R,D,w3,S,T
! The discriminant is: -4p^3 - 27q^2
    Dt=-4._dp*p_coeff**3-27._dp*q_coeff**2
! Roots found using Cardano's formula
! https://proofwiki.org/wiki/Cardano%27s_Formula
    if (Dt < 0._dp) then
! One real root y1 and a complex conjugate pair y2 = y3*
        Q=THIRD*p_coeff
        R=-0.5_dp*q_coeff
        D=Q**3+R**2 !NOT the usual discriminant D=Dt/(-27*4)
        ! S=(R+sqrt(D))**THIRD
        ! T=(R-sqrt(D))**THIRD
        w3=R+sign(1._dp,R)*sqrt(D)
        S=cbrt(w3)
        T=-Q/S
! The roots
        ry1=S+T
        ry2=-0.5_dp*ry1
        iy2=0.5_dp*SQRT3*(S-T)
        ry3=ry2
        iy3=-iy2
    else if (Dt > 0.0_dp) then ! three distinct real roots
        iy2=0._dp
        iy3=0._dp
        tht=THIRD*acos(1.5_dp * q_coeff/p_coeff * sqrt(-3._dp/p_coeff))
        ry1=2._dp*sqrt(-p_coeff*THIRD)*cos(tht)
        ry2=2._dp*sqrt(-p_coeff*THIRD)*cos(tht-2._dp*PIO3)
        ry3=2._dp*sqrt(-p_coeff*THIRD)*cos(tht-4._dp*PIO3)
    else ! Dt = 0 ! Two real distinct roots: y1, y2 = y3  
        R=-0.5_dp*q_coeff
        S=cbrt(R)
        ry1=2._dp*S
        ry2=-S
        iy2=0._dp
        ry3=ry2
        iy3=0._dp
    end if
    roots(1)=dcmplx(ry1,0._dp)
    roots(2)=dcmplx(ry2,iy2)
    roots(3)=dcmplx(ry3,iy3)
  end subroutine fcubic_roots

  PURE FUNCTION cbrt(x)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: cbrt
    REAL(DP), PARAMETER :: THIRD=1._dp/3._dp
    cbrt=sign(abs(x)**THIRD,x)
  END FUNCTION CBRT

  SUBROUTINE laguer(a,x,its)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: its
    COMPLEX(DP), INTENT(INOUT) :: x
    COMPLEX(DP), DIMENSION(:), INTENT(IN) :: a
    REAL(DP), PARAMETER :: EPS=epsilon(1.0_dp)
    INTEGER, PARAMETER :: MR=8,MT=10,MAXIT=MT*MR
    INTEGER :: iter,m
    REAL(DP) :: abx,abp,abm,err
    COMPLEX(DP) :: dx,x1,f,g,h,sq,gp,gm,g2,abx_c
    COMPLEX(DP), DIMENSION(size(a)) :: b,d
    REAL(DP), DIMENSION(MR) :: frac = &
            (/ 0.5_dp,0.25_dp,0.75_dp,0.13_dp,0.38_dp,0.62_dp,0.88_dp,1.0_dp /)
    m=size(a)-1
    do iter=1,MAXIT
        its=iter
        abx=abs(x)
        b(m+1:1:-1)=poly_term(a(m+1:1:-1),x)
        d(m:1:-1)=poly_term(b(m+1:2:-1),x)
        f=poly(x,d(2:m))
        abx_c=dcmplx(abx)
        err=EPS*real(poly(abx_c,dcmplx(abs(b(1:m+1)))))
        if (abs(b(1)) <= err) RETURN
        g=d(1)/b(1)
        g2=g*g
        h=g2-2.0_dp*f/b(1)
        sq=sqrt((m-1)*(m*h-g2))
        gp=g+sq
        gm=g-sq
        abp=abs(gp)
        abm=abs(gm)
        if (abp < abm) gp=gm
        if (max(abp,abm) > 0.0) then
            dx=m/gp
        else
            dx=exp(dcmplx(log(1.0_dp+abx),iter))
        end if
        x1=x-dx
        if (x == x1) RETURN
        if (mod(iter,MT) /= 0) then
            x=x1
        else
            x=x-dx*frac(iter/MT)
        end if
    end do
    print *, 'laguer: too many iterations'
  END SUBROUTINE laguer

  SUBROUTINE zroots(a,roots,polish)
    IMPLICIT NONE
    COMPLEX(DP), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(DP), DIMENSION(:), INTENT(OUT) :: roots
    LOGICAL, INTENT(IN) :: polish
    REAL(DP), PARAMETER :: EPS=1.0e-8_dp
    INTEGER :: j,its,m
    INTEGER, DIMENSION(size(roots)) :: indx
    COMPLEX(DP) :: x
    COMPLEX(DP), DIMENSION(size(a)) :: ad
    m=size(a)-1
    ad(:)=a(:)
    do j=m,1,-1
        x=dcmplx(0.0_dp)
        call laguer(ad(1:j+1),x,its)
        if (abs(aimag(x)) <= 2.0_dp*EPS**2*abs(real(x))) &
                x=dcmplx(real(x))
        roots(j)=x
        ad(j:1:-1)=poly_term(ad(j+1:2:-1),x)
    end do
    if (polish) then
        do j=1,m
            call laguer(a(:),roots(j),its)
        end do
    end if
  END SUBROUTINE zroots

  FUNCTION poly(x,coeffs)
    COMPLEX(DP), INTENT(IN) :: x
    COMPLEX(DP), DIMENSION(:), INTENT(IN) :: coeffs
    COMPLEX(DP) :: poly
    COMPLEX(DP) :: pow
    COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER, PARAMETER :: NPAR_POLY=8
    INTEGER :: i,n,nn
    n=size(coeffs)
    if (n <= 0) then
        poly=0.0_dp
    else if (n < NPAR_POLY) then
        poly=coeffs(n)
        do i=n-1,1,-1
            poly=x*poly+coeffs(i)
        end do
    else
        allocate(vec(n+1))
        pow=x
        vec(1:n)=coeffs
        do
            vec(n+1)=0.0_dp
            nn=ishft(n+1,-1)
            vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
            if (nn == 1) exit
            pow=pow*pow
            n=nn
        end do
        poly=vec(1)
        deallocate(vec)
    end if
  END FUNCTION poly

  RECURSIVE FUNCTION poly_term(a,b) RESULT(u)
    COMPLEX(DP), DIMENSION(:), INTENT(IN) :: a
    COMPLEX(DP), INTENT(IN) :: b
    COMPLEX(DP), DIMENSION(size(a)) :: u
    INTEGER, PARAMETER :: NPAR_POLYTERM=8
    INTEGER :: n,j
    n=size(a)
    if (n <= 0) RETURN
    u(1)=a(1)
    if (n < NPAR_POLYTERM) then
        do j=2,n
            u(j)=a(j)+b*u(j-1)
        end do
    else
        u(2:n:2)=poly_term(a(2:n:2)+a(1:n-1:2)*b,b*b)
        u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
    end if
  END FUNCTION poly_term

  FUNCTION dfridr(func,x,h,err)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x,h
    REAL(DP), INTENT(OUT) :: err
    REAL(DP) :: dfridr
    INTERFACE
       FUNCTION func(x)
         USE precision1
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x
         REAL(DP) :: func
       END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: NTAB=10
    REAL(DP), PARAMETER :: CON=1.4_dp,CON2=CON*CON,BIG=huge(x),SAFE=2.0
    INTEGER :: ierrmin,i,j
    REAL(DP) :: hh
    REAL(DP), DIMENSION(NTAB-1) :: errt,fac
    REAL(DP), DIMENSION(NTAB,NTAB) :: a
    hh=h
    a(1,1)=(func(x+hh)-func(x-hh))/(2.0_dp*hh)
    err=BIG
    fac(1:NTAB-1)=geop(CON2,CON2,NTAB-1)
    do i=2,NTAB
       hh=hh/CON
       a(1,i)=(func(x+hh)-func(x-hh))/(2.0_dp*hh)
       do j=2,i
          a(j,i)=(a(j-1,i)*fac(j-1)-a(j-1,i-1))/(fac(j-1)-1.0_dp)
       end do
       errt(1:i-1)=max(abs(a(2:i,i)-a(1:i-1,i)),abs(a(2:i,i)-a(1:i-1,i-1)))
       ierrmin=iminloc(errt(1:i-1))
       if (errt(ierrmin) <= err) then
          err=errt(ierrmin)
          dfridr=a(1+ierrmin,i)
       end if
       if (abs(a(i,i)-a(i-1,i-1)) >= SAFE*err) RETURN
    end do
  END FUNCTION dfridr

  FUNCTION iminloc(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER, DIMENSION(1) :: imin
    INTEGER :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc

  FUNCTION geop(first,factor,n)
    REAL(DP), INTENT(IN) :: first,factor
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: geop
    INTEGER, PARAMETER :: NPAR_GEOP=4
    INTEGER, PARAMETER :: NPAR2_GEOP=2
    INTEGER :: k,k2
    REAL(DP) :: temp
    if (n > 0) geop(1)=first
    if (n <= NPAR_GEOP) then
       do k=2,n
          geop(k)=geop(k-1)*factor
       end do
    else
       do k=2,NPAR2_GEOP
          geop(k)=geop(k-1)*factor
       end do
       temp=factor**NPAR2_GEOP
       k=NPAR2_GEOP
       do
          if (k >= n) exit
          k2=k+k
          geop(k+1:min(k2,n))=temp*geop(1:min(k,n-k))
          temp=temp*temp
          k=k2
       end do
    end if
  END FUNCTION geop
  
  function rc_c(x,y)
      implicit none
      complex(dp), intent(in) :: x,y
      complex(dp) :: rc_c
      real(dp), parameter :: r=1.0e-4_dp,&
              C1=0.3_dp,C2=1._dp/7._dp,C3=3._dp/8._dp,&
              C4=9._dp/22._dp,C5=159._dp/208._dp,C6=9._dp/8._dp
      real(dp) :: Q,m1
      complex(dp) :: xm,ym,y0,w,A0,sqrtx,sqrty,lambm,Am,s
      if ((real(y) < 0.) .and. (aimag(y) == 0.)) then
          xm=x-y
          ym=-y
          w=sqrt(x/(x-y))
      else
          xm=x
          ym=y
          w=dcmplx(1._dp,0._dp)
      end if
      y0=ym
      A0=(xm+2._dp*ym)/3._dp
      Q=(r*3._dp)**(-0.125_dp)*abs(A0-xm)
      Am=A0
      m1=1._dp
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          lambm=2._dp*sqrtx*sqrty+ym
          Am=0.25_dp*(Am+lambm)
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          m1=m1*0.25_dp        
          if (m1*Q < abs(Am)) exit
      end do
      s=m1*(y0-A0)/Am
      rc_c=w*(1._dp + C1*s**2 + C2*s**3 + C3*s**4 &
              + C4*s**5 + C5*s**6 + C6*s**7)/sqrt(Am)
  end function rc_c

  function rc_r(x,y)
      implicit none
      real(dp), intent(in) :: x,y
      real(dp) :: rc_r
      real(dp), parameter :: r=1.0e-4_dp,&
              C1=0.3_dp,C2=1._dp/7._dp,C3=3._dp/8._dp,&
              C4=9._dp/22._dp,C5=159._dp/208._dp,C6=9._dp/8._dp
      real(dp) :: Q,m1
      real(dp) :: xm,ym,y0,w,A0,sqrtx,sqrty,lambm,Am,s
      if (real(y) < 0.) then
          xm=x-y
          ym=-y
          w=sqrt(x/(x-y))
      else
          xm=x
          ym=y
          w=1._dp
      end if
      y0=ym
      A0=(xm+2._dp*ym)/3._dp
      Q=(r*3._dp)**(-0.125_dp)*abs(A0-xm)
      Am=A0
      m1=1._dp
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          lambm=2._dp*sqrtx*sqrty+ym
          Am=0.25_dp*(Am+lambm)
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          m1=m1*0.25_dp        
          if (m1*Q < abs(Am)) exit
      end do
      s=m1*(y0-A0)/Am
      rc_r=w*(1._dp + C1*s**2 + C2*s**3 + C3*s**4 &
              + C4*s**5 + C5*s**6 + C6*s**7)/sqrt(Am)
  end function rc_r

  function rj_c(x,y,z,p)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Future work: implement special case when x,y,z real and nonnegative   !
      ! with at most one allowed to be zero, and real(p)<=0.                  !
      ! Note case when x,y,z are complex has not been established (not dealt  !
      ! with in arxiv:math/9409227v1)                                         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      complex(dp), intent(in) :: x,y,z,p
      complex(dp) :: rj_c
      real(dp), parameter :: r=1.0e-4_dp
      real(dp), parameter :: C1=3._dp/14._dp, C2=1._dp/6._dp, C3=9._dp/88._dp
      real(dp), parameter :: C4=3._dp/22._dp, C5=9._dp/52._dp, C6=3._dp/26._dp
      real(dp) :: Q,m1
      complex(dp) :: xm,ym,zm,pm,A0,delta,Am,lambm,dm,em
      complex(dp) :: sqrtx,sqrty,sqrtz,sqrtp
      complex(dp) :: X1,Y1,Z1,P1,XYZ,E2,E3,E4,E5
      complex(dp) :: sum1,xarg,yarg
      ! Future work: implement special case when x,y,z real and nonnegative 
      ! with at most one allowed to be zero, and real(p)<=0.
      ! Note case when x,y,z are complex has not been established (not dealt
      ! with in arxiv:math/9409227v1)
      if (real(p) <= 0.) then
          STOP 'need real(p)>0 for this implementation of Rj to work'        
      else
          xm=x
          ym=y
          zm=z
          pm=p
      end if
      A0=(xm+ym+zm+2._dp*pm)*0.2_dp
      delta=(pm-xm)*(pm-ym)*(pm-zm)
      Q=(r*0.25_dp)**(-1._dp/6._dp) &
              * max(abs(A0-xm),abs(A0-ym),abs(A0-zm),abs(A0-pm))
      Am=A0
      m1=1._dp
      sum1=0._dp
      xarg=1._dp
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          sqrtz=sqrt(zm)
          sqrtp=sqrt(pm)
          lambm=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          dm=(sqrtp+sqrtx)*(sqrtp+sqrty)*(sqrtp+sqrtz)
          em=m1**3*delta/dm**2
          yarg=1._dp+em
          sum1=sum1+m1/dm*rc(xarg,yarg)
          Am=0.25_dp*(Am+lambm)
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          zm=0.25_dp*(zm+lambm)
          pm=0.25_dp*(pm+lambm)
          m1=m1*0.25_dp
          if (m1*Q < abs(Am)) exit
      end do
      X1=(A0-x)*m1/Am
      Y1=(A0-y)*m1/Am
      Z1=(A0-z)*m1/Am
      P1=-0.5_dp*(X1+Y1+Z1)
      XYZ=X1*Y1*Z1
      e2=X1*Y1+X1*Z1+Y1*Z1-3._dp*P1**2
      e3=XYZ+2._dp*(e2+2._dp*P1**2)*P1
      e4=(2._dp*XYZ+e2*P1+3._dp*P1**3)*P1
      e5=XYZ*P1**2
      rj_c=6._dp*sum1+m1*Am**(-1.5_dp) &
              *(1._dp-e2*(C1+C5*e2-C3*E2)+C2*e3-C4*e4+C6*e5)
  end function rj_c

  function rj_r(x,y,z,p)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Future work: implement special case when x,y,z real and nonnegative   !
      ! with at most one allowed to be zero, and real(p)<=0.                  !
      ! Note case when x,y,z are complex has not been established (not dealt  !
      ! with in arxiv:math/9409227v1)                                         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      real(dp), intent(in) :: x,y,z,p
      real(dp) :: rj_r
      real(dp), parameter :: r=1.0e-4_dp
      real(dp), parameter :: C1=3._dp/14._dp, C2=1._dp/6._dp, C3=9._dp/88._dp
      real(dp), parameter :: C4=3._dp/22._dp, C5=9._dp/52._dp, C6=3._dp/26._dp
      real(dp) :: Q,m1
      real(dp) :: xm,ym,zm,pm,A0,delta,Am,lambm,dm,em
      real(dp) :: sqrtx,sqrty,sqrtz,sqrtp
      real(dp) :: X1,Y1,Z1,P1,XYZ,E2,E3,E4,E5
      real(dp) :: sum1,xarg,yarg
      ! Future work: implement special case when x,y,z real and nonnegative 
      ! with at most one allowed to be zero, and real(p)<=0.
      ! Note case when x,y,z are complex has not been established (not dealt
      ! with in arxiv:math/9409227v1)
      if (real(p) <= 0.) then
          STOP 'need real(p)>0 for this implementation of Rj to work'        
      else
          xm=x
          ym=y
          zm=z
          pm=p
      end if
      A0=(xm+ym+zm+2._dp*pm)*0.2_dp
      delta=(pm-xm)*(pm-ym)*(pm-zm)
      Q=(r*0.25_dp)**(-1._dp/6._dp) &
              * max(abs(A0-xm),abs(A0-ym),abs(A0-zm),abs(A0-pm))
      Am=A0
      m1=1._dp
      sum1=0._dp
      xarg=1._dp
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          sqrtz=sqrt(zm)
          sqrtp=sqrt(pm)
          lambm=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          dm=(sqrtp+sqrtx)*(sqrtp+sqrty)*(sqrtp+sqrtz)
          em=m1**3*delta/dm**2
          yarg=1._dp+em
          sum1=sum1+m1/dm*rc(xarg,yarg)
          Am=0.25_dp*(Am+lambm)
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          zm=0.25_dp*(zm+lambm)
          pm=0.25_dp*(pm+lambm)
          m1=m1*0.25_dp
          if (m1*Q < abs(Am)) exit
      end do
      X1=(A0-x)*m1/Am
      Y1=(A0-y)*m1/Am
      Z1=(A0-z)*m1/Am
      P1=-0.5_dp*(X1+Y1+Z1)
      XYZ=X1*Y1*Z1
      e2=X1*Y1+X1*Z1+Y1*Z1-3._dp*P1**2
      e3=XYZ+2._dp*(e2+2._dp*P1**2)*P1
      e4=(2._dp*XYZ+e2*P1+3._dp*P1**3)*P1
      e5=XYZ*P1**2
      rj_r=6._dp*sum1+m1*Am**(-1.5_dp) &
              *(1._dp-e2*(C1+C5*e2-C3*E2)+C2*e3-C4*e4+C6*e5)
  end function rj_r

  function rd_c(x,y,z)
      implicit none
      complex(dp), intent(in) :: x,y,z
      complex(dp) :: rd_c
      real(dp), parameter :: r=1.0e-4_dp
      real(dp), parameter :: C1=3._dp/14._dp, C2=1._dp/6._dp, C3=9._dp/88._dp
      real(dp), parameter :: C4=3._dp/22._dp, C5=9._dp/52._dp, C6=3._dp/26._dp
      real(dp), parameter :: THIRD=1._dp/3._dp
      real(dp) :: Q,m1
      complex(dp) :: xm,ym,zm,A0,Am,lambm
      complex(dp) :: sqrtx,sqrty,sqrtz
      complex(dp) :: X1,Y1,Z1,e2,e3,e4,e5
      complex(dp) :: sum1
      A0=(x+y+3._dp*z)*0.2_dp
      xm=x
      ym=y
      zm=z
      Q=1._dp/(r*0.25_dp)**(1._dp/6._dp) &
              * max(abs(A0-xm),abs(A0-ym),abs(A0-zm))
      Am=A0
      m1=1._dp
      sum1=0._dp
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          sqrtz=sqrt(zm)
          lambm=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          sum1=sum1+m1/(sqrtz*(zm+lambm))
          Am=0.25_dp*(Am+lambm)
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          zm=0.25_dp*(zm+lambm)
          m1=m1*0.25_dp
          if (m1*Q < abs(Am)) exit
      end do
      X1=(A0-x)*m1/Am
      Y1=(A0-y)*m1/Am
      Z1=(-X1-Y1)*THIRD
      e2=X1*Y1-6._dp*Z1**2
      e3=(3._dp*X1*Y1-8._dp*Z1**2)*Z1
      e4=3._dp*(X1*Y1-Z1**2)*Z1**2
      e5=X1*Y1*Z1**3
      rd_c=3._dp*sum1+m1*Am**(-1.5_dp) &
              *(1._dp-e2*(C1+C5*e3-C3*e2)+C2*e3-C4*e4+C6*e5)
  end function rd_c

  function rd_r(x,y,z)
      implicit none
      real(dp), intent(in) :: x,y,z
      real(dp) :: rd_r
      real(dp), parameter :: r=1.0e-4_dp
      real(dp), parameter :: C1=3._dp/14._dp, C2=1._dp/6._dp, C3=9._dp/88._dp
      real(dp), parameter :: C4=3._dp/22._dp, C5=9._dp/52._dp, C6=3._dp/26._dp
      real(dp), parameter :: THIRD=1._dp/3._dp
      real(dp) :: Q,m1
      real(dp) :: xm,ym,zm,A0,Am,lambm
      real(dp) :: sqrtx,sqrty,sqrtz
      real(dp) :: X1,Y1,Z1,e2,e3,e4,e5
      real(dp) :: sum1
      A0=(x+y+3._dp*z)*0.2_dp
      xm=x
      ym=y
      zm=z
      Q=1._dp/(r*0.25_dp)**(1._dp/6._dp) &
              * max(abs(A0-xm),abs(A0-ym),abs(A0-zm))
      Am=A0
      m1=1._dp
      sum1=0._dp
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          sqrtz=sqrt(zm)
          lambm=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          sum1=sum1+m1/(sqrtz*(zm+lambm))
          Am=0.25_dp*(Am+lambm)
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          zm=0.25_dp*(zm+lambm)
          m1=m1*0.25_dp
          if (m1*Q < abs(Am)) exit
      end do
      X1=(A0-x)*m1/Am
      Y1=(A0-y)*m1/Am
      Z1=(-X1-Y1)*THIRD
      e2=X1*Y1-6._dp*Z1**2
      e3=(3._dp*X1*Y1-8._dp*Z1**2)*Z1
      e4=3._dp*(X1*Y1-Z1**2)*Z1**2
      e5=X1*Y1*Z1**3
      rd_r=3._dp*sum1+m1*Am**(-1.5_dp) &
              *(1._dp-e2*(C1+C5*e3-C3*e2)+C2*e3-C4*e4+C6*e5)
  end function rd_r

  function rf_c(x,y,z)
      implicit none
      complex(dp), intent(in) :: x,y,z
      complex(dp) :: rf_c
      real(dp), parameter :: ERRTOL=0.0025_dp
      real(dp), parameter :: a1=1._dp/24._dp,THIRD=1._dp/3._dp
      real(dp), parameter :: C2=0.1_dp,C3=3._dp/44._dp
      real(dp), parameter :: C4=1._dp/14._dp
      complex(dp) :: xm,ym,zm,Am,lambm
      complex(dp) :: sqrtx,sqrty,sqrtz
      complex(dp) :: delx,dely,delz,e2,e3
      xm=x
      ym=y
      zm=z
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          sqrtz=sqrt(zm)
          lambm=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          zm=0.25_dp*(zm+lambm)
          Am=THIRD*(xm+ym+zm)
          if (Am == 0._dp) then
              delx=0._dp
              dely=0._dp
              delz=0._dp
          else
              delx=(Am-xm)/Am
              dely=(Am-ym)/Am
              delz=(Am-zm)/Am
          endif
          if (max(abs(delx),abs(dely),abs(delz)) < ERRTOL) exit
      end do
      e2=delx*dely-delz*delz
      e3=delx*dely*delz
      rf_c=(1._dp+(a1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(Am)
  end function rf_c

  function rf_r(x,y,z)
      implicit none
      real(dp), intent(in) :: x,y,z
      real(dp) :: rf_r
      real(dp), parameter :: ERRTOL=0.0025_dp
      real(dp), parameter :: a1=1._dp/24._dp,THIRD=1._dp/3._dp
      real(dp), parameter :: C2=0.1_dp,C3=3._dp/44._dp
      real(dp), parameter :: C4=1._dp/14._dp
      real(dp) :: xm,ym,zm,Am,lambm
      real(dp) :: sqrtx,sqrty,sqrtz
      real(dp) :: delx,dely,delz,e2,e3
      xm=x
      ym=y
      zm=z
      do
          sqrtx=sqrt(xm)
          sqrty=sqrt(ym)
          sqrtz=sqrt(zm)
          lambm=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          xm=0.25_dp*(xm+lambm)
          ym=0.25_dp*(ym+lambm)
          zm=0.25_dp*(zm+lambm)
          Am=THIRD*(xm+ym+zm)
          if (Am == 0._dp) then
              delx=0._dp
              dely=0._dp
              delz=0._dp
          else
              delx=(Am-xm)/Am
              dely=(Am-ym)/Am
              delz=(Am-zm)/Am
          endif
          if (max(abs(delx),abs(dely),abs(delz)) < ERRTOL) exit
      end do
      e2=delx*dely-delz*delz
      e3=delx*dely*delz
      rf_r=(1._dp+(a1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(Am)
  end function rf_r

  double precision FUNCTION rf_2(x,y,z)
      ! ************************************************************************
      ! *     PURPOSE: Compute Carlson fundamental integral RF
      ! *              R_F=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2)
      ! *     ARGUMENTS: Symmetric arguments x,y,z
      ! *     ROUTINES CALLED:  None.
      ! *     ALGORITHM: Due to B.C. Carlson.
      ! *     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
      ! *     REMARKS: Accepts real arguments only
      ! *     AUTHOR:  Press et al (1992).
      ! *     DATE WRITTEN:  25 Mar 91.
      ! *     REVISIONS:
      ! ***********************************************************************
      implicit none
      double precision x,y,z,ERRTOL,THIRD,a1,C2,C3,C4
      PARAMETER(ERRTOL=0.0025d0,THIRD=1.d0/3.d0, &
                    a1=1.d0/24.d0,C2=0.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      double precision alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty, &
                    sqrtz,sqrt,xt,yt,zt
      xt=x
      yt=y
      zt=z
1     continue
      sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt=0.25d0*(xt+alamb)
      yt=0.25d0*(yt+alamb)
      zt=0.25d0*(zt+alamb)
      ave=THIRD*(xt+yt+zt)
      if(ave.eq.0.d0) then
          delx=0.d0
          dely=0.d0
          delz=0.d0
      else
          delx=(ave-xt)/ave
          dely=(ave-yt)/ave
          delz=(ave-zt)/ave
      endif
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)go to 1
      e2=delx*dely-delz*delz
      e3=delx*dely*delz
      rf_2=(1.d0+(a1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
  END FUNCTION rf_2

  FUNCTION rd_2(x,y,z)
      ! ***********************************************************************
      ! *     PURPOSE: Compute Carlson degenerate integral RD
      ! *              R_D(x,y,z)=3/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2)
      ! *     ARGUMENTS: x,y,z
      ! *     ROUTINES CALLED:  None.
      ! *     ALGORITHM: Due to B.C. Carlson.
      ! *     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
      ! *     REMARKS: An earlier algorithm of Carlson 1994
      ! *     AUTHOR:  Press et al (1992)
      ! *     DATE WRITTEN:  25 Mar 91.
      ! *     REVISIONS:
      ! ***********************************************************************
      complex(dp), intent(in) :: x,y,z
      complex(dp) :: rd_2
      real(dp), parameter :: ERRTOL=0.0015_dp,a1=3._dp/14._dp,C2=1._dp/6._dp, &
              C3=9._dp/22._dp,C4=3._dp/26._dp,C5=0.25_dp*C3,C6=1.5_dp*C4
      complex(dp) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty, &
              sqrtz,sum,xt,yt,zt
      xt=x
      yt=y
      zt=z
      sum=0._dp
      fac=1._dp
      do
          sqrtx=sqrt(xt)
          sqrty=sqrt(yt)
          sqrtz=sqrt(zt)
          alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
          sum=sum+fac/(sqrtz*(zt+alamb))
          fac=0.25_dp*fac
          xt=.25_dp*(xt+alamb)
          yt=.25_dp*(yt+alamb)
          zt=.25_dp*(zt+alamb)
          ave=.2_dp*(xt+yt+3._dp*zt)
          delx=(ave-xt)/ave
          dely=(ave-yt)/ave
          delz=(ave-zt)/ave
          if(max(abs(delx),abs(dely),abs(delz)).lt.ERRTOL) exit
      end do
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6._dp*eb
      ee=ed+ec+ec
      rd_2=3._dp*sum+fac*(1._dp+ed*(-a1+C5*ed-C6*delz*ee) &
              +delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
  END FUNCTION rd_2

  subroutine dopri_stepper(ode,h,ti,to,k,yi,yo,dyi,dyo,err,atol,rtol)
      implicit none
      real(dp), intent(in) :: h,ti,atol,rtol
      real(dp), dimension(:), intent(in) :: yi,dyi
      real(dp), dimension(size(yi)), intent(out) :: yo,dyo
      real(dp), dimension(7,size(yi)), intent(out) :: k
      real(dp), dimension(size(yi)) :: y,dy,del,del0
      real(dp), intent(out) :: to,err
      real(dp) :: t
      integer :: Neq,j
      interface
          subroutine ode(t,y,dy)
              use precision1
              implicit none
              real(dp), intent(in) :: t
              real(dp), dimension(:), intent(in) :: y
              real(dp), dimension(size(y)), intent(out) :: dy
          end subroutine ode
      end interface
      real(dp), parameter :: &
              C2=0.2_dp,&
              C3=0.3_dp,&
              C4=0.8_dp,&
              C5=8._dp/9._dp,&
              A21=0.2_dp,&
              A31=3._dp/40._dp,&
              A32=9._dp/40._dp,&
              A41=44._dp/45._dp,&
              A42=-56._dp/15._dp,&
              A43=32._dp/9._dp,&
              A51=19372._dp/6561._dp,&
              A52=-25360._dp/2187._dp,&
              A53=64448._dp/6561._dp,&
              A54=-212._dp/729._dp,&
              A61=9017._dp/3168._dp,&
              A62=-355._dp/33._dp,&
              A63=46732._dp/5247._dp,&
              A64=49._dp/176._dp,&
              A65=-5103._dp/18656._dp,&
              A71=35._dp/384._dp,&
              A73=500._dp/1113._dp,&
              A74=125._dp/192._dp,&
              A75=-2187._dp/6784._dp,&
              A76=11._dp/84._dp,&
              E1=71._dp/57600._dp,&
              E3=-71._dp/16695._dp,&
              E4=71._dp/1920._dp,&
              E5=-17253._dp/339200._dp,&
              E6=22._dp/525._dp,&
              E7=-1._dp/40._dp
      Neq=size(yi)
      k(1,:)=h*dyi
! second piece      
      t=ti+C2*h
      y=yi+A21*k(1,:)
      call ode(t,y,dy)
      k(2,:)=h*dy
! third
      t=ti+C3*h
      y=yi+A31*k(1,:)+A32*k(2,:)
      call ode(t,y,dy)
      k(3,:)=h*dy
! fourth
      t=ti+C4*h
      y=yi+A41*k(1,:)+A42*k(2,:)+A43*k(3,:)
      call ode(t,y,dy)
      k(4,:)=h*dy
! fifth
      t=ti+C5*h
      y=yi+A51*k(1,:)+A52*k(2,:)+A53*k(3,:) &
              +A54*k(4,:)
      call ode(t,y,dy)
      k(5,:)=h*dy
! sixth
      t=ti+h
      y=yi+A61*k(1,:)+A62*k(2,:)+A63*k(3,:) &
              +A64*k(4,:)+A65*k(5,:)
      call ode(t,y,dy)
      k(6,:)=h*dy
! seventh and last
      to=ti+h
      yo=yi+A71*k(1,:)+A73*k(3,:)+A74*k(4,:) &
              +A75*k(5,:)+A76*k(6,:) !final
! estimate the local error 
      call ode(to,yo,dyo) !FSAL
      k(7,:)=h*dyo
      del=E1*k(1,:)+E3*k(3,:)+E4*k(4,:) &
              +E5*k(5,:)+E6*k(6,:)+E7*k(7,:)
      do j=1,Neq
          del0(j)=atol+rtol*max(abs(yo(j)),abs(yi(j)))
      end do
      err=0._dp
      do j=1,Neq
          err=err+(del(j)/del0(j))**2
      end do
      err=sqrt(err/(1._dp*Neq))
  end subroutine dopri_stepper
  
  function ydout(h,t1,tout,k,y1,y2)
      ! Dense output, error is O(h^5)
      implicit none
      real(dp) :: h,t1,tout,theta,theta1
      real(dp), dimension(:), intent(in) :: y1,y2
      real(dp), dimension(7,size(y1)), intent(in) :: k
      real(dp), dimension(size(y1)) :: con1,con2,con3,con4,con5
      real(dp), dimension(size(y1)) :: dy,ydout
      real(dp), parameter :: &
              D1=-12715105075._dp/11282082432._dp,     &
              D3=87487479700._dp/32700410799._dp,      &
              D4=-10690763975._dp/1880347072._dp,      &
              D5=701980252875._dp/199316789632._dp,    &
              D6=-1453857185._dp/822651844._dp,        &
              D7=69997945._dp/29380423._dp
      !Interpolating polynomial. Returns y between y1 and y2 at tout
      !   NOTE: our k1,k2 etc include h in them hence why
      !         our con5 is not multiplied by h while in
      !         Numerical Recipes (2007) it is.
      dy=y2-y1
      con1=y1
      con2=dy
      con3=k(1,:)-dy
      con4=dy-k(7,:)-con3
      con5=D1*k(1,:)+D3*k(3,:)+D4*k(4,:) &
              +D5*k(5,:)+D6*k(6,:)+D7*k(7,:)
      theta=(tout-t1)/h
      theta1=1._dp-theta
      ydout=con1+theta*(con2+theta1*(con3 &
              +theta*(con4+theta1*con5)))
  end function ydout

  subroutine controller(err,errold,reject,h)
      implicit none
      logical :: reject
      real(dp) :: err,errold,h,scale,safe
      real(dp) :: maxscale,minscale,alpha,beta
      safe=0.9_dp
      minscale=0.2_dp
      maxscale=10._dp
      beta=0.02_dp !0._dp
      alpha=0.2_dp - beta*0.75_dp
      if (err <= 1._dp) then
          if (err == 0._dp) then
              scale=maxscale
          else
              scale=safe*(1._dp/(err**alpha))*errold**beta
              if(scale < minscale) scale=minscale
              if(scale > maxscale) scale=maxscale
          endif
          if (reject) then
              h=h*min(scale, 1._dp)
          else
              h=h*scale
          endif
          errold=max(err, 1e-4_dp)
          reject=.false.
      else
          scale=max(safe*(1._dp/err)**alpha, minscale)
          h=h*scale
          reject=.true.
      endif
  end subroutine controller

  subroutine ang2vec(theta,phi,vector)
!***********************************************************************
!     PURPOSE: renders the vector (x,y,z) corresponding to angles
!              theta (co-latitude measured from North pole, in [0,Pi] radians)
!              and phi (longitude measured eastward, in radians)
!              North pole is (x,y,z)=(0,0,1)
!     ARGUMENTS: theta,phi
!     ROUTINES CALLED:
!     ALGORITHM:
!     ACCURACY:
!     REMARKS: adapted from Healpix library v3.30
!     AUTHOR:  Gorski et al. (2013)
!     DATE WRITTEN:  Feb 2000
!     REVISIONS:
!**********************************************************************
      use constants, only : PI
      implicit none
      real(dp), intent(in) :: theta,phi
      real(dp), intent(out), dimension(3) :: vector
      if (theta<0.0_dp .or. theta>PI) then
          print *, 'ang2vec: theta : ',theta,' is out of range [0, Pi]'
          STOP
      end if
      vector(1)=sin(theta)*cos(phi)
      vector(2)=sin(theta)*sin(phi)
      vector(3)=cos(theta)
      return
  end subroutine ang2vec

  subroutine vec2ang(vector,theta,phi)
!***********************************************************************
!     PURPOSE: renders the angles theta, phi corresponding to vector (x,y,z)
!              theta (co-latitude measured from North pole, in [0,Pi] radians)
!              and phi (longitude measured eastward, in [0,2Pi[ radians)
!              North pole is (x,y,z)=(0,0,1)
!     ARGUMENTS: vector
!     ROUTINES CALLED:
!     ALGORITHM:
!     ACCURACY:
!     REMARKS: adapted from Healpix library v3.30
!     AUTHOR:  Gorski et al. (2013)
!     DATE WRITTEN:  Feb 2000
!     REVISIONS: 2011-08
!**********************************************************************
      use constants, only : PI
      implicit none
      real(dp), intent(in), dimension(3) :: vector
      real(dp), intent(out) :: theta,phi
      theta=atan2(sqrt(vector(1)**2+vector(2)**2), vector(3))
      phi=0._dp
      if (vector(1) /= 0._dp .or. vector(2) /= 0._dp) &
              phi=atan2(vector(2),vector(1)) ! phi in ]-pi,pi]
      if (phi < 0._dp) phi=phi+2._dp*PI ! phi in [0,2pi[
      return
  end subroutine vec2ang

  subroutine rotate(vector,psi,axis)
      implicit none
      real(dp), intent(inout) :: vector(3)
      real(dp), intent(in) :: psi
      character(len=1), intent(in) :: axis
      real(dp), parameter :: ze=0._dp,on=1._dp
      real(dp), dimension(3,3) :: m1,m2,m3
      real(dp) :: c1,s1
      c1=cos(psi)
      s1=sin(psi)
      select case (axis)
      case ('x')
          m3(:,1)=(/ on, ze, ze /)
          m3(:,2)=(/ ze, c1, s1 /)
          m3(:,3)=(/ ze,-s1, c1 /)
          vector=matmul(m3,vector)
      case ('y')
          m2(:,1)=(/ c1, ze, s1 /)
          m2(:,2)=(/ ze, on, ze /)
          m2(:,3)=(/-s1, ze, c1 /)
          vector=matmul(m2,vector)
      case ('z')
          m1(:,1)=(/ c1, s1, ze /)
          m1(:,2)=(/-s1, c1, ze /)
          m1(:,3)=(/ ze, ze, on /)
          vector=matmul(m1,vector)
      case default
          STOP 'axis must be either ''x'',''y'' or ''z'''
      end select
  end subroutine rotate

END MODULE utils
