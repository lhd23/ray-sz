module test_module
  use constants
  use szelliptic
  use szlocal
  use szfuncs
  implicit none

contains
  
    subroutine test_get_k
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Compute the age of the universe two ways: !
        ! (i) carlson symmetric forms               !
        ! (ii) splines                              !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use cosmo_params, only : ct2
        implicit none
        real(dp) :: r,k,f,fint,p,q,ul,kint,err1
        complex(dp), dimension(3) :: rts
        type(szshell) :: shell
        integer i
        open(unit=3,file='k.txt',form='formatted',status='replace')
        write(3,fmt=*) 'reference (correct) t_age [Mpc]: ', ct2
        write(3,fmt='(9X,A,6X,A,7X,A,3X,A)') &
                'r', 'k(carlson)', 'k(spline)', 'rel. error'
        do i=0,40
            r=1._dp + (200._dp-1._dp)/40 * i
            ul=r
            !new method (carlson)
            call init_shell(r,shell)
            f=sqrt(C0)*szellip_1111(ul,shell%szc)
            !NB. last saved roots does not necessarily
            !correspond to the output of zbrent
            call get_k(r,k)
            
            !spline method
            if (allocated(kspl%ddf)) then
                call splint(kspl%x,kspl%f,kspl%ddf,kspl%n,r,kint)
                p=-kint/OH2
                q=2._dp*Mfunc(r)/OH2
                call fcubic_roots(p,q,rts)
                fint=sqrt(C0)*szellip_1111(ul,rts)
                err1=abs(1._dp-fint/f)
            end if
            write(3,fmt='(f10.4,2f16.9,E13.3)') r, f, fint, err1
        end do
        close(3)
    end subroutine test_get_k
  
    subroutine test_get_R
        implicit none
        type(szshell) :: shell 
        real(dp) :: t,r,RR
        integer i
        t=-40._dp        
        write(*,*) 'Fix t =', t, ', vary r'
        do i=0,40
            r=1._dp + (200._dp-1._dp)/40 * i
            call init_shell(r,shell)
            call get_R(shell,t,RR)
            write(*,'(f16.9,f16.9)') r,RR
        end do
        r=20._dp         
        write(*,*) '*********************************'
        write(*,*) 'Fix r =', r, ', vary t'
        call init_shell(r,shell)
        do i =0,40
            t=0._dp + (-400._dp)/40 * i
            call get_R(shell,t,RR)
            write(*,'(f16.9,f16.9)') t,RR
        end do
    end subroutine test_get_R
        
    subroutine test_get_kprime
        use utils, only : dfridr
        implicit none
        type(szshell) :: shell
        real(dp) :: h,err,r,kp,kp1,err1
        integer i
        h = 0.1_dp
        open(unit=13,file='kp.txt',form='formatted',status='replace')
        write(13,'(9X,A,5X,A,7X,A,3X,A)') &
                'r','k''(numerical)','k''(carlson)','rel. error'
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            kp=dfridr(k_of_r,r,h,err)
            call init_shell(r,shell)       
            call get_kprime(shell,kp1)
            err1=abs(1._dp-kp/kp1)
            write(13,'(F10.4,E18.8,E18.8,E13.3)') r,kp,kp1,err1
        end do
        close(13)
    contains
        real(dp) function k_of_r(r)
            implicit none
            real(dp), intent(in) :: r
            call get_k(r,k_of_r)
        end function k_of_r
    end subroutine test_get_kprime
  
    subroutine test_get_Rprime
        use utils, only : dfridr
        implicit none
        type(szshell) :: shell
        real(dp) :: t,h,r,err,err1
        real(dp) :: Rprime,Rprime1
        integer i
        h=0.1_dp
        t=-100._dp
        open(unit=11,file='Rp.txt',form='formatted',status='replace')
        write(11,'(9X,A,3X,A,3X,A,3X,A)') &
                'r', 'R''(numerical)', 'R''(carlson)', 'rel. error'
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            Rprime=dfridr(R_of_r,r,h,err)
            call init_shell(r,shell)
            call get_Rprime(shell,t,Rprime1)
            err1=abs(1._dp-Rprime/Rprime1)
            write(11,fmt='(F10.4,F16.10,F14.10,E13.3)') r, Rprime, Rprime1, err1
        end do
        close(11)
    contains
        real(dp) function R_of_r(r)
            implicit none
            real(dp), intent(in) :: r
            call init_shell(r,shell)
            call get_R(shell,t,R_of_r)
        end function R_of_r
    end subroutine test_get_Rprime

    subroutine test_Rdotprime
        use utils, only : dfridr
        implicit none
        type(szshell) :: shell
        real(dp) :: t,h,r,err,err1
        real(dp) :: Rdp,Rdp1
        integer i
        h=0.1_dp
        t=-100._dp
        open(unit=14,file='Rdp.txt',form='formatted',status='replace')
        write(14,'(9X,A,3X,A,3X,A,3X,A)') &
                'r', 'Rdot''(numerical)', 'Rdot''(carlson)', 'rel. error'
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            Rdp=dfridr(Rdot_of_r,r,h,err)
            call init_shell(r,shell)
            call get_Rdotprime(shell,t,Rdp1)
            err1=abs(1._dp-Rdp/Rdp1)
            write(14,'(F10.4,E19.8,E17.8,E13.3)') r, Rdp, Rdp1, err1
        end do
        close(14)
    contains
        real(dp) function Rdot_of_r(r)
            real(dp), intent(in) :: r
            call init_shell(r,shell)
            call get_Rdot(shell,t,Rdot_of_r)
        end function Rdot_of_r
    end subroutine test_Rdotprime

    subroutine test_Hubble
        use cosmo_params, only : H0
        implicit none
        type(szshell) :: shell
        real(dp) :: t,r,H
        integer :: i,j
        open(unit=33,file='H.txt',form='formatted',status='replace')
        write(33,'(A)') 'Hubble function normalised to present H0'
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            write(33,'(f10.4)',advance='no') r
            do j=0,5
                t=-100._dp/5*j
                call init_shell(r,shell)
                call Hubble(shell,t,H)
                write(33,'(f11.6)',advance='no') H/H0
            end do
            write(33,*)
        end do
    end subroutine test_Hubble
    
    subroutine test_discriminant
        implicit none
        type(szshell) :: shell
        real(dp) :: r,p,q,Dt
        integer i
        open(unit=9,file='discrmt.txt',form='formatted',status='replace')
        write(unit=9,fmt='(9X,A,16X,A,18X,A,7X,A)') 'r','p','q','Discriminant'
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            call init_shell(r,shell)
            p=shell%szc%p
            q=shell%szc%q
            Dt=-4._dp*p**3-27._dp*q**2
            write(unit=9,fmt='(f10.4,f17.9,f19.9,E19.9)') r,p,q,Dt
        end do
        close(9)
    end subroutine test_discriminant

    subroutine test_p3
        implicit none
        type(szshell) :: shell
        real(dp) :: r
        integer i
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            call init_shell(r,shell)
            write(*,*) r,shell%szc%p3
        end do
    end subroutine test_p3

    subroutine test_formula_S9
        use utils, only : rd,gaussint
        implicit none
        real(dp), dimension(4) :: a,b !quartic
        real(dp), dimension(:,:), allocatable :: am,bm
        real(dp), dimension(size(a)) :: X,Y    
        real(dp), dimension(size(a),size(a)) :: U
        real(dp), dimension(size(a),size(a)) :: d
        real(dp), dimension(size(a),size(a)) :: d1
        real(dp) :: r,p,r2,p2 !upper and lower limits
        real(dp) :: yout,yout_integral
        complex(dp) :: a1,a2,a3
        integer :: na,nb
        a=(/0._dp,1._dp,1._dp,1._dp/)
        b=(/1._dp,1._dp,-1._dp,0._dp/)
        r=1._dp/SQRT2 !in u
        p=1._dp/SQRT3
        na=size(a)
        nb=size(b)
        if (na /= nb) STOP 'ensure a,b are same shape'
        r2=r**2 !in t
        p2=p**2
        X=sqrt(a+b*r2)
        Y=sqrt(a+b*p2)
        U=0._dp
        U(1,4)=(X(1)*X(4)*Y(2)*Y(3)+Y(1)*Y(4)*X(2)*X(3))/(r2-p2)
        allocate(am(na,1))
        allocate(bm(1,na))        
        am=reshape(a,(/na,1/))
        bm=reshape(b,(/1,na/))
        d1=matmul(am,bm)
        d=d1-transpose(d1)
        a3=dcmplx(U(1,4)**2)
        a1=a3+1._dp
        a2=a3-1._dp
        yout=THIRD*rd(a1,a2,a3)+r*p/U(1,4)
        yout_integral=gaussint(fn_S9,p,r)
        write(*,*) yout,yout_integral
    contains
        real(dp) function fn_S9(ua)
            implicit none
            real(dp), intent(in) :: ua
            fn_S9=ua**2*1._dp/sqrt(1._dp-ua**4)
        end function fn_S9
    end subroutine test_formula_S9

    subroutine test_formula_213
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Evaluate integral (2.13) in Carlson 1987.                            !
        ! This routine differs from test_formula_213 as the full, general      !
        ! formula is used and quantities have not been expanded and simplified !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use szlocal
        implicit none
        type(szshell) :: shell
        complex(dp), dimension(4) :: a,b
        complex(dp) :: yout,yout2
        real(dp) :: r,ul,ll
        integer, dimension(3) :: p
        r=15._dp
        ul=r
        ll=0._dp
        call init_shell(r,shell)
        a(1) = dcmplx(0._dp)
        a(2:4) = -shell%szc%roots
        b = (/1._dp,1._dp,1._dp,1._dp /)
        p = (/1,2,3/)
        yout = f213(a,b,ul,ll)
        yout2 = f213_2(shell,ul,p)
        write(*,*) yout
        write(*,*) yout2
        ! ul=1._dp/SQRT2
        ! ll=0._dp
        ! a = (/0._dp,1._dp,1._dp,2._dp/)
        ! b = (/1._dp,-1._dp,1._dp,1._dp/)
        ! yout = f213(a,b,ul,ll)
        ! write(*,*) yout
    end subroutine test_formula_213

    complex(dp) function f213(a,b,ulim,llim)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Use the general formula to evaluate eq (2.13) !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use utils, only : rd,rf,rf_2
        implicit none
        integer, parameter :: n=4
        complex(dp), dimension(n), intent(in) :: a,b
        real(dp), intent(in) :: ulim,llim
        complex(dp), dimension(n) :: X,Y
        complex(dp), dimension(n,1) :: am
        complex(dp), dimension(1,n) :: bm
        complex(dp), dimension(n,n) :: d1,d
        complex(dp) :: U12,U13,U14,U122,U132,U142
        complex(dp) :: s1,s2,s3,t1,t2,t3
        X = sqrt(a+b*ulim)
        Y = sqrt(a+b*llim)
        U12 = (X(1)*X(2)*Y(3)*Y(4) + Y(1)*Y(2)*X(3)*X(4)) / (ulim-llim)
        U13 = (X(1)*X(3)*Y(2)*Y(4) + Y(1)*Y(3)*X(2)*X(4)) / (ulim-llim)
        U14 = (X(1)*X(4)*Y(2)*Y(3) + Y(1)*Y(4)*X(2)*X(3)) / (ulim-llim)
        U122 = U12**2
        U132 = U13**2
        U142 = U14**2
        am = reshape(a, (/n,1/) )
        bm = reshape(b, (/1,n/) )
        d1 = matmul(am,bm)
        d = d1 - transpose(d1)
        s1 = b(1) / d(1,4) - 2._dp * b(2) / d(2,4) - 2._dp * b(3) / d(3,4)
        t1 = 2._dp / 9._dp * s1 * d(1,2) * d(1,3) * rd(U122,U132,U142)
        s2 = b(4) * d(1,2) * d(1,3) / (d(1,4) * d(2,4) * d(3,4))
        t2 = -2._dp / 3._dp * s2 * rf(U122,U132,U142)
        s3 = -2._dp / 3._dp * b(4) / (d(2,4) * d(3,4))
        t3 = s3 * X(1) * X(2) * X(3) / X(4)**3
        f213 = t1 + t2 + t3
    end function f213

    complex(dp) function f213_2(shell,ulim,perm)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Only cyclic permutations of (1,2,3) allowed         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use szelliptic
        use utils, only : rd,rf,rd_2
        implicit none
        type(szshell), intent(in) :: shell
        real(dp), intent(in) :: ulim        
        integer, dimension(3), intent(in) :: perm
        complex(dp) :: x2,x3,x4,a2,a3,a4
        complex(dp) :: s,t1,t2,t3
        real(dp) :: xi
        x2 = shell%szc%roots(perm(1))
        x3 = shell%szc%roots(perm(2))
        x4 = shell%szc%roots(perm(3))
        a2 = 1._dp/ulim - 1._dp/x2
        a3 = 1._dp/ulim - 1._dp/x3
        a4 = 1._dp/ulim - 1._dp/x4
        xi = shell%szc%xi !sqrt(-x1.x2.x3)
        s = TWOTHIRD * xi * shell%szc%f(perm(3)) / x4**2
        t1 = THIRD * (1._dp / x2 + 1._dp / x3 + 1._dp / x4 &
                      - 3._dp * x4 / (x2 * x3)) * rd(a2,a3,a4)
        t2 = rf(a2,a3,a4)
        ! t3 = -sqrt(a2 * a3 / a4**3)
        t3 = -sqrt(a2) * sqrt(a3) / sqrt(a4)**3
        f213_2 = s * (t1 + t2 + t3)
    end function f213_2

    complex(dp) function f207(shell,ulim,perm)
        use szelliptic
        use utils, only : rd
        implicit none
        real(dp), intent(in) :: ulim
        type(szshell), intent(in) :: shell
        integer, dimension(3), intent(in) :: perm
        complex(dp) :: x2,x3,x4,xi,a2,a3,a4
        x2 = shell%szc%roots(perm(1))
        x3 = shell%szc%roots(perm(2))
        x4 = shell%szc%roots(perm(3))
        xi = shell%szc%xi
        a2 = 1._dp/ulim - 1._dp/x2
        a3 = 1._dp/ulim - 1._dp/x3
        a4 = 1._dp/ulim - 1._dp/x4
        f207 = -TWOTHIRD / (x4 * xi) * rd(a2,a3,a4)
    end function f207

    subroutine test_pX555
        use utils, only : gaussint
        implicit none
        type(szshell) :: shell
        real(dp) :: yout,yout1
        real(dp) :: r,ll,ul
        reaL(dp) :: a1,a2,b2,z2
        integer i
        ll=0._dp
        open(unit=17,file='p5555.txt',form='formatted',status='replace')
        open(unit=18,file='p3555.txt',form='formatted',status='replace')
        open(unit=19,file='p1555.txt',form='formatted',status='replace')
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            ul=r
            call init_shell(r,shell)
            a1=real(shell%szc%roots(1))
            a2=real(shell%szc%roots(2))
            b2=aimag(shell%szc%roots(2))
            z2=a2**2+b2**2
! All integrals ought to be real valued
            yout=szellip_5555(shell%szc,ul)
            yout1=gaussint(f5555,ll,ul)
            write(17,fmt='(F10.4,E17.9,E17.9)') r,yout,yout1
            yout=szellip_3555(shell%szc,ul) 
            yout1=gaussint(f3555,ll,ul)
            write(18,fmt='(F10.4,E17.9,E17.9)') r,yout,yout1
            yout=szellip_1555(shell%szc,ul) 
            yout1=gaussint(f1555,ll,ul)
            write(19,fmt='(F10.4,E17.9,E17.9)') r,yout,yout1
        end do
        close(19)
        close(18)
        close(17)
! The first root is real, the others complex conjugate pair
    contains
        real(dp) function f5555(x)
            implicit none
            real(dp), intent(in) :: x
            real(dp) q
! Check branch ambiguity if complex valued
            q=(x-a1)*(x**2-2._dp*a2*x+z2)
            if (q < 0._dp) STOP 'f5555: integrand not real'
            f5555=sqrt(x/q)**5
        end function f5555
        real(dp) function f3555(x)
            implicit none
            real(dp), intent(in) :: x
            real(dp) q
! Check branch ambiguity if complex valued
            q=(x-a1)*(x**2-2._dp*a2*x+z2)
            if (q < 0._dp) STOP 'f3555: integrand not real'
            f3555=sqrt(x)**3/sqrt(q)**5
        end function f3555
        real(dp) function f1555(x)
            implicit none
            real(dp), intent(in) :: x
            real(dp) q
! Check branch ambiguity if complex valued
            q=(x-a1)*(x**2-2._dp*a2*x+z2)
            if (q < 0._dp) STOP 'f1555: integrand not real'
            f1555=sqrt(x)/sqrt(q)**5
        end function f1555        
    end subroutine test_pX555
    
    subroutine test_get_kpprime
        use utils, only : dfridr
        implicit none
        type(szshell) :: shell
        real(dp) :: h,err,r,kpp,kpp1,err1
        integer i
        h = 0.1_dp
        open(unit=13,file='kpp.txt',form='formatted',status='replace')
        write(13,'(9X,A,4X,A,6X,A,3X,A)') &
                'r','k''''(numerical)','k''''(carlson)','rel. error'
        do i=0,40
            r=1._dp+(200._dp-1._dp)/40*i
            kpp=dfridr(kp_of_r,r,h,err)
            call init_shell(r,shell)       
            call get_kpprime(shell,kpp1)
            err1=abs(1._dp-kpp/kpp1)
            write(13,'(F10.4,E18.8,E18.8,E13.3)') r,kpp,kpp1,err1
        end do
        close(13)
    contains
        real(dp) function kp_of_r(r)
            implicit none
            real(dp), intent(in) :: r
            call init_shell(r,shell)
            call get_kprime(shell,kp_of_r)
        end function kp_of_r
    end subroutine test_get_kpprime

    subroutine test_propagation
        use cosmo_params, only : r_obs,theta_obs
        use sznullgeo
        implicit none
        real(dp) :: RA,DEC,red,obs_pos(3)
        RA=0.78539816339744828_dp
        DEC=0.10210642238260403_dp
        obs_pos=(/r_obs,theta_obs,0.543*PI/)
        call propagation(obs_pos,RA,DEC,red)
        print *,'redshift', red
    end subroutine test_propagation

end module test_module
