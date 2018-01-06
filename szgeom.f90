module szgeom
  use precision1, only : dp
  use cosmo_params, only : model
  use szfuncs, only : S_sz,Pp_sz,Qp_sz,Sp_sz,Ppp_sz,Qpp_sz,Spp_sz
  use szlocal
  implicit none

  type four_position
      real(dp), dimension(4) :: x
      type(szlocal_class) :: szloc
  end type four_position

  type(four_position) :: pos

  private
  public :: init_four_position,del_four_position,four_position, &
            metric,tetrad,christoffel,density

contains

    function init_four_position(four_pos_arr) result(f)
        implicit none
        real(dp), dimension(4), intent(in) :: four_pos_arr
        type(four_position) :: f
        f%x=four_pos_arr
        f%szloc=init_szlocal_class(f%x(1),f%x(2))
    end function init_four_position

    subroutine del_four_position(xc)
! Pointers are allocated so need deallocate (not nullify)
! to prevent memory from becoming inaccessible
        implicit none
        type(four_position), intent(inout) :: xc
        call del_szlocal_class(xc%szloc)
    end subroutine del_four_position

    subroutine metric(x,g)
        implicit none
        real(dp), dimension(4), intent(in) :: x
        real(dp), dimension(4,4), intent(out) :: g
        type(szlocal_class) :: szloc
        real(dp) :: t,rc,theta,phi,s3,c3,s4,c4,mc3
        real(dp) :: Pp,Qp,Sp,S,R,Rp,k,ppr,p2
        real(dp) :: f1,f3,f4,f5,fp,grr,gr3,gr4
! Here rc is coordinate r not areal R
        t=x(1); rc=x(2); theta=x(3); phi=x(4)
        s3=sin(theta); s4=sin(phi)
        c3=cos(theta); c4=cos(phi)
        mc3=1._dp-c3
        Pp=Pp_sz(rc); Qp=Qp_sz(rc); Sp=Sp_sz(rc)
        S=S_sz(rc)
        szloc=init_szlocal_class(t,rc)
        call get_R(szloc)
        call get_Rprime(szloc)
        R=szloc%f%R
        Rp=szloc%f%Rp
        k=szloc%f%k
        ppr=Rp/R
        p2=R**2
        f1=(Pp*c4+Qp*s4)/S
        f3=(Qp*c4-Pp*s4)/S
        f4=f1*s3+(Sp/S)*c3
        f5=(Sp/S*s3+f1*mc3)**2+(f3*mc3)**2
        fp=ppr+f4
        grr=f5+fp**2/(1._dp-k) !metric only
        gr3=f1*mc3+(Sp/S)*s3
        gr4=mc3*s3*f3 !metric only
        g=0._dp
        g(1,1)= 1._dp
        g(2,2)=-p2*grr
        g(2,3)= p2*gr3
        g(3,2)= g(2,3)
        g(2,4)=-p2*gr4
        g(4,2)= g(2,4)
        g(3,3)=-p2
        g(4,4)= g(3,3)*s3**2
    end subroutine metric

    subroutine tetrad(x,Ec)
        implicit none
        real(dp), dimension(4), intent(in) :: x
        real(dp), dimension(4,4), intent(out) :: Ec
        type(szlocal_class) :: szloc
        real(dp) :: t,rc,theta,phi,s3,c3,s4,c4,mc3
        real(dp) :: Pp,Qp,Sp,S,R,Rp,k,ppr
        real(dp) :: f1,f3,f4,fp,gr3
! Here rc is coordinate r not areal R
        t=x(1); rc=x(2); theta=x(3); phi=x(4)
        s3=sin(theta); s4=sin(phi)
        c3=cos(theta); c4=cos(phi)        
        mc3=1._dp-c3
        Pp=Pp_sz(rc); Qp=Qp_sz(rc); Sp=Sp_sz(rc)
        S=S_sz(rc)
        szloc=init_szlocal_class(t,rc)
        call get_R(szloc)
        call get_Rprime(szloc)
        R=szloc%f%R
        Rp=szloc%f%Rp
        k=szloc%f%k
        ppr=Rp/R
        f1=(Pp*c4+Qp*s4)/S
        f3=(Qp*c4-Pp*s4)/S
        f4=f1*s3+(Sp/S)*c3
        fp=ppr+f4
        gr3=f1*mc3+(Sp/S)*s3
        Ec=0._dp
        Ec(1,1)=1._dp
! sqrt(1-k)/(R'+(S'/S)Rcos(theta))
        Ec(2,2)=sqrt(1._dp-k)/(R*fp)
        Ec(3,3)=1._dp/R
        Ec(4,4)=Ec(3,3)/s3
        Ec(3,2)=Ec(2,2)*gr3
        Ec(4,2)=-Ec(2,2)*mc3*f3/s3
    end subroutine tetrad

    subroutine christoffel(xc,gam)
        implicit none
        type(four_position), intent(inout) :: xc
        real(dp), dimension(4,4,4), intent(out) :: gam
        real(dp) :: t,rc,theta,phi,s3,c3,s4,c4,mc3
        real(dp) :: R,Rp,kp,Rd,Rdp,Rpp,mk
        real(dp) :: Pp_S,Qp_S,Sp_S,Ppp_S,Qpp_S,Spp_S
        real(dp) :: p2,ppt,ppr,pprt,pprr
        real(dp) :: fr3,f4,f5,fp
        real(dp) :: gr3,gr4,rr3,rp3,frrr
        real(dp) :: N_S,Np_S,N3_S,g23,g24,g33,g44
! Here rc is coordinate r not areal R
        t=xc%x(1); rc=xc%x(2); theta=xc%x(3); phi=xc%x(4)
! define shorthands
        s3=sin(theta); s4=sin(phi)
        c3=cos(theta); c4=cos(phi)        
        mc3=1._dp-c3
        call get_R(xc%szloc)
        call get_Rdot(xc%szloc)
        call get_kprime(xc%szloc)
        call get_Rprime(xc%szloc)
        call get_Rdotprime(xc%szloc)
        call get_Rpprime(xc%szloc)
        R=xc%szloc%f%R
        Rd=xc%szloc%f%Rd
        kp=xc%szloc%f%kp
        Rp=xc%szloc%f%Rp
        Rdp=xc%szloc%f%Rdp
        Rpp=xc%szloc%f%Rpp
        mk=1._dp-xc%szloc%f%k
! Phi(t,r) and its partial derivatives
        p2=R**2
! 1st derivatives of Phi(t,r) divided by Phi(t,r)
        ppt=Rd/R
        ppr=Rp/R
! 2nd derivatives of Phi(t,r) divided by Phi(t,r)
        pprt=Rdp/R
        pprr=Rpp/R
        Sp_S=Sp_sz(rc)/S_sz(rc)
        Qp_S=Qp_sz(rc)/S_sz(rc)
        Pp_S=Pp_sz(rc)/S_sz(rc)
        Spp_S=Spp_sz(rc)/S_sz(rc)
        Qpp_S=Qpp_sz(rc)/S_sz(rc)
        Ppp_S=Ppp_sz(rc)/S_sz(rc)
! N(r,phi) = P' cos(phi) + Q' sin(phi)
        N_S=Pp_S*c4+Qp_S*s4
! (dN/dr)/S = ( P'' cos(phi) + Q'' sin(phi) ) /S
        Np_S=Ppp_S*c4+Qpp_S*s4
! (dN/dphi)/S = ( -P' sin(phi) + Q' cos(phi) ) / S
        N3_S=Qp_S*c4-Pp_S*s4
! d/dr (dN/dphi /S) = (dN'/dphi)/S - dN/dphi S'/S^2
        fr3=Qpp_S*c4-Ppp_S*s4-N3_S*Sp_S
! (4.11) Bolejko 2016 (up to sign): N/S sin(th) + S'/S cos(th)
        f4=N_S*s3+Sp_S*c3
        fp=ppr+f4
! (4.7) Bolejko 2016: third and fourth square bracket terms
        f5=(Sp_S*s3+N_S*mc3)**2+(N3_S*mc3)**2
! g_{r\theta} / R^2
        gr3=N_S*mc3+Sp_S*s3 !all
! g_{r\phi} / R^2
        gr4=mc3*s3*N3_S
        rr3=(Np_S-N_S*Sp_S)*mc3+(Spp_S-Sp_S**2)*s3 !gamma only
        rp3=N_S*c3-Sp_S*s3 !gamma only
        frrr=0.5_dp*kp/mk &
         +(pprr-ppr**2+Np_S*s3+Spp_S*c3-Sp_S*f4 &
           -gr3*rp3+mc3*N3_S**2-f5*mk)/fp !gamma only
        g23=+p2*gr3
        g24=-p2*gr4
        g33=-p2
        g44=-p2*s3**2
        gam=0._dp
        select case (model)
            case('Szekeres','szekeres')
! Nonzero t components
                gam(1,2,2)=p2*(ppt*f5+(pprt+ppt*f4)*fp/mk)
                gam(1,2,3)=-ppt*g23
                gam(1,3,2)=gam(1,2,3)
                gam(1,3,3)=-ppt*g33
                gam(1,2,4)=-ppt*g24
                gam(1,4,2)=gam(1,2,4)
                gam(1,4,4)=-ppt*g44
! Nonzero r components
                gam(2,1,2)=(pprt+ppt*f4)/fp
                gam(2,2,1)=gam(2,1,2)
                gam(2,2,2)=ppr+ frrr
                gam(2,2,3)=(gr3*mk+rp3)/fp
                gam(2,3,2)=gam(2,2,3)
                gam(2,2,4)=(1._dp-mc3*mk)*s3*N3_S/fp
                gam(2,4,2)=gam(2,2,4)
                gam(2,3,3)=-mk/fp
                gam(2,4,4)=-mk*s3**2/fp
! Nonzero theta components
                gam(3,1,2)=(pprt-ppr*ppt)*gr3/fp
                gam(3,2,1)=gam(3,1,2)
                gam(3,1,3)=ppt
                gam(3,3,1)=gam(3,1,3)
                gam(3,2,2)=gr3*frrr -fp*(gr3+rp3/mk) -rr3 -mc3*s3*N3_S**2
                gam(3,2,3)=gr3*gam(2,2,3) + ppr
                gam(3,3,2)=gam(3,2,3)
                gam(3,2,4)=gr3*gam(2,2,4) - s3**2*N3_S
                gam(3,4,2)=gam(3,2,4)
                gam(3,3,3)=-mk*gr3/fp
                gam(3,4,4)=gr3*gam(2,4,4) - s3*c3
! Nonzero phi components
                gam(4,1,2)=N3_S*mc3/s3*(ppt-gam(2,1,2))
                gam(4,2,1)=gam(4,1,2)
                gam(4,1,4)=ppt
                gam(4,4,1)=gam(4,1,4)
                gam(4,2,2)=(mc3*(N3_S*(ppr-frrr-Sp_S) +fr3)-fp/mk*N3_S)/s3
                gam(4,2,3)=N3_S-mc3*N3_S/s3*gam(2,2,3)
                gam(4,3,2)=gam(4,2,3)
                gam(4,2,4)=ppr-mc3*N3_S/s3*gam(2,2,4)
                gam(4,4,2)=gam(4,2,4)
                gam(4,3,3)=mc3*mk*N3_S/fp/s3
                gam(4,3,4)=c3/s3
                gam(4,4,3)=gam(4,3,4)
                gam(4,4,4)=-mc3*s3*N3_S*gam(2,3,3)
            case('LTB','ltb')
                gam(3,3,1)=ppt
                gam(3,1,3)=gam(3,3,1)
                gam(4,4,1)=gam(3,3,1)
                gam(4,1,4)=gam(3,3,1)
                gam(4,4,2)=ppr
                gam(4,2,4)=gam(4,4,2)
                gam(1,4,4)=ppt*p2*s3**2
                gam(2,1,2)=pprt
                gam(2,2,1)=gam(2,1,2)
                gam(1,2,2)=p2*pprt/mk
                gam(2,2,2)=pprr+0.5_dp*kp/mk
                gam(2,4,4)=-mk*s3**2/ppr
                gam(1,3,3)=p2*ppt
                gam(2,3,3)=-mk/ppr
                gam(3,2,3)=ppr
                gam(3,3,2)=gam(3,2,3)
                gam(3,4,4)=-s3*c3
                gam(4,3,4)=c3/s3
                gam(4,4,3)=gam(4,3,4)
        end select
    end subroutine christoffel

    subroutine density(xc,rho)
! Computes kappa \times rho = 8\pi G rho
        implicit none
        type(four_position), intent(inout) :: xc
        real(dp), intent(out) :: rho
        real(dp) :: Pp,Qp,Sp,S,N,ep_e
        real(dp) :: t,rc,th,ph,R,R2,Rp,M,Mp
        t=xc%x(1); rc=xc%x(2); th=xc%x(3); ph=xc%x(4)
        S=S_sz(rc)
        Pp=Pp_sz(rc); Qp=Qp_sz(rc); Sp=Sp_sz(rc)
        N=Pp*cos(ph)+Qp*sin(ph)
! BNW 2016 (4.11): epsilon'/epsilon
        ep_e=-(N*sin(th)+Sp*cos(th))/S
        Mp=Mprime(rc)
        call get_R(xc%szloc)
        call get_Rprime(xc%szloc)
        R=xc%szloc%f%R
        R2=R**2
        Rp=xc%szloc%f%Rp
        M=xc%szloc%f%M
        rho=2._dp*(Mp-3._dp*M*ep_e)/(R2*(Rp-R*ep_e))
    end subroutine density

end module szgeom
