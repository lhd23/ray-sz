program test_carlson
    use precision1, only : dp
    use utils, only : rc,rd,rd_2,rj,rf,rf_r
    implicit none
    complex(dp) :: x,y,z,p
    complex(dp) :: ref
    
    x=0._dp
    y=1._dp
    z=2._dp
    p=3._dp
    ref=0.77688623778582_dp
    write(*,'(2f16.12)') real(rj(x,y,z,p)), real(ref)
    x=2._dp
    y=3._dp
    z=4._dp
    p=5._dp
    ref=0.14297579667157_dp
    write(*,'(2f16.12)') real(rj(x,y,z,p)), real(ref)
    x=dcmplx(0._dp,1._dp)
    y=dcmplx(0._dp,-1._dp)
    z=0._dp
    p=2._dp
    ref=1.6490011662711_dp
    write(*,'(2f16.12)') real(rj(x,y,z,p)), real(ref)
    x=dcmplx(-1._dp,1._dp)
    y=dcmplx(-1._dp,-1._dp)
    z=1._dp
    p=2._dp
    ref=0.94148358841220_dp
    write(*,'(2f16.12)') real(rj(x,y,z,p)), real(ref)
    x=dcmplx(0._dp,1._dp)
    y=dcmplx(0._dp,-1._dp)
    z=0._dp
    p=dcmplx(1._dp,-1._dp)
    ref=dcmplx(1.8260115229009_dp,1.2290661908643_dp)
    print *, rj(x,y,z,p), ref
    
    ! x=dcmplx(0._dp)
    ! y=dcmplx(2._dp)
    ! z=dcmplx(1._dp)
    ! print *, rd(x,y,z)
    ! x=dcmplx(2._dp)
    ! y=dcmplx(3._dp)
    ! z=dcmplx(4._dp)
    ! print *, rd(x,y,z)
    ! x=dcmplx(0._dp,1._dp)
    ! y=dcmplx(0._dp,-1._dp)
    ! z=dcmplx(2._dp)
    ! print *, rd(x,y,z)
    ! x=dcmplx(0._dp)
    ! y=dcmplx(0._dp,1._dp)
    ! z=dcmplx(0._dp,-1._dp)
    ! print *, rd(x,y,z)
    ! x=dcmplx(0._dp)
    ! y=dcmplx(-1._dp,1._dp)
    ! z=dcmplx(0._dp,1._dp)
    ! print *, rd(x,y,z)
    ! x=dcmplx(-2._dp,-1._dp)
    ! y=dcmplx(0._dp,-1._dp)
    ! z=dcmplx(-1._dp,1._dp)
    ! print *, rd(x,y,z)

    !rf
    x=1._dp
    y=2._dp
    z=0._dp
    ref=1.3110287771461_dp
    print *, rf(x,y,z), ref
    x=dcmplx(0._dp,1._dp)
    y=dcmplx(0._dp,-1._dp)
    z=dcmplx(0._dp)
    ref=1.8540746773014_dp
    print *, rf(x,y,z), ref
    x=dcmplx(-1._dp,1._dp)
    y=dcmplx(0._dp,1._dp)
    z=dcmplx(0._dp)
    ref=dcmplx(0.79612586584234_dp,-1.2138566698365_dp)
    print *, rf(x,y,z), ref
    x=2._dp
    y=3._dp
    z=4._dp
    ref=0.58408284167715_dp
    print *, rf(x,y,z), ref
    x=dcmplx(0._dp,1._dp)
    y=dcmplx(0._dp,-1._dp)
    z=dcmplx(2._dp)
    ref=1.0441445654064_dp
    print *, rf(x,y,z), ref
    x=dcmplx(-1._dp,1._dp)
    y=dcmplx(0._dp,1._dp)
    x=dcmplx(1._dp,-1._dp)
    ref=dcmplx(0.93912050218619_dp,-0.53296252018635_dp)
    print *, rf(x,y,z), ref
end program test_carlson

