module tfsf_gausian
    use constants
    implicit none
    real(kind=8) :: dd,zbk
    real(kind=8) :: theta,phi,gam
    real(kind=8) :: qx,qy
    real(kind=8) :: uk

    !球座標から直交座標への変換係数
    real(kind=8) :: vxthe,vythe,vzthe
    real(kind=8) :: uxthe,uythe
    real(kind=8) :: vxphi,vyphi
    real(kind=8) :: uxphi,uyphi,uzphi
    !入射波到来方向の単位ベクトル
    real(kind=8) :: r0x,r0y
    real(kind=8) :: cogam,sigam
    real(kind=8) :: dis,vbk

    !全電磁界と散乱界の境界面
    integer :: ibd0,jbd0
    integer :: ibd1,jbd1
    integer :: ie0(2),ie1(2),je0(2),je1(2)
    integer :: ih0(2),ih1(2),jh0(2),jh1(2)
    integer :: ie2(2),ie3(2),je2(2),je3(2)
    integer :: ih2(2),ih3(2),jh2(2),jh3(2)
    contains
    subroutine init_ts_gausian()
        use constants
        implicit none
        alpha = 16.0d0/tau0/tau0
        theta0 = 90.0d0
        theta = theta0*radi0
        
        !散乱領域での速度
        vbk = c/sqrt(epsbk*mubk)
        zbk = z0*sqrt(mubk/epsbk)

        !入射角度の変換
        theta = theta0*radi0
        phi   = phi0*radi0
        gam   = gamma0*radi0

        !入射方向の単位ベクトル
        r0x = cos(phi)
        r0y = sin(phi)

        !極座標から直交座標への変換パラメータ
        vxthe = cos(theta)*cos(phi)
        vxphi =-sin(phi)
        vythe = cos(theta)*sin(phi)
        vyphi = cos(phi)
        vzthe =-sin(theta)
        uxthe =-vxphi/zbk
        uxphi = vythe/zbk
        uythe =-vyphi/zbk
        uyphi = vythe/zbk
        uzphi = vzthe/zbk

        cogam = cos(gam)
        sigam = sin(gam)

        !全電磁界領域
        ibd0 = lpml(1)+lx
        ibd1 = nx-ibd0
        jbd0 = lpml(2)+ly
        jbd1 = ny-jbd0
        if(jbndinout(jbd0)) then
            call set_Irange(ie0(1),ie1(1))
            ie2(1) = ie0(1)
            ie3(1) = ie1(1)
            if(ie1(1).eq.ibd1) then
                ie3(1) = ibd1-1
            endif
        else
            ie0(1) = 0
            ie1(1) = -1
            ie2(1) = 0
            ie3(1) = -1
        endif
        if(jbndinout(jbd1)) then
            call set_Irange(ie0(2),ie1(2))
            ie2(2) = ie0(2)
            ie3(2) = ie1(2)
            if(ie1(2).eq.ibd1) then
                ie3(2) = ibd1-1
            endif
        else
            ie0(2) = 0
            ie1(2) = -1
            ie2(2) = 0
            ie3(2) = -1
        endif
        if(ibndinout(ibd0)) then
            call set_Jrange(je0(1),je1(1))
            je2(1) = je0(1)
            je3(1) = je1(1)
            if(je1(1).eq.jbd1) then
                je3(1) = jbd1-1
            endif
        else
            je0(1) = 0
            je1(1) = -1
            je2(1) = 0
            je3(1) = -1
        endif
        if(ibndinout(ibd1)) then
            call set_Jrange(je0(2),je1(2))
            je2(2) = je0(2)
            je3(2) = je1(2)
            if(je1(2).eq.jbd1) then
                je3(2) = jbd1-1
            endif
        else
            je0(2) = 0
            je1(2) = -1
            je2(2) = 0
            je3(2) = -1
        endif

        if(jbndinout(jbd0)) then
            call set_Irange(ih0(1),ih1(1))
            ih2(1) = ih0(1)
            ih3(1) = ih1(1)
            if(ih1(1).eq.ibd1) then
                ih3(1) = ibd1-1
            endif
        else
            ih0(1) = 0
            ih1(1) = -1
            ih2(1) = 0
            ih3(1) = -1
        endif
        if(jbndinout(jbd1)) then
            call set_Irange(ih0(2),ih1(2))
            ih2(2) = ih0(2)
            ih3(2) = ih1(2)
            if(ih1(2).eq.ibd1) then
                ih3(2) = ibd1-1
            endif
        else
            ih0(2) = 0
            ih1(2) = -1
            ih2(2) = 0
            ih3(2) = -1
        endif
        if(ibndinout(ibd0)) then
            call set_Jrange(jh0(1),jh1(1))
            jh2(1) = jh0(1)
            jh3(1) = jh1(1)
            if(jh1(1).eq.jbd1) then
                jh3(1) = jbd1-1
            endif
        else
            jh0(1) = 0
            jh1(1) = -1
            jh2(1) = 0
            jh3(1) = -1
        endif
        if(ibndinout(ibd1)) then
            call set_Jrange(jh0(2),jh1(2))
            jh2(2) = jh0(2)
            jh3(2) = jh1(2)
            if(jh1(2).eq.jbd1) then
                jh3(2) = jbd1-1
            endif
        else
            jh0(2) = 0
            jh1(2) = -1
            jh2(2) = 0
            jh3(2) = -1
        endif
    
        !波頭の位置
        qx = (nx-lpml(1))*dx*r0x
        qy = (ny-lpml(2))*dy*r0y
        dis = 0.0d0
        !dd = abs(qx)
        dd = qx
        if(dd>dis) dis = dd
        !dd = abs(qy)
        dd = qy 
        if(dd>dis) dis = dd
        !dd = abs(qx+qy)
        dd = qx+qy 
        if(dd>dis) dis = dd
        if(myrank.eq.0) then
            write(30,*)"vxphi:",vxphi,"vyphi",vyphi
            write(30,*)"qx:",qx/dx,"qy:",qy/dy
            write(30,*)"dis:",dis/dx
        endif
    end subroutine

    subroutine e_add_gausian()
        use constants
        implicit none
        call ejbnd_u(ie0(1),ie1(1),ie2(1),ie3(1))
        call ejbnd_o(ie0(2),ie1(2),ie2(2),ie3(2))
        call eibnd_l(je0(1),je1(1),je2(1),je3(1))
        call eibnd_r(je0(2),je1(2),je2(2),je3(2))
    end subroutine

    subroutine h_add_gausian()
        use constants
        implicit none
        call hjbnd_u(ih0(1),ih1(1),ih2(1),ih3(1))
        call hjbnd_o(ih0(2),ih1(2),ih2(2),ih3(2))
        call hibnd_l(jh0(1),jh1(1),jh2(1),jh3(1))
        call hibnd_r(jh0(2),jh1(2),jh2(2),jh3(2))
    end subroutine

    subroutine ejbnd_u(i0,i1,i2,i3)
        use constants
        implicit none
        integer :: i,j
        integer,intent(in) :: i0,i1,i2,i3
        j = jbd0
        !$omp parallel do 
        do i=i0,i1
            ez(i,j) = ez(i,j) + bezy(i,j)*hxinc_gausian(i,j-1)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do i=i2,i3
            ex(i,j) = ex(i,j) - bexy(i,j)*hzinc_gausian(i,j-1)
        enddo
        !$omp end parallel do
    end subroutine

    subroutine ejbnd_o(i0,i1,i2,i3)
        use constants
        implicit none
        integer :: i,j 
        integer,intent(in) :: i0,i1,i2,i3
        j = jbd1
        !$omp parallel do 
        do i=i0,i1
            ez(i,j) = ez(i,j) - bezy(i,j)*hxinc_gausian(i,j)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do i=i2,i3
            ex(i,j) = ex(i,j) + bexy(i,j)*hzinc_gausian(i,j)
        enddo
        !$omp end parallel do
    end subroutine

    subroutine eibnd_l(j0,j1,j2,j3)
        use constants
        implicit none
        integer :: i,j 
        integer,intent(in) :: j0,j1,j2,j3
        i = ibd0 
        !$omp parallel do 
        do j=j0,j1
            ez(i,j) = ez(i,j) - bezx(i,j)*hyinc_gausian(i-1,j)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do j=j2,j3
            ey(i,j) = ey(i,j) + beyx(i,j)*hzinc_gausian(i-1,j)
        enddo
        !$omp end parallel do
    end subroutine

    subroutine eibnd_r(j0,j1,j2,j3)
        use constants
        implicit none
        integer :: i,j 
        integer,intent(in) :: j0,j1,j2,j3
        i = ibd1 
        !$omp parallel do 
        do j=j0,j1
            ez(i,j) = ez(i,j) + bezx(i,j)*hyinc_gausian(i,j)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do j=j2,j3
            ey(i,j) = ey(i,j) - beyx(i,j)*hzinc_gausian(i,j)
        enddo
        !$omp end parallel do
    end subroutine

    subroutine hjbnd_u(i0,i1,i2,i3)
        use constants
        implicit none
        integer :: i,j 
        integer,intent(in) :: i0,i1,i2,i3
        j=jbd0-1
        !$omp parallel do 
        do i=i0,i1
            hx(i,j) = hx(i,j) + bmxy(i,j)*ezinc_gausian(i,j+1)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do i=i2,i3
            hz(i,j) = hz(i,j) - bmzy(i,j)*exinc_gausian(i,j+1)
        enddo
        !$omp end parallel do
    end subroutine

    subroutine hjbnd_o(i0,i1,i2,i3)
        use constants
        implicit none
        integer :: i,j 
        integer,intent(in) :: i0,i1,i2,i3
        j=jbd1
        !$omp parallel do 
        do i=i0,i1
            hx(i,j) = hx(i,j) - bmxy(i,j)*ezinc_gausian(i,j)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do i=i2,i3
            hz(i,j) = hz(i,j) + bmzy(i,j)*exinc_gausian(i,j)
        enddo
        !$omp end parallel do
    end subroutine
    
    subroutine hibnd_l(j0,j1,j2,j3)
        use constants 
        implicit none
        integer :: i,j 
        integer,intent(in) :: j0,j1,j2,j3
        i=ibd0-1
        !$omp parallel do 
        do j=j0,j1
            hy(i,j) = hy(i,j) - bmyx(i,j)*ezinc_gausian(i+1,j)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do j=j2,j3
            hz(i,j) = hz(i,j) + bmzx(i,j)*eyinc_gausian(i+1,j)
        enddo
        !$omp end parallel do
    end subroutine

    subroutine hibnd_r(j0,j1,j2,j3)
        use constants 
        implicit none
        integer :: i,j 
        integer,intent(in) :: j0,j1,j2,j3
        i=ibd1
        !$omp parallel do 
        do j=j0,j1
            hy(i,j) = hy(i,j) + bmyx(i,j)*ezinc_gausian(i,j)
        enddo
        !$omp end parallel do

        !$omp parallel do 
        do j=j2,j3
            hz(i,j) = hz(i,j) - bmzx(i,j)*eyinc_gausian(i,j)
        enddo
        !$omp end parallel do
    end subroutine
    
    logical function ibndinout(ibnd)
        use constants
        implicit none
        integer,intent(in) :: ibnd
        if((ibnd.gt.istart).and.(ibnd.lt.iend)) then
            ibndinout = .true.
        else
            ibndinout = .false.
        endif
    end function

    logical function jbndinout(jbnd)
        use constants
        implicit none
        integer,intent(in) :: jbnd
        if((jbnd.gt.jstart).and.(jbnd.lt.jend)) then
            jbndinout = .true.
        else 
            jbndinout = .false.
        endif
    end function

    subroutine set_Irange(i0,i1)
        use constants
        implicit none
        integer,intent(inout) :: i0,i1
        if(istart.gt.ibd0.and.iend.lt.ibd1) then 
            i0 = istart 
            i1 = iend
        else if(istart.lt.ibd0.and.iend.gt.ibd0) then
            i0 = ibd0 
            i1 = iend 
        else if(istart.lt.ibd1.and.iend.gt.ibd1) then
            i0 = istart 
            i1 = ibd1 
        else if(iend.lt.ibd0) then
            i0 = 0
            i1 = -1
        else if(istart.gt.ibd1) then
            i0 = 0 
            i1 = -1
        endif
    end subroutine

    subroutine set_Jrange(j0,j1)
        use constants
        implicit none
        integer,intent(inout) :: j0,j1
        if(jstart.gt.jbd0.and.jend.lt.jbd1) then 
            j0 = jstart 
            j1 = jend
        else if(jstart.lt.jbd0.and.jend.gt.jbd0) then
            j0 = jbd0 
            j1 = jend 
        else if(jstart.lt.jbd1.and.jend.gt.jbd1) then
            j0 = jstart 
            j1 = jbd1 
        else if(jend.lt.jbd0) then
            j0 = 0
            j1 = -1
        else if(jstart.gt.jbd1) then
            j0 = 0 
            j1 = -1
        endif
    end subroutine 

    real(kind=8) function exinc_gausian(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = (i+0.5d0)*dx
        y = j*dy
        eth = cogam*einc_gausian(x,y)
        eph = sigam*einc_gausian(x,y)
        exinc_gausian = vxthe*eth+vxphi*eph
    end function

    real(kind=8) function eyinc_gausian(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = i*dx
        y = (j+0.5d0)*dy
        eth = cogam*einc_gausian(x,y)
        eph = sigam*einc_gausian(x,y)
        eyinc_gausian =vythe*eth+vyphi*eph
    end function

    real(kind=8) function ezinc_gausian(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth
        x = i*dx
        y = j*dy
        eth = cogam*einc_gausian(x,y)
        ezinc_gausian = vzthe*eth
    end function

    real(kind=8) function hxinc_gausian(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = i*dx
        y = (j+0.5d0)*dy
        eth = cogam*hinc_gausian(x,y)
        eph = sigam*hinc_gausian(x,y)
        hxinc_gausian = uxthe*eth+uxphi*eph
    end function

    real(kind=8) function hyinc_gausian(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eth,eph
        x = (i+0.5d0)*dx
        y = j*dy
        eth = cogam*hinc_gausian(x,y)
        eph = sigam*hinc_gausian(x,y)
        hyinc_gausian = uythe*eth+uyphi*eph
    end function

    real(kind=8) function hzinc_gausian(i,j)
        use constants
        implicit none
        integer,intent(in) :: i,j
        real(kind=8) :: x,y,eph
        x = (i+0.5d0)*dx
        y = (j+0.5d0)*dy
        eph = sigam*hinc_gausian(x,y)
        hzinc_gausian = uzphi*eph
    end function


    real(kind=8) function einc_gausian(x,y)
        use constants
        implicit none
        real(kind=8) :: tau
        real(kind=8),intent(in) :: x,y
        tau = t+(r0x*x+r0y*y-dis)/vbk
        einc_gausian = gausian(tau)
    end function

    real(kind=8) function hinc_gausian(x,y)
        use constants
        implicit none
        real(kind=8) :: tau
        real(kind=8),intent(in) :: x,y
        tau = t+(r0x*x+r0y*y-dis)/vbk
        hinc_gausian = gausian(tau)
    end function

    real(kind=8) function gausian(tau)
        use constants
        implicit none
        real(kind=8) :: tt
        real(kind=8),intent(in) :: tau
        tt = tau-tau0
        tt = tt*tt
        gausian = amp*exp(-tt*alpha)
    end function gausian

end module