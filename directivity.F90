module directivity 
    use constants
    implicit none
    integer :: i0,i1,j0,j1 !積分線の位置
    real(kind=8) :: ic0,jc0 !積分線の中心

    integer,parameter :: np=73
    real(kind=8) :: theta,phi
    real(kind=8),allocatable :: hatx(:),haty(:),hatz(:)
    real(kind=8),allocatable :: px(:),py(:),pz(:),sx(:),sy(:),sz(:),D(:)
    
    complex(kind=8),allocatable :: js1(:,:),js2(:,:),js3(:,:),js4(:,:)
    complex(kind=8),allocatable :: wx(:),wy(:),wz(:),wzz(:)
    complex(kind=8),allocatable :: ux(:),uy(:),uz(:),uzz(:)
    complex(kind=8),allocatable :: wphi(:),uphi(:)
    complex(kind=8),allocatable :: dphi(:),dz(:)
    
    integer :: itgstr,itgend,nintg,ndt
    real(kind=8) :: ak0,tintg,wi
    complex(kind=8),parameter :: cj=(0.0d0,1.0d0)
    complex(kind=8) :: cexpe,cexph,cofe,cofh
    contains

    subroutine init_dir()
        use constants
        implicit none
        integer ::test
        integer :: p
        write(30,'(a16)')"set directivity"
        allocate(hatx(np),haty(np),hatz(np))
        allocate(px(np),py(np),pz(np))
        allocate(sx(np),sy(np),sz(np))
        allocate(wx(np),wy(np),wz(np),wzz(np))
        allocate(ux(np),uy(np),uz(np),uzz(np))
        allocate(wphi(np),uphi(np),dphi(np),dz(np))
        allocate(D(np))

        wx = (0.0d0,0.0d0)
        wy = (0.0d0,0.0d0)
        wz = (0.0d0,0.0d0)
        
        ux = (0.0d0,0.0d0)
        uy = (0.0d0,0.0d0)
        uz = (0.0d0,0.0d0)
        

        !閉曲線の位置
        i0 = lpml(1)+isx 
        i1 = nx-lpml(1)-isx
        j0 = lpml(2)+isy 
        j1 = ny-lpml(2)-isy
        !閉曲線の中心
        ic0 = (i0+i1)*0.5d0
        jc0 = (j0+j1)*0.5d0
        allocate(js1(j0:j1,4),js2(j0:j1,4))
        allocate(js3(i0:i1,4),js4(i0:i1,4))
        js1 = (0.0d0,0.0d0)
        js2 = (0.0d0,0.0d0)
        js3 = (0.0d0,0.0d0)
        js4 = (0.0d0,0.0d0)

        !角度のラジアン変換
        theta1 = 90.0d0
        theta  = theta1*radi0
        do p=1,np
            phi = 5.0d0*(p-1)*radi0
            
            hatx(p) = sin(theta)*cos(phi)
            haty(p) = sin(theta)*sin(phi)
            hatz(p) = cos(theta)

            sx(p)   = cos(theta)*cos(phi)
            sy(p)   = cos(theta)*sin(phi)
            sz(p)   =-sin(theta)
            
            px(p)   =-sin(phi)
            py(p)   = cos(phi)
            pz(p)   = 0.0d0
        enddo

        !積分周期
        nintg = 3
        !積分時間
        tintg = nintg/freq
        wi = 0.333333333d0*dt/tintg
        !積分step数
        ndt = int(tintg/dt)+1
        if(mod(ndt,2).ne.0.0d0) ndt = ndt + 1
        !波数
        ak0 = 2.0d0*pai/lambda
        !積分スタート
        itgstr = nstep-ndt
        !積分終了
        itgend = nstep
        write(30,'(a4,e10.3)')"k0:",ak0
        write(30,'(a8,i5)')"itgstr:",itgstr
        write(30,'(a8,i5)')"itgend:",itgend
        write(30,'(a7,e10.3,/)')"omega:",omega
        
    end subroutine
            
    subroutine out_dir()
        use HDF5
        use hdfio
        use constants
        implicit none
        integer :: p
        integer(kind=HSIZE_T),parameter :: dim3(1) = np
        call sur1f
        call sur2f 
        call sur3f 
        call sur4f
        do p=1,np
            wzz(p)  = wx(p)*sx(p)+wy(p)*sy(p)+wz(p)*sz(p)
            wphi(p) = wx(p)*px(p)+wy(p)*py(p)+wz(p)*pz(p)
            uzz(p)  = ux(p)*sx(p)+uy(p)*sy(p)+uz(p)*sz(p)
            uphi(p) = ux(p)*px(p)+uy(p)*py(p)+uz(p)*pz(p)
            dphi(p) =-z0*wphi(p)+uzz(p)
            dz(p)   =-z0*wzz(p)-uphi(p)
            D(p) = cdabs(dphi(p))**2+cdabs(dz(p))**2
        enddo

        call hdfopen(filename(17),datasetname(7),file_id(17),group_id(17),0)
        call wrt1d(file_id(17),group_id(17),datasetname(7),dim3,D(:),istat1(17),istat2(17))
        call hdfclose(file_id(17),group_id(17),error(17))
    end subroutine 

    !i=i0面の遠方界
    subroutine sur1f()
        use constants
        implicit none
        integer :: i,j,p
        real(kind=8) :: x,y
        real(kind=8) :: dl,rr
        complex(kind=8) :: ss

        i=i0
        dl = dy 
        x = (i-ic0)*dx
        do p=1,np
            do j=j0,j1-1
                y=(j-jc0+0.5d0)*dy
                rr = hatx(p)*x+haty(p)*y
                ss = cdexp(cj*ak0*rr)*dl
                wy(p) = wy(p) + js1(j,4)*ss
                wz(p) = wz(p) - js1(j,3)*ss
                uy(p) = uy(p) - js1(j,2)*ss
                uz(p) = uz(p) + js1(j,1)*ss
            enddo
        enddo
    end subroutine

    !i=i1面の遠方界
    subroutine sur2f()
        use constants
        implicit none
        integer :: i,j,p
        real(kind=8) :: x,y
        real(kind=8) :: dl,rr
        complex(kind=8) :: ss

        i=i1
        dl = dy 
        x = (i-ic0)*dx
        do p=1,np
            do j=j0,j1-1
                y=(j-jc0+0.5d0)*dy
                rr = hatx(p)*x+haty(p)*y
                ss = cdexp(cj*ak0*rr)*dl
                wy(p) = wy(p) - js2(j,4)*ss
                wz(p) = wz(p) + js2(j,3)*ss
                uy(p) = uy(p) + js2(j,2)*ss
                uz(p) = uz(p) - js2(j,1)*ss
            enddo
        enddo
    end subroutine

    !j=j0面の遠方界
    subroutine sur3f()
        use constants
        implicit none
        integer :: i,j,p
        real(kind=8) :: x,y
        real(kind=8) :: dl,rr
        complex(kind=8) :: ss

        j=j0
        dl = dx 
        y = (j-jc0)*dy
        do p=1,np
            do i=i0,i1-1
                x=(i-ic0+0.5d0)*dx
                rr = hatx(p)*x+haty(p)*y
                ss = cdexp(cj*ak0*rr)*dl
                wx(p) = wx(p) - js3(i,4)*ss
                wz(p) = wz(p) + js3(i,3)*ss
                ux(p) = ux(p) + js3(i,2)*ss
                uz(p) = uz(p) - js3(i,1)*ss
            enddo
        enddo
    end subroutine

    !j=j1面の遠方界
    subroutine sur4f()
        use constants
        implicit none
        integer :: i,j,p
        real(kind=8) :: x,y
        real(kind=8) :: dl,rr
        complex(kind=8) :: ss

        j=j1
        dl = dx
        y = (j-jc0)*dy
        do p=1,np
            do i=i0,i1-1
                x=(i-ic0+0.5d0)*dx
                rr = hatx(p)*x+haty(p)*y
                ss = cdexp(cj*ak0*rr)*dl
                wx(p) = wx(p) + js4(i,4)*ss
                wz(p) = wz(p) - js4(i,3)*ss
                ux(p) = ux(p) - js4(i,2)*ss
                uz(p) = uz(p) + js4(i,1)*ss
            enddo
        enddo
    end subroutine

    subroutine k_e()
        use constants
        implicit none
        cexph = cdexp(-cj*omega*t)
        if(step.ge.itgstr.and.step.le.itgend) then
            if(step.eq.itgstr.or.step.eq.itgend) then
                cofh = wi*cexph
            else if(mod(step-itgstr,2).eq.0.0d0) then
                cofh = 2.0d0*wi*cexph
            else
                cofh = 4.0d0*wi*cexph
            endif
            call sur1h
            call sur2h
            call sur3h
            call sur4h
        endif
    end subroutine

    subroutine k_m()
        use constants
        implicit none
        cexpe = cdexp(-cj*omega*t)
        if((step.ge.itgstr).and.(step.le.itgend)) then
            if((step.eq.itgstr).or.(step.eq.itgend)) then
                cofe = wi*cexpe
            else if(mod(step-itgstr,2).eq.0.0d0) then
                cofe = 2.0d0*wi*cexpe
            else
                cofe = 4.0d0*wi*cexpe
            endif
            call sur1e
            call sur2e
            call sur3e
            call sur4e
        endif
    end subroutine

    !i=i0の電流
    subroutine sur1e()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: eys,ezs
        i=i0
        do j=j0,j1-1
            eys = ey(i,j)
            ezs = 0.5d0*(ez(i,j)+ez(i,j+1))
            js1(j,1) = js1(j,1)+eys*cofe
            js1(j,2) = js1(j,2)+ezs*cofe
        enddo
    endsubroutine
    
    !i=i1の電流
    subroutine sur2e()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: eys,ezs
        i=i1
        do j=j0,j1-1
            eys = ey(i,j)
            ezs = 0.5d0*(ez(i,j)+ez(i,j+1))
            js2(j,1) = js2(j,1)+eys*cofe
            js2(j,2) = js2(j,2)+ezs*cofe
        enddo
    end subroutine

    !j=j0の電流
    subroutine sur3e()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: exs,ezs
        j=j0
        do i=i0,i1-1
            exs = ex(i,j)
            ezs = 0.5d0*(ez(i,j)+ez(i+1,j))
            js3(i,1) = js3(i,1)+exs*cofe
            js3(i,2) = js3(i,2)+ezs*cofe
        enddo
    endsubroutine

    !j=j1の電流
    subroutine sur4e()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: exs,ezs
        j=j1
        do i=i0,i1-1
            exs = ex(i,j)
            ezs = 0.5d0*(ez(i,j)+ez(i+1,j))
            js4(i,1) = js4(i,1)+exs*cofe
            js4(i,2) = js4(i,2)+ezs*cofe
        enddo
    endsubroutine

    !i=i0の磁流
    subroutine sur1h()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: hys,hzs
        i=i0
        do j=j0,j1-1
            hys = 0.25d0*(hy(i,j)+hy(i-1,j)+hy(i,j+1)+hy(i-1,j+1))
            hzs = 0.50d0*(hz(i,j)+hz(i-1,j))
            js1(j,3) = js1(j,3)+hys*cofh
            js1(j,4) = js1(j,4)+hzs*cofh
        enddo
    endsubroutine

    !i=i1の磁流
    subroutine sur2h()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: hys,hzs
        i=i1
        do j=j0,j1-1
            hys = 0.25d0*(hy(i,j)+hy(i-1,j)+hy(i,j+1)+hy(i-1,j+1))
            hzs = 0.50d0*(hz(i,j)+hz(i-1,j))
            js2(j,3) = js2(j,3)+hys*cofh
            js2(j,4) = js2(j,4)+hzs*cofh
        enddo
    endsubroutine

    !j=j0の磁流
    subroutine sur3h()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: hxs,hzs
        j=j0
        do i=i0,i1-1
            hxs = 0.25d0*(hx(i,j)+hx(i,j-1)+hx(i+1,j)+hx(i+1,j-1))
            hzs = 0.50d0*(hz(i,j)+hz(i,j-1))
            js3(i,3) = js3(i,3)+hxs*cofh
            js3(i,4) = js3(i,4)+hzs*cofh
        enddo
    endsubroutine
    
    !j=j1の磁流
    subroutine sur4h()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: hxs,hzs
        j=j1
        do i=i0,i1-1
            hxs = 0.25d0*(hx(i,j)+hx(i,j-1)+hx(i+1,j)+hx(i+1,j-1))
            hzs = 0.50d0*(hz(i,j)+hz(i,j-1))
            js4(i,3) = js4(i,3)+hxs*cofh
            js4(i,4) = js4(i,4)+hzs*cofh
        enddo
    end subroutine
end module

    