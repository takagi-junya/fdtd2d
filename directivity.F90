module directivity 
    use constants
    implicit none
    integer :: i0,i1,j0,j1 !積分線の位置
    integer :: emi0,emi1,emj0,emj1
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
    complex(kind=8),allocatable :: pwx(:),pwy(:),pwz(:)
    complex(kind=8),allocatable :: pux(:),puy(:),puz(:)
    
    integer :: itgstr,itgend,nintg,ndt
    real(kind=8) :: ak0,tintg,wi
    complex(kind=8),parameter :: cj=(0.0d0,1.0d0)
    complex(kind=8) :: cexpe,cexph,cofe,cofh
    contains

    subroutine init_dir()
        use constants
        implicit none
        integer :: p
        allocate(hatx(1:np),haty(1:np),hatz(1:np))
        allocate(px(1:np),py(1:np),pz(1:np))
        allocate(sx(1:np),sy(1:np),sz(1:np))
        allocate(wx(1:np),wy(1:np),wz(1:np))
        allocate(ux(1:np),uy(1:np),uz(1:np))

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

        emj0 = jstart 
        emj1 = jend 
        if((coords(1).eq.0.or.coords(1).eq.ndims(1)-1).and.coords(2).eq.0) then
            emj0 = j0
        endif
        if((coords(1).eq.0.or.coords(1).eq.ndims(1)-1).and.coords(2).eq.ndims(2)-1) then
            emj1 = j1 
        endif  

        emi0 = istart 
        emi1 = iend 
        if((coords(2).eq.0.or.coords(2).eq.ndims(2)-1).and.coords(1).eq.0) then
            emi0 = i0
        endif
        if((coords(2).eq.0.or.coords(2).eq.ndims(2)-1).and.coords(1).eq.ndims(1)-1) then
            emi1 = i1 
        endif  

        allocate(js1(emj0:emj1,4),js2(emj0:emj1,4))
        allocate(js3(emi0:emi1,4),js4(emi0:emi1,4))
        js1 = (0.0d0,0.0d0)
        js2 = (0.0d0,0.0d0)
        js3 = (0.0d0,0.0d0)
        js4 = (0.0d0,0.0d0)

        if(myrank.eq.0) then
            allocate(pwx(1:np),pwy(1:np),pwz(1:np))
            allocate(pux(1:np),puy(1:np),puz(1:np))
            allocate(wzz(1:np),uzz(1:np))
            allocate(wphi(1:np),uphi(1:np),dphi(1:np),dz(1:np))
            allocate(D(1:np))
            pwx = (0.0d0,0.0d0)
            pwy = (0.0d0,0.0d0)
            pwz = (0.0d0,0.0d0)
            pux = (0.0d0,0.0d0)
            puy = (0.0d0,0.0d0)
            puz = (0.0d0,0.0d0)
            
        endif

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
        
        if(myrank.eq.0) then
            write(30,'(a16)')"set directivity"
            write(30,'(a4,e10.3)')"k0:",ak0
            write(30,'(a8,i5)')"itgstr:",itgstr
            write(30,'(a8,i5)')"itgend:",itgend
            write(30,'(a7,e10.3,/)')"omega:",omega
        endif
    end subroutine
            
    subroutine out_dir()
        use HDF5
        use MPI
        use hdfio
        use constants
        implicit none
        integer :: p
        integer(kind=HSIZE_T),parameter :: dim3(1) = np
        if(coords(1).eq.0) then
            call sur1f
        endif
        if(coords(1).eq.ndims(1)-1) then
            call sur2f
        endif 
        if(coords(2).eq.0) then
            call sur3f
        endif
        if(coords(2).eq.ndims(2)-1) then 
            call sur4f
        endif
        call mpi_reduce(wx(1),pwx(1),np,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm2d,mpierr)
        call mpi_reduce(wy(1),pwy(1),np,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm2d,mpierr)
        call mpi_reduce(wz(1),pwz(1),np,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm2d,mpierr)
        call mpi_reduce(ux(1),pux(1),np,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm2d,mpierr)
        call mpi_reduce(uy(1),puy(1),np,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm2d,mpierr)
        call mpi_reduce(uz(1),puz(1),np,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm2d,mpierr)
        call mpi_barrier(comm2d,mpierr)
        if(myrank.eq.0) then 
            do p=1,np
                wzz(p)  = pwx(p)*sx(p)+pwy(p)*sy(p)+pwz(p)*sz(p)
                wphi(p) = pwx(p)*px(p)+pwy(p)*py(p)+pwz(p)*pz(p)
                uzz(p)  = pux(p)*sx(p)+puy(p)*sy(p)+puz(p)*sz(p)
                uphi(p) = pux(p)*px(p)+puy(p)*py(p)+puz(p)*pz(p)
                dphi(p) =-z0*wphi(p)+uzz(p)
                dz(p)   =-z0*wzz(p)-uphi(p)
                D(p) = cdabs(dphi(p))**2+cdabs(dz(p))**2
            enddo
            call hdfopen(filename(17),datasetname(7),file_id(17),group_id(17),0)
            call wrt1d(file_id(17),group_id(17),datasetname(7),dim3,D(:),istat1(17),istat2(17))
            call hdfclose(file_id(17),group_id(17),error(17))
        endif
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
            do j=emj0,emj1-1
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
            do j=emj0,emj1-1
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
            do i=emi0,emi1-1
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
            do i=emi0,emi1-1
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
            if(coords(1).eq.0) then
                call exchg2dj(hy(istart-1:iend+1,jstart-1:jend+1),lx0)
                call sur1h
            endif
            if(coords(1).eq.ndims(1)-1) then
                call exchg2dj(hy(istart-1:iend+1,jstart-1:jend+1),lx0)
                call sur2h
            endif
            if(coords(2).eq.0) then
                call exchg2di(hx(istart-1:iend+1,jstart-1:jend+1))
                call sur3h
            endif
            if(coords(2).eq.ndims(1)-1) then
                call exchg2di(hx(istart-1:iend+1,jstart-1:jend+1))
                call sur4h
            endif
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
            if(coords(1).eq.0) then
                call exchg2dj(ez(istart-1:iend+1,jstart-1:jend+1),lx0)
                call sur1e
            endif
            if(coords(1).eq.ndims(1)-1) then
                call exchg2dj(ez(istart-1:iend+1,jstart-1:jend+1),lx0)
                call sur2e
            endif
            if(coords(2).eq.0) then
                call exchg2di(ez(istart-1:iend+1,jstart-1:jend+1))
                call sur3e
            endif
            if(coords(2).eq.ndims(1)-1) then
                call exchg2di(ez(istart-1:iend+1,jstart-1:jend+1))
                call sur4e
            endif
        endif
    end subroutine

    !i=i0の電流
    subroutine sur1e()
        use constants
        implicit none
        integer :: i,j
        real(kind=8) :: eys,ezs
        i=i0
        do j=emj0,emj1-1
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
        do j=emj0,emj1-1
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
        do i=emi0,emi1-1
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
        do i=emi0,emi1-1
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
        do j=emj0,emj1-1
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
        do j=emj0,emj1-1
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
        do i=emi0,emi1-1
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
        do i=emi0,emi1-1
            hxs = 0.25d0*(hx(i,j)+hx(i,j-1)+hx(i+1,j)+hx(i+1,j-1))
            hzs = 0.50d0*(hz(i,j)+hz(i,j-1))
            js4(i,3) = js4(i,3)+hxs*cofh
            js4(i,4) = js4(i,4)+hzs*cofh
        enddo
    end subroutine

    subroutine finalize_dir()
        deallocate(hatx,haty,hatz)
        deallocate(px,py,pz)
        deallocate(sx,sy,sz)
        deallocate(wx,wy,wz)
        deallocate(ux,uy,uz)
        deallocate(js1,js2,js3,js4)
        if(myrank.eq.0) then
            deallocate(pwx,pwy,pwz)
            deallocate(pux,puy,puz)
            deallocate(wzz,uzz)
            deallocate(wphi,uphi,dphi,dz,D)
        endif
    end subroutine 
end module

    