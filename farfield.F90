module farfield
    use constants
    implicit none
    integer :: i0,i1,j0,j1  !積分面の位置
    real(kind=8) :: ic0,jc0         !閉曲線の中心
    integer :: ms,me,mf        !計算スタート、エンド時刻
    real(kind=8) :: rrmax
    real(kind=8) :: theta,phi
    real(kind=8),allocatable :: wx(:),wy(:),wz(:)
    real(kind=8),allocatable :: ux(:),uy(:),uz(:)
    real(kind=8),allocatable :: wphi(:),wzz(:),uphi(:),uzz(:)
    real(kind=8),allocatable :: dphi(:),dzz(:)
    real(kind=8) :: hatx,haty,hatz
    real(kind=8) :: sx,sy,sz
    real(kind=8) :: px,py
    real(kind=8) :: ct
    contains
    subroutine init_far()
        use constants
        implicit none
        write(30,'(a14)')"set farfield"
        !観測点角度
        theta1 = 90.0d0
        theta = theta1*radi0
        phi   = phi1*radi0 

        !積分閉曲線
        i0 = lpml(1)+isx
        i1 = nx-lpml(1)-isx
        j0 = lpml(2)+isy
        j1 = ny-lpml(2)-isy

        !閉曲線の中心
        ic0 = (i0+i1)*0.5d0
        jc0 = (j0+j1)*0.5d0

        ct = c*dt
        rrmax = sqrt(((i0-ic0)*dx)**2.0d0+((j0-jc0)*dy)**2.0d0)
        ms = intf(1.0-rrmax/ct)
        me = intf(1.5+rrmax/ct)+nstep
        mf = ms+nstep

        allocate(wx(ms:me),wy(ms:me),wz(ms:me))
        allocate(ux(ms:me),uy(ms:me),uz(ms:me))
        allocate(wzz(ms:me),wphi(ms:me),uzz(ms:me),uphi(ms:me))
        allocate(dphi(ms:me),dzz(ms:me))
        wx = 0.0d0
        wy = 0.0d0
        wz = 0.0d0
        ux = 0.0d0
        uy = 0.0d0
        uz = 0.0d0

        hatx = sin(theta)*cos(phi)
        haty = sin(theta)*sin(phi)
        hatz = cos(theta)
        sx = cos(theta)*cos(phi)
        sy = cos(theta)*sin(phi)
        sz =-sin(theta)
        px =-sin(phi)
        py = cos(phi)

        write(30,'(a20,i5,a5,i5)')"Far field start:",ms,"end:",me
        write(30,'(a30,f10.2,a5,f10.2)')"Observation  angle theta:",theta1,"phi:",phi1
        
    end subroutine

    subroutine out_far()
        use constants
        implicit none
        integer ::m
        do m=ms,me
            wphi(m)= wx(m)*px+wy(m)*py
            wzz(m) = wx(m)*sx+wy(m)*sy+wz(m)*sz
            uphi(m)= ux(m)*px+uy(m)*py
            uzz(m) = ux(m)*sx+uy(m)*sy+uz(m)*sz
            dphi(m)=-z0*wphi(m)+uzz(m)
            dzz(m) =-z0*wzz(m) -uphi(m)
        end do
        call empotential
    end subroutine out_far
    
    subroutine  empotential()
        use HDF5
        use hdfio
        use constants
        implicit none
        integer :: m
        dim2(1) = nstep+1
        call hdfopen(filename(11),datasetname(1),file_id(11),group_id(11),0)
        call wrt1d(file_id(11),group_id(11),datasetname(1),dim2,wphi(:),istat1(11),istat2(11))
        call hdfclose(file_id(11),group_id(11),error(11))

        call hdfopen(filename(12),datasetname(2),file_id(12),group_id(12),0)
        call wrt1d(file_id(12),group_id(12),datasetname(2),dim2,wzz(:),istat1(12),istat2(12))
        call hdfclose(file_id(12),group_id(12),error(12))

        call hdfopen(filename(13),datasetname(3),file_id(13),group_id(13),0)
        call wrt1d(file_id(13),group_id(13),datasetname(3),dim2,uphi(:),istat1(13),istat2(13))
        call hdfclose(file_id(13),group_id(13),error(13))

        call hdfopen(filename(14),datasetname(4),file_id(14),group_id(14),0)
        call wrt1d(file_id(14),group_id(14),datasetname(4),dim2,uzz(:),istat1(14),istat2(14))
        call hdfclose(file_id(14),group_id(14),error(14))

        call hdfopen(filename(15),datasetname(5),file_id(15),group_id(15),0)
        call wrt1d(file_id(15),group_id(15),datasetname(5),dim2,dphi(:),istat1(15),istat2(15))
        call hdfclose(file_id(15),group_id(15),error(15))

        call hdfopen(filename(16),datasetname(6),file_id(16),group_id(16),0)
        call wrt1d(file_id(16),group_id(16),datasetname(6),dim2,dzz(:),istat1(16),istat2(16))
        call hdfclose(file_id(16),group_id(16),error(16))
    end subroutine empotential
!-----------------------------------------------------------------------
!     電流の寄与
!-----------------------------------------------------------------------
    subroutine jfarfld()
        use constants
        implicit none
        call jsur1()
        call jsur2()
        call jsur3()
        call jsur4()
    end subroutine jfarfld
!----------------------------------------------------------------------
!     磁流による寄与
!----------------------------------------------------------------------
    subroutine mfarfld()
        use constants
        implicit none
        call msur1()
        call msur2()
        call msur3()
        call msur4()
    end subroutine mfarfld

    subroutine jsur1()
        use constants
        use omp_lib
        implicit none
        real(kind=8) :: ds
        real(kind=8) :: x,y
        real(kind=8) :: eta,nt,tn
        real(kind=8) :: hyavg,hzavg 
        integer :: i,j,m
        integer :: myid
        i = i0
        ds = dy
        x = (i-ic0)*dx
        nt = t/dt 

        do j=j0,j1-1
            y = (j-jc0+0.5d0)*dy
            hyavg = 0.25d0*(hy(i,j)+hy(i-1,j)+hy(i,j+1)+hy(i-1,j+1))
            hzavg = 0.50d0*(hz(i,j)+hz(i-1,j))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m
            wz(m) = wz(m)-(1.0d0-eta)*hyavg*ds
            wy(m) = wy(m)+(1.0d0-eta)*hzavg*ds
            wz(m+1)=wz(m+1)-eta*hyavg*ds
            wy(m+1)=wy(m+1)+eta*hzavg*ds
        enddo
    end subroutine

    subroutine msur1()
        use constants
        implicit none
        real(kind=8) :: ds 
        real(kind=8) :: x,y 
        real(kind=8) :: eyavg,ezavg
        real(kind=8) :: nt,eta,tn
        integer :: i,j,m

        i = i0 
        ds = dy 
        x = (i-ic0)*dx
        nt = t/dt
        do j=j0,j1-1
            y = (j-jc0+0.5d0)*dy
            eyavg = ey(i,j)
            ezavg = 0.5d0*(ez(i,j)+ez(i,j+1))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m
            uz(m) = uz(m) + (1.0d0-eta)*eyavg*ds 
            uy(m) = uy(m) - (1.0d0-eta)*ezavg*ds
            uy(m+1)=uy(m+1) - eta*ezavg*ds
            uz(m+1)=uz(m+1) + eta*eyavg*ds 
        enddo
    end subroutine

    subroutine jsur2()
        use constants
        implicit none
        real(kind=8) :: ds 
        real(kind=8) :: x,y 
        real(kind=8) ::eta,nt,tn 
        real(kind=8) :: hyavg,hzavg
        integer :: i,j,m
        i = i1
        ds = dy
        x = (i-ic0)*dx
        nt = t/dt
        do j=j0,j1-1
            y = (j-jc0+0.5d0)*dy
            hyavg = 0.25d0*(hy(i,j)+hy(i-1,j)+hy(i,j+1)+hy(i-1,j+1))
            hzavg = 0.50d0*(hz(i,j)+hz(i-1,j))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m 
            wy(m) = wy(m)-(1.0d0-eta)*hzavg*ds 
            wz(m) = wz(m)+(1.0d0-eta)*hyavg*ds
            wy(m+1)=wy(m+1)-eta*hzavg*ds 
            wz(m+1)=wz(m+1)+eta*hyavg*ds
        enddo
    end subroutine

    subroutine msur2()
        use constants
        implicit none
        real(kind=8) :: ds 
        real(kind=8) :: x,y 
        real(kind=8) :: eyavg,ezavg
        real(kind=8) :: nt,eta,tn
        integer :: i,j,m

        i = i1 
        ds = dy 
        x = (i-ic0)*dx
        nt = t/dt
        do j=j0,j1-1
            y = (j-jc0+0.5d0)*dy
            eyavg = ey(i,j)
            ezavg = 0.5d0*(ez(i,j)+ez(i,j+1))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m
            uy(m) = uy(m) + (1.0d0-eta)*ezavg*ds
            uz(m) = uz(m) - (1.0d0-eta)*eyavg*ds
            uy(m+1)=uy(m+1) + eta*ezavg*ds 
            uz(m+1)=uz(m+1) - eta*eyavg*ds
        enddo
    end subroutine

    subroutine jsur3()
        use constants
        implicit none
        real(kind=8) :: ds 
        real(kind=8) :: x,y 
        real(kind=8) :: eta,nt,tn 
        real(kind=8) :: hxavg,hzavg 
        integer :: i,j,m 
        j = j0 
        ds = dx 
        y = (j-jc0)*dy 
        nt = t/dt
        do i=i0,i1-1
            x = (i-ic0+0.5d0)*dx 
            hxavg = 0.25d0*(hx(i,j)+hx(i,j-1)+hx(i+1,j)+hx(i+1,j-1))
            hzavg = 0.50d0*(hz(i,j)+hz(i,j-1))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m  
            wz(m) = wz(m)+(1.0d0-eta)*hxavg*ds
            wx(m) = wx(m)-(1.0d0-eta)*hzavg*ds 
            wz(m+1)=wz(m+1)+eta*hxavg*ds 
            wx(m+1)=wx(m+1)-eta*hzavg*ds
        enddo
    end subroutine

    subroutine msur3()
        use constants
        implicit none
        real(kind=8) :: ds
        real(kind=8) :: x,y 
        real(kind=8) :: exavg,ezavg 
        real(kind=8) :: nt,eta,tn
        integer :: i,j,m

        j = j0 
        ds = dx 
        y = (j-jc0)*dy
        nt = t/dt
        do i=i0,i1-1
            x = (i-ic0+0.5d0)*dx
            exavg = ex(i,j)
            ezavg = 0.5d0*(ez(i,j)+ez(i+1,j))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m
            uz(m) = uz(m) - (1.0d0-eta)*exavg*ds 
            ux(m) = ux(m) + (1.0d0-eta)*ezavg*ds 
            uz(m+1) = uz(m+1) - eta*exavg*ds 
            ux(m+1) = ux(m+1) + eta*ezavg*ds 
        enddo
    end subroutine 
    
    subroutine jsur4()
        use constants
        implicit none
        real(kind=8) :: ds 
        real(kind=8) :: x,y 
        real(kind=8) :: eta,nt,tn 
        real(kind=8) :: hxavg,hzavg 
        integer :: i,j,m 
        j = j1 
        ds = dx 
        y = (j-jc0)*dy 
        nt = t/dt 
        do i=i0,i1-1 
            x = (i-ic0+0.5d0)*dx 
            hxavg = 0.25d0*(hx(i,j)+hx(i,j-1)+hx(i+1,j)+hx(i+1,j-1))
            hzavg = 0.5d0*(hz(i,j)+hz(i,j-1))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m 
            wz(m) = wz(m)-(1.0d0-eta)*hxavg*ds 
            wx(m) = wx(m)+(1.0d0-eta)*hzavg*ds 
            wz(m+1)=wz(m+1)-eta*hxavg*ds
            wx(m+1)=wx(m+1)+eta*hzavg*ds 
        enddo
    end subroutine

    subroutine msur4()
        use constants
        implicit none
        real(kind=8) :: ds
        real(kind=8) :: x,y 
        real(kind=8) :: exavg,ezavg 
        real(kind=8) :: nt,eta,tn
        integer :: i,j,m

        j = j1 
        ds = dx 
        y = (j-jc0)*dy
        nt = t/dt
        do i=i0,i1-1
            x = (i-ic0+0.5d0)*dx
            exavg = ex(i,j)
            ezavg = 0.5d0*(ez(i,j)+ez(i+1,j))
            tn = nt-(hatx*x+haty*y)/ct
            m = intf(tn)
            eta = tn-m 
            uz(m) = uz(m) + (1.0d0-eta)*exavg*ds 
            ux(m) = ux(m) - (1.0d0-eta)*ezavg*ds 
            uz(m+1) = uz(m+1) + eta*exavg*ds 
            ux(m+1) = ux(m+1) - eta*ezavg*ds 
        enddo
    end subroutine 

    integer function intf(x)
        use constants
        implicit none
        real(kind=8) :: x
        if(x>=0) then
             intf = int(x)
        else 
            intf = int(x)-1
        endif
    end function
end module