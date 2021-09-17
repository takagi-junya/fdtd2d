subroutine setup()
    use HDF5
    use constants
    use hdfio
    implicit none
    integer :: i,j,k

    stride = 1
    ostart=0
    erad = 0.0
    odom = (/0,0,0,0/)

    namelist /space/nxx,nyy,dx,dy,abc,pbc,lpml
    namelist /time/deltat,nstep
    namelist /output/out,ostart,odom,stride,comp,io,jo
    namelist /scatt/mode,lx,ly,gamma0,theta0,phi0,amp,freq,tau0!lambda,tau0
    namelist /far/isx,isy,theta1,phi1
    namelist /object/obj,med,ic,jc,lx2,ly2,epsr,radius
    namelist /wave/kwave,amps,orgs,ang,pt,pw
    namelist /plasma/pls,prad,erad,nu,wp

    open(10,file="param.inp",action="read")
    read(10,nml=space)
    read(10,nml=time)
    read(10,nml=output)
    read(10,nml=scatt)
    read(10,nml=far)
    read(10,nml=object)
    read(10,nml=wave)    
    read(10,nml=plasma)
    close(10)

    dx = dble(dx)
    dy = dble(dy)
    deltat=dble(deltat)
    gamma0 = dble(gamma0)
    theta0 = dble(theta0)
    phi0 = dble(phi0)
    amp = dble(amp)
    !lambda = dble(lambda)*dx
    freq = dble(freq)
    tau0=dble(tau0)
    theta1=dble(theta1)
    phi1=dble(phi1)
    epsr = dble(epsr)
    radius = dble(radius)
    orgs = dble(orgs)
    prad = dble(prad)
    erad = dble(erad)
    nu = dble(nu)
    wp = dble(wp)

    !吸収境界
    if(abc==0) then
        nx = nxx
        ny = nyy
    else
        nx = nxx+2*lpml(1)
        ny = nyy+2*lpml(2)
    endif
    ic = int(0.5*nx)
    jc = int(0.5*ny)
    jo = jc

    write(30,'(a4,i4,a4,i4)')"nx:",nx,"ny:",ny
    write(30,'(a4,i4,a4,i4,/)')"ic:",ic,"jc:",jc

    allocate(ex(0:nx,0:ny),ey(0:nx,0:ny),ez(0:nx,0:ny))
    allocate(hx(0:nx,0:ny),hy(0:nx,0:ny),hz(0:nx,0:ny))
    allocate(jx(0:nx,0:ny),jy(0:nx,0:ny),jz(0:nx,0:ny))
    
    allocate(nd(0:nx,0:ny))

    allocate(aex(0:nx,0:ny),aey(0:nx,0:ny),aez(0:nx,0:ny))
    allocate(bexy(0:nx,0:ny),beyx(0:nx,0:ny))
    allocate(bezx(0:nx,0:ny),bezy(0:nx,0:ny))
    
    allocate(amx(0:nx,0:ny),amy(0:nx,0:ny),amz(0:nx,0:ny))
    allocate(bmxy(0:nx,0:ny),bmyx(0:nx,0:ny))
    allocate(bmzx(0:nx,0:ny),bmzy(0:nx,0:ny))
    
    
    allocate(ajx(0:nx,0:ny),ajy(0:nx,0:ny),ajz(0:nx,0:ny))
    allocate(ajex(0:nx,0:ny),ajey(0:nx,0:ny),ajez(0:nx,0:ny))

    allocate(epsd(-1:nx,-1:ny),sgmed(-1:nx,-1:ny))
    allocate(mud(-1:nx,-1:ny),sgmmd(-1:nx,-1:ny))

    if(pls.ge.5) then
        allocate(vx(0:nx,0:ny),vy(0:nx,0:ny),vz(0:nx,0:ny))
        allocate(avx(0:nx,0:ny),avy(0:nx,0:ny),avz(0:nx,0:ny))
        avx = 0.0d0
        avy = 0.0d0
        avz = 0.0d0
    endif
    
    filename(1:6) = (/"ex.h5","ey.h5","ez.h5","hx.h5","hy.h5","hz.h5"/)
    filename(7:10) = (/"jx.h5","jy.h5","jz.h5","nd.h5"/)
    filename(11:14) = (/"wphi.h5","wz.h5","uphi.h5","uz.h5"/)
    filename(15:17)=(/"dphi.h5","dz.h5","D.h5"/)
    filename(18:20)=(/"vx.h5","vy.h5","vz.h5"/)

    groupname(1:6) = (/"ex","ey","ez","hx","hy","hz"/)
    groupname(7:10) = (/"jx","jy","jz","nd"/)
    groupname(11:14)=(/"wphi","wz","uphi","uz"/)
    groupname(15:17)=(/"dphi","dz","D"/)
    groupname(18:20)=(/"vx","vy","vz"/)

    datasetname(1:4)=(/"wphi","wz","uphi","uz"/)
    datasetname(5:7)=(/"dphi","dz","D"/)

    dt=deltat/(c*sqrt(1.0d0/(dx*dx)+1.0d0/(dy*dy)))
    !freq = c/(lambda*dx)
    omega = 2*pai*freq
    lambda = c/freq
    jx = 0.0d0
    jy = 0.0d0
    jz = 0.0d0
    ajx = 0.0d0
    ajy = 0.0d0
    ajz = 0.0d0
    ajex = 0.0d0
    ajey = 0.0d0 
    ajez = 0.0d0

    if(kwave==0) then
        pc = orgs(1)*dx
        pw = pw*dx
    endif

    if(mode==1.or.mode==2.or.mode==3) then
        write(30,'(a11,e10.3)')"frequency:",freq 
        write(30,'(a3,e10.3)')"T:",1/freq
        write(30,'(a4,f10.3)')"ka:",2*pai*(radius*dx)/lambda
    endif

    write(30,'(a4,e10.3)')"dt:",dt
    write(30,'(a4,e10.3,a4,e10.3,/)')"dx:",dx,"dy:",dy
    
    !背景媒質
    do j=-1,ny
        do i=-1,nx
            epsd(i,j)=epsbk
            mud(i,j)=mubk
            sgmed(i,j)=sigebk
            sgmmd(i,j)=sigmbk
        end do
    end do

    !係数の計算
    do j=0,ny
        do i=0,nx
            epsx = 0.5d0*(epsd(i,j)+epsd(i,j-1))*eps0
            sgex = 0.5d0*(sgmed(i,j)+sgmed(i,j-1))
            a = 0.5d0*sgex*dt/epsx
            aex(i,j) = (1.0d0-a)/(1.0d0+a)
            bexy(i,j)= dt/epsx/(1.0d0+a)/dy

            epsy = 0.5d0*(epsd(i,j)+epsd(i-1,j))*eps0
            sgey = 0.5d0*(sgmed(i,j)+sgmed(i-1,j))
            a = 0.5d0*sgey*dt/epsy
            aey(i,j) = (1.0d0-a)/(1.0d0+a)
            beyx(i,j)= dt/epsy/(1.0d0+a)/dx
            
            epsz = 0.25d0*(epsd(i,j)+epsd(i-1,j)+epsd(i,j-1)+epsd(i-1,j-1))*eps0
            sgez = 0.25d0*(sgmed(i,j)+sgmed(i-1,j)+sgmed(i,j-1)+sgmed(i-1,j-1))
            a = 0.5d0*sgez*dt/epsz
            aez(i,j) = (1.0d0-a)/(1.0d0+a)
            bezx(i,j) = dt/epsz/(1.0d0+a)/dx 
            bezy(i,j) = dt/epsz/(1.0d0+a)/dy
            
            mux = 0.5d0*(mud(i,j)+mud(i-1,j))*mu0
            sgmx= 0.5d0*(sgmmd(i,j)+sgmmd(i-1,j))
            a = 0.5d0*sgmx*dt/mux
            amx(i,j) = (1.0d0-a)/(1.0d0+a)
            bmxy(i,j)= dt/mux/(1.0d0+a)/dy

            muy = 0.5d0*(mud(i,j)+mud(i,j-1))*mu0
            sgmy= 0.5d0*(sgmmd(i,j)+sgmmd(i,j-1))
            a = 0.5d0*sgmy*dt/muy
            amy(i,j) = (1.0d0-a)/(1.0d0+a)
            bmyx(i,j)= dt/muy/(1.0d0+a)/dx

            muz = mud(i,j)*mu0
            sgmz= sgmmd(i,j)
            a = 0.5d0*sgmz*dt/muz
            amz(i,j) = (1.0d0-a)/(1.0d0+a)
            bmzx(i,j) = dt/muz/(1.0d0+a)/dx
            bmzy(i,j) = dt/muz/(1.0d0+a)/dy
        enddo
    enddo

    !完全導体設置
    if(med==2) then
        write(30,'(a8)')"set PEC"
        call PEC()
    endif

    !プラズマの配置
    if(pls.ge.1) then
        write(30,'(a11)')"set plasma"
        if(pls.eq.1) then
            write(30,'(a12)')"uniform ADE"
            call ADE()
        else if(pls.eq.2) then
            write(30,'(a15)')"nonuniform ADE"
            call ADEg()
        else if(pls.eq.3) then
            write(30,'(a12)')"uniform JEC"
            call JEC()
        else if(pls.eq.4) then
            write(30,'(a15)')"nonuniform JEC"
            call JECg()
        else if(pls.eq.5) then
            write(30,'(a12)')"uniform EOM"
            call EOM()
        else if(pls.eq.6) then
            write(30,'(a15)')"nonuniform EOM"
            call EOMg()
        endif
    endif

    !誘電体設置
    if((med.eq.1).or.(mode.eq.6)) then
        write(30,'(a15)')"set Dielectric"
        call epsmu()
    endif

end subroutine setup