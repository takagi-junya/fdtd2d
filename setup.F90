subroutine setup()
    use HDF5
    use MPI
    use constants
    use hdfio
    implicit none
    integer :: i,j,k

    stride = 1
    ostart=0
    erad = 0.0
    odom = (/0,0,0,0/)

    namelist /space/nxx,nyy,dx,dy,pbc,abc,lpml
    namelist /time/deltat,nstep
    namelist /output/out,ostart,oend,odom,stride,comp,io,jo
    namelist /scatt/mode,wshape,lx,ly,gamma0,theta0,phi0,amp,freq,tau0!lambda,tau0
    namelist /far/isx,isy,theta1,phi1
    namelist /object/obj,med,ic,jc,lx2,ly2,epsr,radius
    namelist /wave/kwave,amps,orgs,ang,pt,pw
    namelist /plasma/pls,prad,erad,nu,wp,wc
    namelist /mpip/prx,pry

    open(10,file="param.inp",action="read")
    read(10,nml=space)
    read(10,nml=time)
    read(10,nml=output)
    read(10,nml=scatt)
    read(10,nml=far)
    read(10,nml=object)
    read(10,nml=wave)    
    read(10,nml=plasma)
    read(10,nml=mpip)
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
    ang  = dble(ang)
    pw = dble(pw)
    pt = dble(pt)
    prad = dble(prad)
    erad = dble(erad)
    nu = dble(nu)
    wp = dble(wp)
    wc(1) = dble(wc(1))
    wc(2) = dble(wc(2))
    wc(3) = dble(wc(3))
    
    !吸収境界
    if(abc.eq.0) then
        nx = nxx
        ny = nyy
    else
        nx = nxx+2*lpml(1)
        ny = nyy+2*lpml(2)
    endif

    ndims(1) = prx
    ndims(2) = pry

    call mpi_cart_create(comm,ndim,ndims,isperiodic,reorder,comm2d,mpierr)
    call mpi_cart_shift(comm2d,0,1,west,east,mpierr)
    call mpi_cart_shift(comm2d,1,1,south,north,mpierr)
    call mpi_cart_coords(comm2d,myrank,ndim,coords,mpierr)

    istart = int(nx/ndims(1))*coords(1)+1
    iend   = int(nx/ndims(1))*(coords(1)+1)
    jstart = int(ny/ndims(2))*coords(2)+1
    jend   = int(ny/ndims(2))*(coords(2)+1)

    lx1 = (iend+1)-(istart-1)+1
    ly1 = jend-jstart+1
    lx0 = iend-istart+1

    if(pbc.eq.1) then 
        call mpi_cart_create(comm,ndim,ndims,pisperiodic,reorder,pcomm,mpierr)
        call mpi_cart_shift(pcomm,0,1,pwest,peast,mpierr)
        call mpi_cart_shift(pcomm,1,1,psouth,pnorth,mpierr)
        allocate(tmpyd(istart:iend),tmpyu(istart:iend))
        allocate(tmpyl(jstart:jend),tmpyr(jstart:jend))
    endif

    call mpi_type_vector(ly1,1,lx1,MPI_DOUBLE_PRECISION,edge,mpierr)
    call mpi_type_commit(edge,mpierr)

    call mpi_type_vector(lpml(2),1,lx1,MPI_DOUBLE_PRECISION,pml_edge,mpierr)
    call mpi_type_commit(pml_edge,mpierr)

    dimsf2d(1) = int(nx/stride)
    dimsf2d(2) = int(ny/stride)
    !dimsf2d = (/nx,ny/)
    
    chunk_dims2d(1) = int((iend-istart+1)/stride)
    chunk_dims2d(2) = int((jend-jstart+1)/stride)
    !chunk_dims2d(1) = iend-istart+1
    !chunk_dims2d(2) = jend-jstart+1

    ic = int(0.5*nx)
    jc = int(0.5*ny)
    jo = jc


    allocate(ex(istart-1:iend+1,jstart-1:jend+1),ey(istart-1:iend+1,jstart-1:jend+1),ez(istart-1:iend+1,jstart-1:jend+1))
    allocate(hx(istart-1:iend+1,jstart-1:jend+1),hy(istart-1:iend+1,jstart-1:jend+1),hz(istart-1:iend+1,jstart-1:jend+1))
    
    allocate(aex(istart-1:iend+1,jstart-1:jend+1),aey(istart-1:iend+1,jstart-1:jend+1),aez(istart-1:iend+1,jstart-1:jend+1))
    allocate(bexy(istart-1:iend+1,jstart-1:jend+1),beyx(istart-1:iend+1,jstart-1:jend+1))
    allocate(bezx(istart-1:iend+1,jstart-1:jend+1),bezy(istart-1:iend+1,jstart-1:jend+1))
    
    allocate(amx(istart-1:iend+1,jstart-1:jend+1),amy(istart-1:iend+1,jstart-1:jend+1),amz(istart-1:iend+1,jstart-1:jend+1))
    allocate(bmxy(istart-1:iend+1,jstart-1:jend+1),bmyx(istart-1:iend+1,jstart-1:jend+1))
    allocate(bmzx(istart-1:iend+1,jstart-1:jend+1),bmzy(istart-1:iend+1,jstart-1:jend+1))
    
    allocate(jx(istart-1:iend+1,jstart-1:jend+1),jy(istart-1:iend+1,jstart-1:jend+1),jz(istart-1:iend+1,jstart-1:jend+1))
    
    allocate(epsd(istart-1:iend+1,jstart-1:jend+1),mud(istart-1:iend+1,jstart-1:jend+1))
    allocate(sgmed(istart-1:iend+1,jstart-1:jend+1),sgmmd(istart-1:iend+1,jstart-1:jend+1))

    if(pls.ge.1) then
        allocate(nd(istart-1:iend+1,jstart-1:jend+1))
        allocate(ajx(istart-1:iend+1,jstart-1:jend+1),ajy(istart-1:iend+1,jstart-1:jend+1),ajz(istart-1:iend+1,jstart-1:jend+1))
        allocate(ajex(istart-1:iend+1,jstart-1:jend+1),ajey(istart-1:iend+1,jstart-1:jend+1),ajez(istart-1:iend+1,jstart-1:jend+1))
        if(pls.ge.5) then
            allocate(vx(istart-1:iend+1,jstart-1:jend+1),vy(istart-1:iend+1,jstart-1:jend+1),vz(istart-1:iend+1,jstart-1:jend+1))
            allocate(avx(istart-1:iend+1,jstart-1:jend+1),avy(istart-1:iend+1,jstart-1:jend+1),avz(istart-1:iend+1,jstart-1:jend+1))
        endif
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
    omega = 2*pai*freq 
    lambda = c/freq 

    pc = orgs(1)*dx
    pw = pw*dx

    if(myrank.eq.0) then
        open(30,file="sim.out")
        write(30,'(a4,i4,a4,i4)')"nx:",nx,"ny:",ny
        write(30,'(a4,i4,a4,i4,/)')"ic:",ic,"jc:",jc
        write(30,'(a4,e10.3)')"dt:",dt
        write(30,'(a7,e10.3,a7,e10.3)')"wp^-1:",1/wp,"wc^-1:",1/wc(3) 
        write(30,'(a4,e10.3,a4,e10.3,/)')"dx:",dx,"dy:",dy
        write(30,'(a11,e10.3)')"frequency:",freq 
        write(30,'(a3,e10.3)')"T:",1/freq
        write(30,'(a8,f10.2)')"lambda:",lambda
        write(30,'(a4,f10.3)')"ka:",2*pai*(radius*dx)/lambda
    endif
    
    !背景媒質
    do j=jstart,jend
        do i=istart,iend
            epsd(i,j)=epsbk
            mud(i,j)=mubk
            sgmed(i,j)=sigebk
            sgmmd(i,j)=sigmbk
        end do
    end do
    
    !係数の計算
    call exchg2di(epsd(istart-1:iend+1,jstart-1:jend+1))
    call exchg2di(mud(istart-1:iend+1,jstart-1:jend+1))
    call exchg2di(sgmed(istart-1:iend+1,jstart-1:jend+1))
    call exchg2di(sgmmd(istart-1:iend+1,jstart-1:jend+1))

    call exchg2dj(epsd(istart-1:iend+1,jstart-1:jend+1),lx0)
    call exchg2dj(mud(istart-1:iend+1,jstart-1:jend+1),lx0)
    call exchg2dj(sgmed(istart-1:iend+1,jstart-1:jend+1),lx0)
    call exchg2dj(sgmmd(istart-1:iend+1,jstart-1:jend+1),lx0)

    do j=jstart,jend
        do i=istart,iend
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
    if(med.ge.2) then
        if(myrank.eq.0) then
            write(30,'(a8)')"set PEC"
        endif
        call PEC()
    endif
    
    imat(1,1) = 1
    imat(2,2) = 1
    imat(3,3) = 1

    
    omat(1,1) = nu
    omat(1,2) = wc(3)
    omat(1,3) =-wc(2)

    omat(2,1) =-wc(3)
    omat(2,2) = nu 
    omat(2,3) = wc(1)

    omat(3,1) = wc(2)
    omat(3,2) =-wc(1)
    omat(3,3) = nu 
    
    !プラズマの配置
    if(pls.ge.1) then
        if(myrank.eq.0) then
            write(30,'(a11)')"set plasma"
        endif
        if(pls.eq.1) then
            call ADE()
        else if(pls.eq.2) then
            call ADEg()
        else if(pls.eq.3) then
            call JEC()
        else if(pls.eq.4) then
            call JECg()
        else if(pls.eq.5) then
            call EOM()
        else if(pls.eq.6) then
            call EOMg()
        endif
    endif

    !誘電体設置
    if((med.ge.1)) then
        call epsmu()
    endif
end subroutine setup