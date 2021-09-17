subroutine JECg()
    use HDF5
    use hdfio
    use constants
    implicit none
    integer :: i,j
    real(kind=8) :: wpp,outw,rdr
    integer :: is,ie,js,je
    real(kind=8) cy_rad

    is = int(nx*0.5d0)-int(prad*1.1)
    ie = int(nx*0.5d0)+int(prad*1.1)
    js = int(ny*0.5d0)-int(prad*1.1)
    je = int(ny*0.5d0)+int(prad*1.1)
    dims(1) = ie-is+1
    dims(2) = je-js+1
    write(30,'(a8,e10.3)')"omegap:",wp 
    write(30,'(a16,e10.3)')"collision freq:",nu
    write(30,'(a21,f10.2,/)')"thickness of plasma:",prad-radius
    ajj = dt/eps0
    do j=0,ny
        do i=0,nx
            if(cy_rad(i,j).le.prad.and.cy_rad(i+1,j).le.prad) then
                if(cy_rad(i,j).ge.radius.and.cy_rad(i+1,j).ge.radius) then
                    rdr = ((cy_rad(i+1,j)+cy_rad(i,j))/2-radius)/(prad-radius)
                    wpp = wp*sqrt(bessel_j0(2.40d0*rdr))
                    if(isnan(wpp)) then
                        wpp = 0.0d0
                    endif
                    ajx(i,j) = exp(-nu*dt)
                    ajex(i,j) = eps0*(wpp**2.0d0)*exp(-nu*dt*0.50d0)*dt
                endif
            endif
            if(cy_rad(i,j).le.prad.and.cy_rad(i,j+1).le.prad) then
                if(cy_rad(i,j).ge.radius.and.cy_rad(i,j+1).ge.radius) then
                    rdr = ((cy_rad(i+1,j)+cy_rad(i,j))/2-radius)/(prad-radius)
                    wpp = wp*sqrt(bessel_j0(2.40d0*rdr))
                    if(isnan(wpp)) then
                        wpp = 0.0d0
                    endif
                    ajy(i,j) = exp(-nu*dt)
                    ajey(i,j) = eps0*(wpp**2.0d0)*exp(-nu*dt*0.50d0)*dt
                endif
            endif
            if(cy_rad(i,j).le.prad) then
                if(cy_rad(i,j).ge.radius) then
                    rdr = ((cy_rad(i+1,j)+cy_rad(i,j))/2-radius)/(prad-radius)
                    wpp = wp*sqrt(bessel_j0(2.40d0*rdr))
                    if(isnan(wpp)) then
                        wpp = 0.0d0
                    endif
                    ajz(i,j) = exp(-nu*dt)
                    ajez(i,j) = eps0*(wpp**2.0d0)*exp(-nu*dt*0.50d0)*dt
                endif
            endif
        enddo
    enddo
    
    open(31,file="wp.txt")
    do j=0,ny
        do i=0,nx
            if(cy_rad(i,j).le.prad) then
                if(cy_rad(i,j).ge.radius) then
                    rdr = ((cy_rad(i+1,j)+cy_rad(i,j))/2-radius)/(prad-radius)
                    wpp = wp*sqrt(bessel_j0(2.40d0*rdr))
                    if(isnan(wpp)) then
                        wpp = 0.0d0
                    endif
                    nd(i,j) = mel*eps0*(wpp**2.0d0)/(qe**2.0d0)
                    if(j.eq.int(ny/2)) then
                        write(31,*) wpp
                    endif
                else
                    nd(i,j) = 0.0d0
                    if(j.eq.int(ny/2)) then
                        write(31,*) 0.0d0
                    endif
                endif
            else
                nd(i,j) = 0.0d0
                if(j.eq.int(ny/2)) then
                    write(31,*) 0.0d0
                endif
            endif
        enddo
    enddo
    close(31)

    write(30,'(a15)')"output density"
    call hdfopen(filename(10),groupname(10),file_id(10),group_id(10),0)
    call wrt2d(file_id(10),group_id(10),"0000",dims,nd(is:ie,js:je),istat1(10),istat2(10))
    call hdfclose(file_id(10),group_id(10),error(10))
end subroutine