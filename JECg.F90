subroutine JECg()
    use HDF5
    use hdfio
    use constants
    implicit none
    integer :: i,j
    real(kind=8) :: wpp,outw,rdr
    real(kind=8) cy_rad

    if(myrank.eq.0) then
        write(30,'(a15)')"nonuniform JEC"
        write(30,'(a8,e10.3)')"omegap:",wp 
        write(30,'(a16,e10.3)')"collision freq:",nu
        write(30,'(a21,f10.2,/)')"thickness of plasma:",prad-radius
    endif
    ajj = dt/eps0
    do j=jstart,jend
        do i=istart,iend
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
    
    do j=jstart,jend
        do i=istart,iend
            if(cy_rad(i,j).le.prad) then
                if(cy_rad(i,j).ge.radius) then
                    rdr = ((cy_rad(i+1,j)+cy_rad(i,j))/2-radius)/(prad-radius)
                    wpp = wp*sqrt(bessel_j0(2.40d0*rdr))
                    if(isnan(wpp)) then
                        wpp = 0.0d0
                    endif
                    nd(i,j) = mel*eps0*(wpp**2.0d0)/(qe**2.0d0)
                else
                    nd(i,j) = 0.0d0
                endif
            else
                nd(i,j) = 0.0d0
            endif
        enddo
    enddo

end subroutine