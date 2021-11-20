!-----------------------------------------------------------------------
!   誘電体
!-----------------------------------------------------------------------
subroutine epsmu()   
    use constants 
    implicit none
    integer :: i,j
    real(kind=8) sp_rad,cy_rad
    
    if(med.eq.3) then
        if(myrank.eq.0) then
            write(30,'(a15)')"set Dielectric"
            write(30,*)"object:cylinder"
            write(30,*)"thickness of dielec:",erad-prad
            write(30,*)""
        endif
        do j=jstart,jend
            do i=istart,iend
                if(cy_rad(i,j).le.erad.and.cy_rad(i+1,j).le.erad) then
                    if(cy_rad(i,j).ge.prad.and.cy_rad(i+1,j).ge.prad) then
                        epsd(i,j) = epsr
                        mud(i,j) = 1.0d0
                        sgmed(i,j) = 0.0d0
                        sgmmd(i,j) = 0.0d0

                        epsx = epsd(i,j)*eps0
                        sgex = sgmed(i,j)
                        a = 0.5d0*sgex*dt/epsx
                        aex(i,j) = (1.0d0-a)/(1.0d0+a)
                        bexy(i,j)= dt/epsx/(1.0d0+a)/dy
                    endif
                endif
                if(cy_rad(i,j).le.erad.and.cy_rad(i,j+1).le.erad) then
                    if(cy_rad(i,j).ge.prad.and.cy_rad(i,j+1).ge.prad) then
                        epsd(i,j) = epsr
                        mud(i,j) = 1.0d0
                        sgmed(i,j) = 0.0d0
                        sgmmd(i,j) = 0.0d0

                        epsy = epsd(i,j)*eps0
                        sgey = sgmed(i,j)
                        a = 0.5d0*sgey*dt/epsy
                        aey(i,j) = (1.0d0-a)/(1.0d0+a)
                        beyx(i,j)= dt/epsy/(1.0d0+a)/dx
                    endif
                endif
                if(cy_rad(i,j).le.erad) then 
                    if(cy_rad(i,j).ge.prad) then
                        epsd(i,j) = epsr
                        mud(i,j) = 1.0d0
                        sgmed(i,j) = 0.0d0
                        sgmmd(i,j) = 0.0d0

                        epsz = epsd(i,j)*eps0
                        sgez = sgmed(i,j)
                        a = 0.5d0*sgez*dt/epsz
                        aez(i,j) = (1.0d0-a)/(1.0d0+a)
                        bezx(i,j) = dt/epsz/(1.0d0+a)/dx 
                        bezy(i,j) = dt/epsz/(1.0d0+a)/dy
                    endif
                endif
            end do
        end do
    else if(obj.eq.2) then
        write(30,*)"object:proper prism"
        do j=jc-ly2,jc+ly2-1
            do i=ic-lx2,ic+lx2-1
                epsd(i,j) = epsr
                mud(i,j) = 1.0d0
                sgmed(i,j) = 0.0d0
                sgmmd(i,j) = 0.0d0
            enddo
        enddo
    else if(obj.eq.3) then
        write(30,*)"object:dielec cylinder"
        write(30,*)""
        do j=jstart,jend
            do i=istart,iend
                if(cy_rad(i,j).le.erad.and.cy_rad(i+1,j).le.erad) then
                    if(cy_rad(i,j).ge.prad.and.cy_rad(i+1,j).ge.prad) then
                        epsd(i,j) = epsr
                        mud(i,j) = 1.0d0
                        sgmed(i,j) = 0.0d0
                        sgmmd(i,j) = 0.0d0
                    endif
                endif
                if(cy_rad(i,j).le.erad.and.cy_rad(i,j+1).le.erad) then
                    if(cy_rad(i,j).ge.prad.and.cy_rad(i,j+1).ge.prad) then
                        epsd(i,j) = epsr
                        mud(i,j) = 1.0d0
                        sgmed(i,j) = 0.0d0
                        sgmmd(i,j) = 0.0d0
                    endif
                endif
                if(cy_rad(i,j).le.erad) then
                    if(cy_rad(i,j).ge.prad) then
                        epsd(i,j) = epsr
                        mud(i,j) = 1.0d0
                        sgmed(i,j) = 0.0d0
                        sgmmd(i,j) = 0.0d0
                    endif
                endif
            end do
            end do
    endif
 end subroutine epsmu