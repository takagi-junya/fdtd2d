module pwave
    implicit none
    real :: orge(3),orgb(3)
contains

    subroutine init_pwave
        use constants
        implicit none
        integer :: i,j
        real :: x
        if(myrank.eq.0) then 
            write(30,'(a10,/)')"set pwave"
        endif
        t=0.0
        do j=jstart,jend
            do i=istart,iend
                x = i*dx
                ey(i,j) = amps(2)*gs_ez(x)
            enddo
        enddo

        do j=jstart,jend
            do i=istart,iend
                x = i*dx
                ez(i,j) = amps(3)*gs_ez(x)
            enddo
        enddo

        t=t+0.5*dt
        do j=jstart,jend
            do i=istart,iend
                x = (i+0.5)*dx
                hy(i,j) = -amps(3)*gs_ez(x)/z0
            enddo
        enddo
        
        do j=jstart,jend
            do i=istart,iend
                x = (i+0.5)*dx
                hz(i,j) = amps(2)*gs_ez(x)/z0
            enddo
        enddo

    end subroutine


    real function gs_ez(x)
        use constants
        implicit none
        real,intent(in) :: x
        real :: xx,aa
        xx = x-c*t-pc
        aa = (xx/pw)*(xx/pw)
        gs_ez = exp(-aa)
    end function gs_ez

    real function gaussian(x,y,org)
        use constants
        implicit none
        real,intent(in) :: x,y,org(:)
        real :: displ

        displ = (x-org(1))*sin(ang(1))*cos(ang(2))+&
        &       (y-org(2))*sin(ang(1))*sin(ang(2))

        gaussian = exp(-displ*displ/pw/pw)
    end function gaussian

end module pwave

