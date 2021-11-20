subroutine EOM()
    use constants
    use HDF5
    use hdfio
    implicit none
    integer :: i,j
    real(kind=8) cy_rad

    if(myrank.eq.0) then
        write(30,'(a12)')"uniform EOM"
        write(30,'(a8,e10.3)')"omegap:",wp 
        write(30,'(a16,e10.3)')"collision freq:",nu
        write(30,'(a21,f10.2,/)')"thickness of plasma:",prad-radius
    endif
    
    ajj = dt/eps0

    omat = dt/2.0d0*omat
    sa = adding_mat(imat(:,:),omat(:,:))
    sa = inversing_mat(sa(:,:))

    sb = adding_mat(imat(:,:),-omat(:,:))

    sab = matmul(sa,sb)
    tc = qe*dt/mel*sa
    do j=jstart,jend
        do i=istart,iend
            if(cy_rad(i,j).le.prad.and.cy_rad(i+1,j).le.prad) then
                if(cy_rad(i,j).ge.radius.and.cy_rad(i+1,j).ge.radius) then
                    avx(i,j)  = (2.0d0-nu*dt)/(2.0d0+nu*dt)
                    ajex(i,j) = -(2.0d0*qe/mel*dt)/(2.0d0+nu*dt)
                endif
            endif
            if(cy_rad(i,j).le.prad.and.cy_rad(i,j+1).le.prad) then
                if(cy_rad(i,j).ge.radius.and.cy_rad(i,j+1).ge.radius) then
                    avy(i,j)  = (2.0d0-nu*dt)/(2.0d0+nu*dt)
                    ajey(i,j) = -(2.0d0*qe/mel*dt)/(2.0d0+nu*dt)
                endif
            endif
            if(cy_rad(i,j).le.prad) then
                if(cy_rad(i,j).ge.radius) then
                    avz(i,j)  = (2.0d0-nu*dt)/(2.0d0+nu*dt)
                    ajez(i,j) = -(2.0d0*qe/mel*dt)/(2.0d0+nu*dt)
                endif
            endif
        enddo
    enddo

    do j=jstart,jend
        do i=istart,iend
            if(cy_rad(i,j).le.prad) then
                if(cy_rad(i,j).ge.radius) then
                    nd(i,j) = mel*eps0*(wp**2.0d0)/(qe**2.0d0)
                else
                    nd(i,j) = 0.0d0
                endif
            else 
                nd(i,j) = 0.0d0
            endif
        enddo
    enddo

    contains 

    function adding_mat(x,y) result(ans)
        use constants
        implicit none 
        integer :: i,j
        real(kind=8),intent(in) :: x(3,3),y(3,3)
        real(kind=8) :: ans(3,3)
        
        do i=1,3
            do j=1,3
                ans(i,j) = x(i,j) + y(i,j)
            enddo
        enddo
        
    end function 

    subroutine show_mat(x)
        use constants 
        implicit none 
        integer :: i,j 
        real(kind=8),intent(in) :: x(3,3)
        real(kind=8) :: y(3,3)
        do i = 1,3
            write(30,*),x(i,1),x(i,2),x(i,3)
        enddo
        write(30,*)""
    end subroutine


    function inversing_mat(x) result(ans)
        use constants
        implicit none
        integer :: n,lwork,ifail,infom
        integer :: ipiv(3),work(192)
        integer :: i,j
        real(kind=8) :: ans(3,3),tmp(3,3)
        real(kind=8),intent(inout) :: x(3,3)
        n = 3
        lwork = n*64
        tmp = x
        call dgetrf(n,n,x,n,ipiv,infom)
        call dgetri(n,x,n,ipiv,work,lwork,infom)
        tmp = matmul(tmp,x)
        ans = x
    end function
end subroutine EOM