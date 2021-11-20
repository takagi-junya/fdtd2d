subroutine efield()
    use constants
    implicit none
    integer :: i,j

    call exchg2dj(hz(istart-1:iend+1,jstart-1:jend+1),lx0)

    if(jend.eq.ny) then
        jend = ny-1
    endif
    if(iend.eq.nx) then
        iend = nx-1
    endif

    !$omp parallel do
    do j=jstart,jend
        !$omp parallel do
        do i=istart,iend
            ex(i,j) = aex(i,j)*ex(i,j) + bexy(i,j)*(hz(i,j)-hz(i,j-1)) - ajj*jx(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    if(jend.eq.ny-1) then
        jend = ny
    endif
    if(iend.eq.nx-1) then
        iend = nx
    endif

    call exchg2di(hz(istart-1:iend+1,jstart-1:jend+1))
    
    if(jend.eq.ny) then
        jend = ny-1
    endif
    if(iend.eq.nx) then
        iend = nx-1
    endif

    !$omp parallel do
    do j=jstart,jend
        !$omp parallel do
        do i=istart,iend
            ey(i,j) = aey(i,j)*ey(i,j) - beyx(i,j)*(hz(i,j)-hz(i-1,j)) - ajj*jy(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    if(jend.eq.ny-1) then
        jend = ny
    endif
    if(iend.eq.nx-1) then
        iend = nx
    endif
    
    call exchg2di(hy(istart-1:iend+1,jstart-1:jend+1))
    call exchg2dj(hx(istart-1:iend+1,jstart-1:jend+1),lx0)

    if(jend.eq.ny) then
        jend = ny-1
    endif
    if(iend.eq.nx) then
        iend = nx-1
    endif

    !$omp parallel do
    do j=jstart,jend
        !$omp parallel do
        do i=istart,iend
            ez(i,j) = aez(i,j)*ez(i,j)&
            &       + bezx(i,j)*(hy(i,j)-hy(i-1,j))&
            &       - bezy(i,j)*(hx(i,j)-hx(i,j-1))&
            &       - ajj*jz(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    if(jend.eq.ny-1) then
        jend = ny
    endif
    if(iend.eq.nx-1) then
        iend = nx
    endif
end subroutine