subroutine hfield()
    use constants
    implicit none
    integer :: i,j

    call exchg2dj(ez(istart-1:iend+1,jstart-1:jend+1),lx0)
    
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
            hx(i,j) = amx(i,j)*hx(i,j)-bmxy(i,j)*(ez(i,j+1)-ez(i,j))
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

    call exchg2di(ez(istart-1:iend+1,jstart-1:jend+1))

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
            hy(i,j) = amy(i,j)*hy(i,j)+bmyx(i,j)*(ez(i+1,j)-ez(i,j))
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
    
    call exchg2di(ey(istart-1:iend+1,jstart-1:jend+1))
    call exchg2dj(ex(istart-1:iend+1,jstart-1:jend+1),lx0)

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
            hz(i,j) = amz(i,j)*hz(i,j)&
            &       -bmzx(i,j)*(ey(i+1,j)-ey(i,j))&
            &       +bmzy(i,j)*(ex(i,j+1)-ex(i,j))
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