subroutine currentEOM()
    use constants
    use omp_lib
    implicit none
    integer :: i,j

    !$omp parallel do
    do j=0,ny
        !$omp parallel do
        do i=0,nx
            jx(i,j) = -qe*nd(i,j)*vx(i,j)
            jy(i,j) = -qe*nd(i,j)*vy(i,j)
            jz(i,j) = -qe*nd(i,j)*vz(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$om end parallel do
            
end subroutine 