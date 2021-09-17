subroutine velocity
    use constants
    implicit none
    integer :: i,j

    !$omp parallel do
    do j=0,ny
        !$omp parallel do
        do i = 0,nx
            vx(i,j) = avx(i,j)*vx(i,j)+ajex(i,j)*ex(i,j)
            vy(i,j) = avy(i,j)*vy(i,j)+ajey(i,j)*ey(i,j)
            vz(i,j) = avz(i,j)*vz(i,j)+ajez(i,j)*ez(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do
end subroutine 
            
