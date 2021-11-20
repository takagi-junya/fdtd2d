subroutine current()
    use constants
    use omp_lib
    implicit none
    integer :: i,j

    !$omp parallel do
    do j=jstart,jend
        !$omp parallel do
        do i=istart,iend
            jx(i,j) = ajx(i,j)*jx(i,j)+ajex(i,j)*ex(i,j)
            jy(i,j) = ajy(i,j)*jy(i,j)+ajey(i,j)*ey(i,j)
            jz(i,j) = ajz(i,j)*jz(i,j)+ajez(i,j)*ez(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$om end parallel do
            
end subroutine 