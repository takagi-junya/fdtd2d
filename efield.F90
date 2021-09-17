subroutine efield()
    use constants
    implicit none
    integer :: i,j

    !$omp parallel do
    do j=1,ny-1
        !$omp parallel do
        do i=0,nx-1
            ex(i,j) = aex(i,j)*ex(i,j) + bexy(i,j)*(hz(i,j)-hz(i,j-1)) - ajj*jx(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    !$omp parallel do
    do j=0,ny-1
        !$omp parallel do
        do i=1,nx-1
            ey(i,j) = aey(i,j)*ey(i,j) - beyx(i,j)*(hz(i,j)-hz(i-1,j)) - ajj*jy(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    !$omp parallel do
    do j=1,ny-1
        !$omp parallel do
        do i=1,nx-1
            ez(i,j) = aez(i,j)*ez(i,j)&
            &       + bezx(i,j)*(hy(i,j)-hy(i-1,j))&
            &       - bezy(i,j)*(hx(i,j)-hx(i,j-1))&
            &       - ajj*jz(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do
     
end subroutine