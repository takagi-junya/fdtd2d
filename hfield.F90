subroutine hfield()
    use constants
    implicit none
    integer :: i,j

    !$omp parallel do
    do j=0,ny-1
        !$omp parallel do
        do i=1,nx-1
            hx(i,j) = amx(i,j)*hx(i,j)-bmxy(i,j)*(ez(i,j+1)-ez(i,j))
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    !$omp parallel do
    do j=1,ny-1
        !$omp parallel do
        do i=0,nx-1
            hy(i,j) = amy(i,j)*hy(i,j)+bmyx(i,j)*(ez(i+1,j)-ez(i,j))
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

    !$omp parallel do
    do j=0,ny-1
        !$omp parallel do
        do i=0,nx-1
            hz(i,j) = amz(i,j)*hz(i,j)&
            &       -bmzx(i,j)*(ey(i+1,j)-ey(i,j))&
            &       +bmzy(i,j)*(ex(i,j+1)-ex(i,j))
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do

end subroutine