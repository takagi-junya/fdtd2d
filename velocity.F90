subroutine velocity
    use constants
    implicit none
    integer :: i,j

    !$omp parallel do
    do j=jstart,jend
        !$omp parallel do
        do i = istart,iend
            vx(i,j) = sab(1,1)*vx(i,j)+sab(1,2)*vy(i,j)+sab(1,3)*vz(i,j)-tc(1,1)*ex(i,j)-tc(1,2)*ey(i,j)-tc(1,3)*ez(i,j)
            vy(i,j) = sab(2,1)*vx(i,j)+sab(2,2)*vy(i,j)+sab(2,3)*vz(i,j)-tc(2,1)*ex(i,j)-tc(2,2)*ey(i,j)-tc(2,3)*ez(i,j)
            vz(i,j) = sab(3,1)*vx(i,j)+sab(3,2)*vy(i,j)+sab(3,3)*vz(i,j)-tc(3,1)*ex(i,j)-tc(3,2)*ey(i,j)-tc(3,3)*ez(i,j)
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do
end subroutine 
            
