subroutine velocity
    use constants
    implicit none
    integer :: i,j

    call exchg2di(nd1(istart-1:iend+1,jstart-1:jend+1))
    call exchg2dj(nd1(istart-1:iend+1,jstart-1:jend+1),lx0)
    !$omp parallel do
    do j=jstart,jend
        !$omp parallel do
        do i = istart,iend
            vx(i,j) = sab(1,1)*vx(i,j)+sab(1,2)*vy(i,j)+sab(1,3)*vz(i,j)-tc(1,1)*ex(i,j)-tc(1,2)*ey(i,j)-tc(1,3)*ez(i,j)-tu(1,1)*(nd1(i+1,j)-nd1(i,j))/dx-tu(1,2)*(nd1(i,j+1)-nd1(i,j))/dy
            vy(i,j) = sab(2,1)*vx(i,j)+sab(2,2)*vy(i,j)+sab(2,3)*vz(i,j)-tc(2,1)*ex(i,j)-tc(2,2)*ey(i,j)-tc(2,3)*ez(i,j)-tu(2,1)*(nd1(i+1,j)-nd1(i,j))/dx-tu(2,2)*(nd1(i,j+1)-nd1(i,j))/dy
            vz(i,j) = sab(3,1)*vx(i,j)+sab(3,2)*vy(i,j)+sab(3,3)*vz(i,j)-tc(3,1)*ex(i,j)-tc(3,2)*ey(i,j)-tc(3,3)*ez(i,j)-tu(3,1)*(nd1(i+1,j)-nd1(i,j))/dx-tu(3,2)*(nd1(i,j+1)-nd1(i,j))/dy
        enddo
        !$omp end parallel do
    enddo
    !$omp end parallel do
end subroutine 
            
