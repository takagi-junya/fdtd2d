subroutine e_pbc()
    use constants
    implicit none
    integer :: i,j
    !for hx
    do i=1,nx-1
        ez(i,ny) = ez(i,1)
    enddo

    !for hy
    do j=1,ny-1
        ez(nx,j) = ez(1,j)
    enddo

    !for hz
    do i=0,nx-1
        ex(i,ny) = ex(i,1)
    enddo

    do j=0,ny-1
        ey(nx,j) = ey(1,j)
    enddo
end subroutine e_pbc