subroutine h_pbc()
    use constants
    implicit none
    integer :: i,j

    !for ex
    do i=0,nx
        hz(i,0) = hz(i,ny-1)
    enddo

    !for ey
    do j=0,ny
        hz(0,j) = hz(nx-1,j)
    enddo

    !for ez
    do i=0,nx
        hx(i,0) = hx(i,ny-1)
    enddo

    do j=0,ny
        hy(0,j) = hy(nx-1,j)
    enddo

end subroutine h_pbc