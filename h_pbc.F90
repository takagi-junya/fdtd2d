subroutine h_pbc()
    use constants
    implicit none
    integer :: i,j

    call exchg2d_periodicy0(hz(istart-1:iend+1,jstart-1:jend+1)) 
    call exchg2d_periodicy0(hx(istart-1:iend+1,jstart-1:jend+1))
    call exchg2d_periodicx0(hz(istart-1:iend+1,jstart-1:jend+1))
    call exchg2d_periodicx0(hy(istart-1:iend+1,jstart-1:jend+1))

    !for ex
    !do i=0,nx
    !    hz(i,0) = hz(i,ny-1)
    !enddo
    !for ez
    !do i=0,nx
    !    hx(i,0) = hx(i,ny-1)
    !enddo

    !for ey
    !do j=0,ny
    !    hz(0,j) = hz(nx-1,j)
    !enddo


    !do j=0,ny
    !    hy(0,j) = hy(nx-1,j)
    !enddo


end subroutine h_pbc