subroutine e_pbc()
    use constants
    implicit none
    integer :: i,j

    call exchg2d_periodicy1(ez(istart-1:iend+1,jstart-1:jend+1))
    call exchg2d_periodicy1(ex(istart-1:iend+1,jstart-1:jend+1))
    call exchg2d_periodicx1(ez(istart-1:iend+1,jstart-1:jend+1))
    call exchg2d_periodicx1(ey(istart-1:iend+1,jstart-1:jend+1))

    !for hx
    !do i=1,nx-1
    !    ez(i,ny) = ez(i,1)
    !enddo

    !for hz
    !do i=0,nx-1
    !    ex(i,ny) = ex(i,1)
    !enddo

    !for hy
    !do j=1,ny-1
    !    ez(nx,j) = ez(1,j)
    !enddo


    !do j=0,ny-1
    !    ey(nx,j) = ey(1,j)
    !enddo
end subroutine e_pbc