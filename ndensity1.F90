subroutine ndensity1
    use constants
    call exchg2dj(vx(istart-1:iend+1,jstart-1:jend+1),lx0)
    call exchg2dj(vy(istart-1:iend+1,jstart-1:jend+1),lx0)
    call exchg2di(vx(istart-1:iend+1,jstart-1:jend+1))
    call exchg2di(vy(istart-1:iend+1,jstart-1:jend+1))
    
    !$omp parallel do 
    do i=istart,iend
        !$omp parallel do
        do j=jstart,jend
            nd1(i,j) = nd1(i,j)-andx(i,j)*(vx(i,j)-vx(i-1,j))-andy(i,j)*(vy(i,j)-vy(i,j-1))
            nd(i,j) = nd(i,j)+nd1(i,j)
        enddo 
        !$omp end parallel do
    enddo
    !$omp end parallel do 

end subroutine