subroutine exchg2di(data)
    use mpi
    use constants
    implicit none
    real(kind=8),intent(in) :: data(istart-1:iend+1,jstart-1:jend+1)
    call mpi_sendrecv(data(iend,jstart),1,edge,east,0,&
                    &data(istart-1,jstart),1,edge,west,0,&
                    &comm2d,mpi_status_ignore,mpierr)
    call mpi_sendrecv(data(istart,jstart),1,edge,west,1,&
                    &data(iend+1,jstart),1,edge,east,1,&
                    &comm2d,mpi_status_ignore,mpierr)
end subroutine

subroutine exchg2dj(data,lxx)
    use mpi
    use constants
    implicit none
    integer,intent(in) :: lxx
    real(kind=8),intent(in) :: data(istart-1:iend+1,jstart-1:jend+1)
    call mpi_sendrecv(data(istart,jend),lxx,MPI_DOUBLE_PRECISION,north,0,&
                    &data(istart,jstart-1),lxx,MPI_DOUBLE_PRECISION,south,0,&
                    &comm2d,mpi_status_ignore,mpierr)
    call mpi_sendrecv(data(istart,jstart),lxx,MPI_DOUBLE_PRECISION,south,1,&
                    &data(istart,jend+1),lxx,MPI_DOUBLE_PRECISION,north,1,&
                    &comm2d,mpi_status_ignore,mpierr)
end subroutine

subroutine exchg2d_periodicy0(fdata)
    use mpi 
    use constants 
    implicit none
    integer :: i
    real(kind=8),intent(inout) :: fdata(istart-1:iend+1,jstart-1:jend+1)

    if(coords(2).eq.0) then
        call mpi_recv(tmpyd(istart),lx0,MPI_DOUBLE_PRECISION,psouth,0,pcomm,mp_status(:,1),mpierr)
        do i=istart,iend 
            fdata(i,0) = tmpyd(i)
        enddo
    else if(coords(2).eq.ndims(2)-1) then
        do i=istart,iend 
            tmpyu(i) = fdata(i,ny-1)
        enddo
        call mpi_send(tmpyu(istart),lx0,MPI_DOUBLE_PRECISION,pnorth,0,pcomm,mpierr)
    endif
end subroutine 

subroutine exchg2d_periodicy1(fdata)
    use mpi 
    use constants 
    implicit none
    integer :: i
    real(kind=8),intent(inout) :: fdata(istart-1:iend+1,jstart-1:jend+1)
    if(coords(2).eq.0) then
        do i=istart,iend 
            tmpyd(i) = fdata(i,1)
        enddo
        call mpi_send(tmpyd(istart),lx0,MPI_DOUBLE_PRECISION,psouth,1,pcomm,mp_status(:,2),mpierr)
    else if(coords(2).eq.ndims(2)-1) then
        call mpi_recv(tmpyu(istart),lx0,MPI_DOUBLE_PRECISION,pnorth,1,pcomm,mpierr)
        do i=istart,iend 
            fdata(i,ny) = tmpyu(i)
        enddo
    endif
end subroutine 

subroutine exchg2d_periodicx0(fdata)
    use mpi 
    use constants 
    implicit none
    integer :: j
    real(kind=8),intent(inout) :: fdata(istart-1:iend+1,jstart-1:jend+1)
    if(coords(1).eq.0) then
        call mpi_recv(tmpyl(jstart),ly1,MPI_DOUBLE_PRECISION,pwest,2,pcomm,mp_status(:,3),mpierr)
        do j=jstart,jend 
            fdata(0,j) = tmpyl(j)
        enddo
    else if(coords(1).eq.ndims(1)-1) then
        do j=jstart,jend 
            tmpyr(j) = fdata(nx-1,j)
        enddo
        call mpi_send(tmpyr(jstart),ly1,MPI_DOUBLE_PRECISION,peast,2,pcomm,mpierr)
    endif
end subroutine 

subroutine exchg2d_periodicx1(fdata)
    use mpi 
    use constants 
    implicit none
    integer :: j
    real(kind=8),intent(inout) :: fdata(istart-1:iend+1,jstart-1:jend+1)
    if(coords(1).eq.0) then
        do j=jstart,jend 
            tmpyl(j) = fdata(1,j)
        enddo
        call mpi_send(tmpyl(jstart),ly1,MPI_DOUBLE_PRECISION,pwest,3,pcomm,mpierr)
    else if(coords(1).eq.ndims(1)-1) then
        call mpi_recv(tmpyr(jstart),ly1,MPI_DOUBLE_PRECISION,peast,3,pcomm,mp_status(:,4),mpierr)
        do j=jstart,jend 
            fdata(nx,j) = tmpyr(j)
        enddo
    endif
end subroutine 

subroutine exchg2d_lpml(data)
    use mpi 
    use constants 
    implicit none 
    real(kind=8),intent(in) :: data(0:lpml(1),jstart-1:jend+1)
    call mpi_sendrecv(data(0,jend),lpml(1),MPI_DOUBLE_PRECISION,north,0,&
                    &data(0,jstart-1),lpml(1),MPI_DOUBLE_PRECISION,south,0,&
                    &comm2d,mpi_status_ignore,mpierr)
    call mpi_sendrecv(data(0,jstart),lpml(1),MPI_DOUBLE_PRECISION,south,1,&
                    &data(0,jend+1),lpml(1),MPI_DOUBLE_PRECISION,north,1,&
                    &comm2d,mpi_status_ignore,mpierr)
end subroutine 

subroutine exchg2d_rpml(data)
    use mpi
    use constants
    implicit none
    real(kind=8),intent(in) :: data(nx-lpml(1):nx+1,jstart-1:jend+1)
    call mpi_sendrecv(data(nx-lpml(1),jend),lpml(1),MPI_DOUBLE_PRECISION,north,0,&
                    &data(nx-lpml(1),jstart-1),lpml(1),MPI_DOUBLE_PRECISION,south,0,&
                    &comm2d,mpi_status_ignore,mpierr)
    call mpi_sendrecv(data(nx-lpml(1),jstart),lpml(1),MPI_DOUBLE_PRECISION,south,1,&
                    &data(nx-lpml(1),jend+1),lpml(1),MPI_DOUBLE_PRECISION,north,1,&
                    &comm2d,mpi_status_ignore,mpierr)
end subroutine

subroutine exchg2d_upml(data)
    use mpi
    use constants
    real(kind=8),intent(in) :: data(istart-1:iend+1,1:lpml(2))
    call mpi_sendrecv(data(iend,1),1,pmlx0_MPI_DOUBLE_PRECISION,east,0,&
                    &data(istart-1,1),1,pmlx0_MPI_DOUBLE_PRECISION,west,0,&
                    &comm2d,mpi_status_ignore,mpierr)
    call mpi_sendrecv(data(istart,1),1,pmlx0_MPI_DOUBLE_PRECISION,west,1,&
                    &data(iend+1,1),1,pmlx0_MPI_DOUBLE_PRECISION,east,1,&
                    &comm2d,mpi_status_ignore,mpierr)
end subroutine

subroutine exchg2d_opml(data)
    use mpi
    use constants
    real(kind=8),intent(in) :: data(istart-1:iend+1,ny-lpml(2):ny)
    call mpi_sendrecv(data(iend,ny-lpml(2)),1,pmlx0_MPI_DOUBLE_PRECISION,east,0,&
                    &data(istart-1,ny-lpml(2)),1,pmlx0_MPI_DOUBLE_PRECISION,west,0,&
                    &comm2d,mpi_status_ignore,mpierr)
    call mpi_sendrecv(data(istart,ny-lpml(2)),1,pmlx0_MPI_DOUBLE_PRECISION,west,1,&
                    &data(iend+1,ny-lpml(2)),1,pmlx0_MPI_DOUBLE_PRECISION,east,1,&
                    &comm2d,mpi_status_ignore,mpierr)
end subroutine