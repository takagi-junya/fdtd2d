!-----------------------------------------------------------------------
!       計算結果の出力
!-----------------------------------------------------------------------
module output
   use HDF5
   use MPI 
   use hdfio 
   implicit none 
   integer :: i,j,k,l 
   integer :: is,ie,js,je 
   integer :: outputrank,outputrank0
   contains

   subroutine out_init()
      use HDF5
      use constants
      use hdfio
      implicit none
      if(myrank.eq.0) then 
         write(30,*)"io:",io 
         write(30,*)"jo:",jo 
      endif
      outputrank0 = 0
      if((io.ge.istart.and.io.le.iend).and.(jo.ge.jstart.and.jo.le.jend)) then 
         outputrank0 = myrank
         write(*,*)"outputrank:",outputrank0
         open(35,file="ex.txt")
         open(36,file="ey.txt")
         open(37,file="ez.txt")
      endif
      call mpi_reduce(outputrank0,outputrank,1,MPI_INTEGER4,MPI_SUM,0,comm2d,mpierr)
      call mpi_bcast(outputrank,1,MPI_INTEGER4,0,comm2d,mpierr)
      if(comp(1).eq.1) then
         call h5popen(filename(1),groupname(1),file_id(1),group_id(1),plist_id(1),0,comm2d,info,hdferr)
      endif
      if(comp(2).eq.1) then
         call h5popen(filename(2),groupname(2),file_id(2),group_id(2),plist_id(2),0,comm2d,info,hdferr)
      endif
      if(comp(3).eq.1) then
         call h5popen(filename(3),groupname(3),file_id(3),group_id(3),plist_id(3),0,comm2d,info,hdferr)
      endif
      if(comp(4).eq.1) then
         call h5popen(filename(4),groupname(4),file_id(4),group_id(4),plist_id(4),0,comm2d,info,hdferr)
      endif
      if(comp(5).eq.1) then
         call h5popen(filename(5),groupname(5),file_id(5),group_id(5),plist_id(5),0,comm2d,info,hdferr)
      endif
      if(comp(6).eq.1) then
         call h5popen(filename(6),groupname(6),file_id(6),group_id(6),plist_id(6),0,comm2d,info,hdferr)
      endif
      if(comp(7).eq.1) then
         call h5popen(filename(7),groupname(7),file_id(7),group_id(7),plist_id(7),0,comm2d,info,hdferr)
      endif
      if(comp(8).eq.1) then
         call h5popen(filename(8),groupname(8),file_id(8),group_id(8),plist_id(8),0,comm2d,info,hdferr)
      endif
      if(comp(9).eq.1) then
         call h5popen(filename(9),groupname(9),file_id(9),group_id(9),plist_id(9),0,comm2d,info,hdferr)
      endif
   end subroutine 

   subroutine out_emf(n)
      use constants 
      use hdf5
      use hdfio 
      use mpi 
      integer,intent(in) :: n
      if(myrank.eq.outputrank) then 
         write(35,*)ex(io,jo)
         write(36,*)ey(io,jo)
         write(37,*)ez(io,jo)
      endif
      if(mod(int(n-ostart),out).eq.0.and.n-ostart.ge.0.and.n.le.oend) then
         write(tag,'(I4.4)') h5count
         if(myrank.eq.0) then
            write(*,'(a8,I3.3)')"h5count:",h5count 
         endif
         if(comp(1).eq.1) then
            call wrt2p(group_id(1),plist_id(1),tag,dimsf2d,chunk_dims2d,coords,ex(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(2).eq.1) then
            call wrt2p(group_id(2),plist_id(2),tag,dimsf2d,chunk_dims2d,coords,ey(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(3).eq.1) then
            call wrt2p(group_id(3),plist_id(3),tag,dimsf2d,chunk_dims2d,coords,ez(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(4).eq.1) then
            call wrt2p(group_id(4),plist_id(4),tag,dimsf2d,chunk_dims2d,coords,hx(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(5).eq.1) then
            call wrt2p(group_id(5),plist_id(5),tag,dimsf2d,chunk_dims2d,coords,hy(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(6).eq.1) then
            call wrt2p(group_id(6),plist_id(6),tag,dimsf2d,chunk_dims2d,coords,hz(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(7).eq.1) then
            call wrt2p(group_id(7),plist_id(7),tag,dimsf2d,chunk_dims2d,coords,jx(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(8).eq.1) then
            call wrt2p(group_id(8),plist_id(8),tag,dimsf2d,chunk_dims2d,coords,jy(istart:iend:stride,jstart:jend:stride))
         endif
         if(comp(9).eq.1) then
            call wrt2p(group_id(9),plist_id(9),tag,dimsf2d,chunk_dims2d,coords,jz(istart:iend:stride,jstart:jend:stride))
         endif
         h5count = h5count + 1  
      endif
      
   end subroutine 
   subroutine out_finailzie()
      use constants
      use hdf5
      use hdfio 
      use mpi
      if(comp(1).eq.1) then
         call h5pclose(file_id(1),group_id(1),plist_id(1),istat1(1))
      endif
      if(comp(2).eq.1) then
         call h5pclose(file_id(2),group_id(2),plist_id(2),istat1(2))
      endif
      if(comp(3).eq.1) then
         call h5pclose(file_id(3),group_id(3),plist_id(3),istat1(3))
      endif
      if(comp(4).eq.1) then
         call h5pclose(file_id(4),group_id(4),plist_id(4),istat1(4))
      endif
      if(comp(5).eq.1) then
         call h5pclose(file_id(5),group_id(5),plist_id(5),istat1(5))
      endif
      if(comp(6).eq.1) then
         call h5pclose(file_id(6),group_id(6),plist_id(6),istat1(6))
      endif
      if(comp(7).eq.1) then
         call h5pclose(file_id(7),group_id(7),plist_id(7),istat1(7))
      endif
      if(comp(8).eq.1) then
         call h5pclose(file_id(8),group_id(8),plist_id(8),istat1(8))
      endif
      if(comp(9).eq.1) then
         call h5pclose(file_id(9),group_id(9),plist_id(8),istat1(9))
      endif
      if(myrank.eq.outputrank) then
         write(*,*)"close file"
         close(35)
         close(36)
         close(37)
      endif
   end subroutine
end module