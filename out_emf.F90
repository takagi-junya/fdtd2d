
!-----------------------------------------------------------------------
!       計算結果の出力
!-----------------------------------------------------------------------
subroutine out_emf(n)
   use HDF5
   use constants
   use hdfio
   implicit none
  
   integer :: i,j,k,l
   integer,intent(in) :: n
   integer :: is,ie,js,je

   is=lpml(1)+odom(1)
   ie=lpml(1)+odom(2)
   js=lpml(2)+odom(3)
   je=lpml(2)+odom(4)
   dims(1)=int((ie-is)/stride)+1
   dims(2)=int((je-js)/stride)+1

   if(n.eq.0) then
     write(30,*)"io:",io," jo:",jo
   endif
   open(100,file="ex.txt")
   write(100,*)ex(io,jo)
   if(obj.eq.0) then
      open(101,file="eyin.txt")
      write(101,*)ey(io,jo)
   else
      open(103,file="ey.txt")
      write(103,*)ey(io,jo)
   endif
   open(102,file="ez.txt")
   write(102,*)ez(io,jo)
   
   if((n-ostart).eq.0) then
      if(comp(1).eq.1) then
         call hdfopen(filename(1),groupname(1),file_id(1),group_id(1),0)
      endif
      if(comp(2).eq.1) then
         call hdfopen(filename(2),groupname(2),file_id(2),group_id(2),0)
      endif
      if(comp(3).eq.1) then
         call hdfopen(filename(3),groupname(3),file_id(3),group_id(3),0)
      endif
      if(comp(4).eq.1) then
         call hdfopen(filename(4),groupname(4),file_id(4),group_id(4),0)
      endif
      if(comp(5).eq.1) then
         call hdfopen(filename(5),groupname(5),file_id(5),group_id(5),0)
      endif
      if(comp(6).eq.1) then
         call hdfopen(filename(6),groupname(6),file_id(6),group_id(6),0)
      endif
      if(comp(7).eq.1) then
         call hdfopen(filename(7),groupname(7),file_id(7),group_id(7),0)
      endif
      if(comp(8).eq.1) then
         call hdfopen(filename(8),groupname(8),file_id(8),group_id(8),0)
      endif
      if(comp(9).eq.1) then
         call hdfopen(filename(9),groupname(9),file_id(9),group_id(9),0)
      endif
   elseif(mod(int(n-ostart),out).eq.0) then
      write(tag,'(I4.4)') h5count
      write(*,'(/,a7,a4,/)')'h5 out:',tag
      if(comp(1).eq.1) then
         call wrt2d(file_id(1),group_id(1),tag,dims,ex(is:ie:stride,js:je:stride),istat1(1),istat2(1))
      endif
      if(comp(2).eq.1) then
         call wrt2d(file_id(2),group_id(2),tag,dims,ey(is:ie:stride,js:je:stride),istat1(2),istat2(2))
      endif
      if(comp(3).eq.1) then
         call wrt2d(file_id(3),group_id(3),tag,dims,ez(is:ie:stride,js:je:stride),istat1(3),istat2(3))
      endif
      if(comp(4).eq.1) then
         call wrt2d(file_id(4),group_id(4),tag,dims,hx(is:ie:stride,js:je:stride),istat1(4),istat2(4))
      endif
      if(comp(5).eq.1) then
         call wrt2d(file_id(5),group_id(5),tag,dims,hy(is:ie:stride,js:je:stride),istat1(5),istat2(5))
      endif
      if(comp(6).eq.1) then
         call wrt2d(file_id(6),group_id(6),tag,dims,hz(is:ie:stride,js:je:stride),istat1(6),istat2(6))
      endif
      if(comp(7).eq.1) then
         call wrt2d(file_id(7),group_id(7),tag,dims,jx(is:ie:stride,js:je:stride),istat1(7),istat2(7))
      endif
      if(comp(8).eq.1) then
         call wrt2d(file_id(8),group_id(8),tag,dims,jy(is:ie:stride,js:je:stride),istat1(8),istat2(8))
      endif
      if(comp(9).eq.1) then
         call wrt2d(file_id(9),group_id(9),tag,dims,jz(is:ie:stride,js:je:stride),istat1(9),istat2(9))
      endif
      
      h5count = h5count + 1  
   elseif(n.eq.nstep) then
      if(comp(1).eq.1) then
         call hdfclose(file_id(1),group_id(1),error(1))
      endif
      if(comp(2).eq.1) then
         call hdfclose(file_id(2),group_id(2),error(2))
      endif
      if(comp(3).eq.1) then
         call hdfclose(file_id(3),group_id(3),error(3))
      endif
      if(comp(4).eq.1) then
         call hdfclose(file_id(4),group_id(4),error(4))
      endif
      if(comp(5).eq.1) then
         call hdfclose(file_id(5),group_id(5),error(5))
      endif
      if(comp(6).eq.1) then
         call hdfclose(file_id(6),group_id(6),error(6))
      endif
      if(comp(7).eq.1) then
         call hdfclose(file_id(7),group_id(7),error(7))
      endif
      if(comp(8).eq.1) then
         call hdfclose(file_id(8),group_id(8),error(8))
      endif
      if(comp(9).eq.1) then
         call hdfclose(file_id(9),group_id(9),error(9))
      endif
      
   endif 
end subroutine out_emf