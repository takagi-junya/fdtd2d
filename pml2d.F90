module pml2d 
   use constants
   implicit none
   type pml                                           
      integer::i0,i1,j0,j1
      real(kind=8),pointer::expml(:,:),eypml(:,:),ezx(:,:),ezy(:,:)
      real(kind=8),pointer::hxpml(:,:),hypml(:,:),hzx(:,:),hzy(:,:) 
      real(kind=8),pointer::aexpml(:,:),aeypml(:,:)
      real(kind=8),pointer::beypml(:,:),bexpml(:,:)
      real(kind=8),pointer::amxpml(:,:),amypml(:,:)
      real(kind=8),pointer::bmypml(:,:),bmxpml(:,:)
   end type pml
   type(pml)::pml_l,pml_r,pml_d,pml_u
   real(kind=8)::copml
   parameter(copml=-1.5280063d-4)
   contains
   
   subroutine initpml()
      use constants
      implicit none
      if(myrank.eq.0) then 
         write(30,*)"init_pml"
         write(30,*)"pmlx:",lpml(1),"pmly:",lpml(2)
         write(30,*)""
      endif
      if(coords(1).eq.0) then
         call init_pml(pml_l,0,lpml(1),jstart-1,jend+1)
      endif
      if(coords(1).eq.ndims(1)-1) then
         call init_pml(pml_r,nx-lpml(1),nx,jstart-1,jend+1)
      endif
      if(coords(2).eq.0) then
         call init_pml(pml_d,istart-1,iend+1,0,lpml(2))
      endif
      if(coords(2).eq.ndims(2)-1) then
         call init_pml(pml_u,istart-1,iend+1,ny-lpml(2),ny)
      endif
   end subroutine initpml

   subroutine init_pml(p,i0,i1,j0,j1)
      use constants
      implicit none
      type(pml)::p
      integer,intent(in) :: i0,i1,j0,j1
      integer :: i,j
      real(kind=8)::a
      real(kind=8)::sigmxm,sigmxe
      real(kind=8)::sigmym,sigmye
      real(kind=8)::smax0x,smax0y
      real(kind=8)::epspml,mupml

      p%i0=i0
      p%i1=i1
      p%j0=j0
      p%j1=j1

      allocate(p%expml(i0:i1,j0:j1))   
      allocate(p%eypml(i0:i1,j0:j1))
      allocate(p%ezx(i0:i1,j0:j1))
      allocate(p%ezy(i0:i1,j0:j1))
      allocate(p%hxpml(i0:i1,j0:j1))
      allocate(p%hypml(i0:i1,j0:j1))
      allocate(p%hzx(i0:i1,j0:j1))
      allocate(p%hzy(i0:i1,j0:j1))

      p%expml=0.0d0                       
      p%eypml=0.0d0
      p%ezx=0.0d0
      p%ezy=0.0d0
      p%hxpml=0.0d0
      p%hypml=0.0d0
      p%hzx=0.0d0
      p%hzy=0.0d0

      allocate(p%aeypml(i0:i1,j0:j1))
      allocate(p%aexpml(i0:i1,j0:j1))
      allocate(p%amypml(i0:i1,j0:j1))
      allocate(p%amxpml(i0:i1,j0:j1))
      allocate(p%beypml(i0:i1,j0:j1))
      allocate(p%bexpml(i0:i1,j0:j1))
      allocate(p%bmypml(i0:i1,j0:j1))
      allocate(p%bmxpml(i0:i1,j0:j1))

      smax0x=copml*rmax*(order+1)/(lpml(1)*dx) 
      smax0y=copml*rmax*(order+1)/(lpml(2)*dy)
      
      mupml=mubk*mu0                                     
      epspml=epsbk*eps0
      do i=i0,i1
         do j=j0,j1          
            if(i<lpml(1)) then                             
               sigmxm=((lpml(1)-i-0.5d0)/lpml(1))**order*smax0x      
               sigmxe=(float(lpml(1)-i)/lpml(1))**order*smax0x
            else if(i>=nx-lpml(1)) then                     
               sigmxm=((i-nx+lpml(1)+0.5d0)/lpml(1))**order*smax0x
               sigmxe=(float(i-nx+lpml(1))/lpml(1))**order*smax0x
            else                                         
               sigmxm=0.0d0
               sigmxe=0.0d0
            end if
         
            if(j<lpml(2)) then                              
               sigmym=((lpml(2)-j-0.5d0)/lpml(2))**order*smax0y
               sigmye=(float(lpml(2)-j)/lpml(2))**order*smax0y
            else if(j>=ny-lpml(2)) then                     
               sigmym=((j-ny+lpml(2)+0.5d0)/lpml(2))**order*smax0y
               sigmye=(float(j-ny+lpml(2))/lpml(2))**order*smax0y
            else                                         
               sigmym=0.0d0
               sigmye=0.0d0
            end if

            sigmxe=sigmxe*epsbk
            a=0.5d0*sigmxe*dt/epspml
            p%aexpml(i,j)=(1.0d0-a)/(1.0d0+a)
            p%bexpml(i,j)=dt/epspml/(1.0d0+a)/dx
            
            sigmye=sigmye*epsbk
            a=0.5d0*sigmye*dt/epspml
            p%aeypml(i,j)=(1.0d0-a)/(1.0d0+a)
            p%beypml(i,j)=dt/epspml/(1.0d0+a)/dy

            sigmxm=sigmxm*epsbk
            a=0.5d0*sigmxm*dt/epspml
            p%amxpml(i,j)=(1.0d0-a)/(1.0d0+a)
            p%bmxpml(i,j)=dt/mupml/(1.0d0+a)/dx

            sigmym=sigmym*epsbk
            a=0.5d0*sigmym*dt/epspml
            p%amypml(i,j)=(1.0d0-a)/(1.0d0+a)
            p%bmypml(i,j)=dt/mupml/(1.0d0+a)/dy
           end do
      end do
   end subroutine init_pml
   
   subroutine epml()
        
      use constants
      implicit none

      if(coords(1).eq.0) then
         call exchg2dj(hz(istart-1:iend+1,jstart-1:jend+1),lx0)
         call exchg2dj(hx(istart-1:iend+1,jstart-1:jend+1),lx0)
         call e_pml(pml_l)
      endif
      if(coords(1).eq.ndims(1)-1) then
         call exchg2dj(hz(istart-1:iend+1,jstart-1:jend+1),lx0)
         call exchg2dj(hx(istart-1:iend+1,jstart-1:jend+1),lx0)
         call e_pml(pml_r)                        !右側のPML
      endif
      if(coords(2).eq.0) then
         call exchg2di(hz(istart-1:iend+1,jstart-1:jend+1))
         call exchg2di(hy(istart-1:iend+1,jstart-1:jend+1))
         call e_pml(pml_d)                        !下部PML
      endif
      if(coords(2).eq.ndims(2)-1) then
         call exchg2di(hz(istart-1:iend+1,jstart-1:jend+1))
         call exchg2di(hy(istart-1:iend+1,jstart-1:jend+1))
         call e_pml(pml_u)                        !上部PML
      endif
   end subroutine epml 

   subroutine e_pml(p)
   
      use constants
      implicit none
      type(pml)::p
      integer::i,j,i0p,i1p,j0p,j1p

      i0p=p%i0               
      i1p=p%i1
      j0p=p%j0                
      j1p=p%j1
      !Ex
      do j=j0p+1,j1p-1
         do i=i0p,i1p-1 
            p%expml(i,j)=p%aeypml(i,j)*p%expml(i,j)&
      &                 +p%beypml(i,j)*(hz(i,j)-hz(i,j-1))
            ex(i,j)=p%expml(i,j)
         end do   
      end do
      !Ey
      do j=j0p,j1p-1   
         do i=i0p+1,i1p-1
            p%eypml(i,j)=p%aexpml(i,j)*p%eypml(i,j)&
      &                 -p%bexpml(i,j)*(hz(i,j)-hz(i-1,j))
            ey(i,j)=p%eypml(i,j)
         end do
      end do
      !Ez
      do j=j0p+1,j1p-1
         do i=i0p+1,i1p-1 
            p%ezx(i,j)=p%aexpml(i,j)*p%ezx(i,j)&
      &                +p%bexpml(i,j)*(hy(i,j)-hy(i-1,j))
            p%ezy(i,j)=p%aeypml(i,j)*p%ezy(i,j)&
      &                -p%beypml(i,j)*(hx(i,j)-hx(i,j-1))
            ez(i,j)=p%ezx(i,j)+p%ezy(i,j)
         end do
      end do
   end subroutine e_pml

   subroutine hpml()     
      use constants
      implicit none
      if(coords(1).eq.0) then
         call exchg2dj(ez(istart-1:iend+1,jstart-1:jend+1),lx0)
         call exchg2dj(ex(istart-1:iend+1,jstart-1:jend+1),lx0)
         call h_pml(pml_l)                      !左側のPML
      endif
      if(coords(1).eq.ndims(1)-1) then
         call exchg2dj(ez(istart-1:iend+1,jstart-1:jend+1),lx0)
         call exchg2dj(ex(istart-1:iend+1,jstart-1:jend+1),lx0)
         call h_pml(pml_r)                      !右側のPML
      endif
      if(coords(2).eq.0) then
         call exchg2di(ez(istart-1:iend+1,jstart-1:jend+1))
         call exchg2di(ey(istart-1:iend+1,jstart-1:jend+1))
         call h_pml(pml_d)                      !下部PML
      endif

      if(coords(2).eq.ndims(2)-1) then
         call exchg2di(ez(istart-1:iend+1,jstart-1:jend+1))
         call exchg2di(ey(istart-1:iend+1,jstart-1:jend+1))
         call h_pml(pml_u)                      !上部PML
      endif
   end subroutine hpml

   subroutine h_pml(p)      
   
      use constants
      implicit none 
      type(pml)::p     
      integer::i,j,i0p,i1p,j0p,j1p
      
      i0p=p%i0           
      i1p=p%i1
      j0p=p%j0             
      j1p=p%j1
      !Hx   
      do j=j0p,j1p-1
         do i=i0p+1,i1p-1
            p%hxpml(i,j)=p%amypml(i,j)*p%hxpml(i,j)&
      &                  -p%bmypml(i,j)*(ez(i,j+1)-ez(i,j))
            hx(i,j)=p%hxpml(i,j)
         end do
      end do
      !Hy
      do j=j0p+1,j1p-1
         do i=i0p,i1p-1
            p%hypml(i,j)=p%amxpml(i,j)*p%hypml(i,j)&
      &                  +p%bmxpml(i,j)*(ez(i+1,j)-ez(i,j))
            hy(i,j)=p%hypml(i,j)
         end do
      end do
      !Hz
      do j=j0p,j1p-1
         do i=i0p,i1p-1     
            p%hzx(i,j)=p%amxpml(i,j)*p%hzx(i,j)&
      &                 -p%bmxpml(i,j)*(ey(i+1,j)-ey(i,j))
            p%hzy(i,j)=p%amypml(i,j)*p%hzy(i,j)&
      &                 +p%bmypml(i,j)*(ex(i,j+1)-ex(i,j))
            hz(i,j)=p%hzx(i,j)+p%hzy(i,j)
         end do
      end do
   end subroutine h_pml
end module
  