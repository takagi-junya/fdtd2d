!-----------------------------------------------------------------------
!   完全導体
!-----------------------------------------------------------------------
subroutine PEC()   
    use constants 
    implicit none
    integer :: i,j
    real(kind=8) cy_rad
    if(obj.eq.1) then
            if(myrank.eq.0) then
                write(30,*)"ic:",ic,"jc:",jc
                write(30,'(a21)')"object type:cylinder"
                write(30,'(a12,f10.2,/)')"PEC radius:",radius
            endif
          do j=jstart,jend
              do i=istart,iend
                  if(cy_rad(i,j).le.radius.and.cy_rad(i+1,j).le.radius) then
                      aex(i,j)  =  0.0d0
                      bexy(i,j) =  0.0d0
                  endif
                  if(cy_rad(i,j).le.radius.and.cy_rad(i,j+1).le.radius) then
                      aey(i,j)  =  0.0d0
                      beyx(i,j) =  0.0d0
                  endif
                  if(cy_rad(i,j).le.radius)then
                      aez(i,j)    =  0.0d0
                      bezx(i,j)   =  0.0d0
                      bezy(i,j)   =  0.0d0
                  endif
              end do
          end do
    else if(obj==2) then
        write(30,'(a20,/)')"object type:prism"
        do j=jc-ly2,jc+ly2
            do i=ic-lx2,ic+lx2-1
              aex(i,j) = 0.0d0
              bexy(i,j) = 0.0d0
            enddo
        enddo
        do j=jc-ly2,jc+ly2-1
            do i=ic-lx2,ic+lx2
              aey(i,j) = 0.0d0
              beyx(i,j) = 0.0d0
            enddo
        enddo
        do j=jc-ly2,jc+ly2
            do i=ic-lx2,ic+lx2
              aez(i,j) = 0.0d0
              bezx(i,j) = 0.0d0
              bezy(i,j) = 0.0d0
            enddo
        enddo
    endif
end subroutine PEC