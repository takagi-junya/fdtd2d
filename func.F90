real(kind=8) function pscat(x)
    use constants
    implicit none
    real(kind=8),intent(in) :: x
    real(kind=8) :: xx,aa
    xx = x-c*t-pc
    aa = (xx/pw)*(xx/pw)
    pscat = amp*exp(-aa)
end function pscat

real(kind=8) function pscat2(x)
    use constants
    implicit none
    real(kind=8),intent(in) :: x
    real(kind=8) :: xx,aa
    xx = x-c*t-pc
    aa = (xx/pw)*(xx/pw)
    pscat2 = amp*exp(-aa)*cos(2.0d0*pai*freq/c*xx)
end function pscat2


real(kind=8) function sp_rad(x,y)
   use constants
   implicit none
   integer,intent(in) :: x,y
   real(kind=8) :: xx,yy
   xx = x
   yy = y
   sp_rad = sqrt((xx-ic)**2+(yy-jc)**2)
end function sp_rad

real(kind=8) function cy_rad(x,y)
   use constants
   implicit none
   integer,intent(in) :: x,y
   real(kind=8) :: xx,yy
   xx = x
   yy = y
   cy_rad = sqrt((xx-ic)**2+(yy-jc)**2)
end function cy_rad