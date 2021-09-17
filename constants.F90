module constants
   use HDF5              
    
   !/space/
   integer :: nxx,nyy           
   real(kind=8) :: dx,dy                 
   integer :: abc
   integer :: pbc(3),lpml(3)

   !/time/
   real(kind=8) :: t
   real(kind=8) :: dt
   real(kind=8) :: deltat                      
   integer :: step,nstep                   

   !/output/
   integer :: out                      
   integer :: ostart
   integer :: odom(4)
   integer :: stride                
   integer :: comp(9)              
   integer :: io,jo             
                
   !/scat/
   integer :: mode            
   integer ::lx,ly              
   real(kind=8) :: gamma0,theta0,phi0,amp
   real(kind=8) :: lambda,tau0,freq,omega

   !/far/
   integer :: isx,isy
   real(kind=8) :: theta1,phi1
   
   !/object/
   integer :: obj              
   integer :: med       
   integer :: ic,jc                    
   integer :: lx2,ly2                    
   real(kind=8) :: epsr,radius                 
   
   !/feed/
   integer :: ip,jp,kp
   real(kind=8) :: duration,t0
   
   !/wave/
   integer :: kwave
   real(kind=8) :: tau
   real(kind=8) :: amps(6)
   real(kind=8) :: orgs(3)
   real(kind=8) :: ang(3)              
   real(kind=8) :: pc,pt,pw
   real(kind=8) :: alpha

   !/plasma
   integer :: pls
   real(kind=8) :: prad,erad
   real(kind=8) :: nu,wp

   integer ::ifed,jfed,kfed
   !PML
   integer,parameter :: order=4
   real(kind=8),parameter :: rmax=-120.0d0

   !全空間セル数
   integer :: nx,ny
   
   real(kind=8),parameter::epsbk=1.0d0,mubk=1.0d0,sigebk=0.0d0,sigmbk=0.0d0

   !HDF5
   integer :: h5count

    
   
   real(kind=8) :: mux,muy,muz
   real(kind=8) :: a,epsx,epsy,epsz,sgmx,sgmy,sgmz,sgmm
   real(kind=8) :: sgex,sgey,sgez
   real(kind=8),parameter::eps0=8.854187817d-12,mu0=1.2566370614d-6 
   real(kind=8),parameter::qe=1.602176462d-19,mel = 9.10938188d-31
   real(kind=8),parameter::c=2.99792458d8,z0=376.73031d0      !
   real(kind=8),parameter::radi0=1.74532925d-2               
   real(kind=8),parameter::pai=3.141592653589                
   real(kind=8),allocatable :: ex(:,:),ey(:,:),ez(:,:)
   real(kind=8),allocatable :: hx(:,:),hy(:,:),hz(:,:)
   real(kind=8),allocatable :: jx(:,:),jy(:,:),jz(:,:)
   real(kind=8),allocatable :: vx(:,:),vy(:,:),vz(:,:)
   real(kind=8),allocatable :: nd(:,:)
    
   real(kind=8),allocatable :: aex(:,:),aey(:,:),aez(:,:)
   real(kind=8),allocatable :: bexy(:,:),bexz(:,:)
   real(kind=8),allocatable :: beyx(:,:),beyz(:,:)
   real(kind=8),allocatable :: bezx(:,:),bezy(:,:)

   real(kind=8),allocatable :: amx(:,:),amy(:,:),amz(:,:)
   real(kind=8),allocatable :: bmxy(:,:),bmxz(:,:)
   real(kind=8),allocatable :: bmyx(:,:),bmyz(:,:)
   real(kind=8),allocatable :: bmzx(:,:),bmzy(:,:)

   real(kind=8),allocatable :: epsd(:,:),sgmed(:,:),mud(:,:),sgmmd(:,:)

   real(kind=8) :: ajj
   real(kind=8),allocatable :: ajx(:,:),ajy(:,:),ajz(:,:)
   real(kind=8),allocatable :: avx(:,:),avy(:,:),avz(:,:)
   real(kind=8),allocatable :: ajex(:,:),ajey(:,:),ajez(:,:)

    
end module constants