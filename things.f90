      module  exNmodule
        real*8, allocatable :: exN(:), psi0(:), rho(:,:)
        complex*16, allocatable :: rhot(:,:)
      end module exNmodule

      module  NAGmodule
        real*8, allocatable :: hmat(:,:),eval(:)
      end module NAGmodule
      
      module matind
      implicit none
        integer*8, allocatable :: bInd1(:),bInd2(:),dInd1(:),dInd2(:),b2Ind1(:),b2Ind2(:)   
        real *8, allocatable :: diagv(:),bv(:),dv(:),b2v(:)   
        integer*8 :: i_of_b,i_of_b2,i_of_d               
      end module matind
      
      module  KrylovH
          real*8,    allocatable :: tridia(:,:)
          real*8,    allocatable :: autoval(:)
      end module KrylovH
      
      module Lanc_span 
          complex*16, allocatable :: x_0(:)
          complex*16, allocatable :: x_km1(:),x_k(:),Hx_k(:),x_kp1(:),y_temp(:)
          complex*16, allocatable :: c_of_t(:),x_of_t(:),A_base(:,:),phases(:)
      end module Lanc_span
      
      module timevolution
          integer*8 :: ndiag,klimit
          real*8 :: ht,ft      
      end module timevolution  


      module things
        real*8 :: PatT
      end module things
      
      module path_names
        character (len=128) :: rpath,paraname,fname
      end module

