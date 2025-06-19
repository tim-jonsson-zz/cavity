module getIndV
  implicit none
contains 
subroutine para2IndV(Nin,Ns,para)
use matind
use convert
  integer, intent(in) :: Nin,Ns
  real*8, intent(in), allocatable :: para(:)    ! win, Deltal, gin,g2; 1,Ns,1,1
  integer*8 :: Ndim,i,jj
  integer :: si,sj,jin,jout,Nd
  integer,allocatable :: v1(:),v2(:),v3(:)     ! isl,iin,iout
  real*8 :: vtemp, gin, g2
  real*8,allocatable :: Deltal(:)
! Hmat values: win b+ b + Delta/2 sum_n sigma_n^z
! Hmat ind & values: gin (b^+ + b) sum_n sigma_n^x
! Hmat ind & values: g2 (b^+ + b) 
! only upper half jj>ii
! s1,s2,...,sn = 0,1; iin = 0,1,...,Nin
  Ndim = dimnsnd(Nin,Ns)

  allocate(v1(Ns+1),v2(Ns+1),v3(Ns+1))
  
  gin = para(Ns+2)
  g2 = para(Ns+3)
  i_of_b = 1
  i_of_d = 1
  do i = 1,Ndim
    v1 = Ind2Setnsnd(i,Nin,Ns)    
    do si = 1,Ns
      ! b+ sigma^x
        v2 = v1
        v2(si) = 1-v1(si)
        jin = v1(Ns+1) + 1
        if (jin <= Nin) then
          v2(Ns+1) = jin
          jj = Set2Indnsnd(Nin,Ns,v2)
          bInd1(i_of_b) = jj
          bInd2(i_of_b) = i
          bv(i_of_b) = sqrt(jin*1.0D0)*gin
          i_of_b = i_of_b + 1
        end if      
    end do
    v3 = v1
    jin = v1(Ns+1) + 1
    if (jin <= Nin) then
      v3(Ns+1) = jin
      jj = Set2Indnsnd(Nin,Ns,v3)
      dInd1(i_of_d) = jj
      dInd2(i_of_d) = i
      dv(i_of_d) = sqrt(jin*1.0D0)*g2
      i_of_d = i_of_d + 1
    end if      
  end do
  deallocate(v1,v2,v3)
  return
end subroutine para2IndV
!--------------------------------------------------------- 
subroutine para2V(Nin,Ns,para)
use matind
use convert
  integer, intent(in) :: Nin,Ns
  real*8, intent(in), allocatable :: para(:)    ! win, Deltal, gin, g2; 1,Ns,1, 1
  integer*8 :: Ndim,i,jj
  integer :: si,sj,jin,jout,Nd
  integer,allocatable :: v1(:),v2(:),v3(:)     ! isl,iin,iout
  real*8 :: vtemp, gin, g2
  real*8,allocatable :: Deltal(:)
! Hmat values: win b+ b + Delta/2 sum_n sigma_n^z 
! Hmat ind & values: gin (b^+ + b) sum_n sigma_n^x
! only upper half jj>ii
! s1,s2,...,sn = 0,1; iin = 0,1,...,Nin
  Ndim = dimnsnd(Nin,Ns)
  allocate(v1(Ns+1),v2(Ns+1),Deltal(Ns))
  Deltal = para(2:Ns+1)
  
  gin = para(Ns+2)

  do i = 1,Ndim
    v1 = Ind2Setnsnd(i,Nin,Ns)    
    vtemp = dot_product(v1(1:Ns)-0.5D0, Deltal)
    vtemp = vtemp + v1(Ns+1) * para(1)
    diagv(i) = vtemp
  end do
  deallocate(v1,v2,Deltal)
  return
end subroutine para2V
!--------------------------------------------------------- 
subroutine getGs(Nin,Ns)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  integer*8 :: Ndim,i,jj
  integer :: Nd
  integer, allocatable :: v1(:)

  Ndim = dimnsnd(Nin,Ns)
  allocate(v1(Ns+1))
  
  exN = 0.0
  do i=1,Ndim
    v1 = Ind2Setnsnd(i,Nin,Ns)
    exN = exN + v1 * hmat(i,1)**2
  end do
  deallocate(v1)
  return
end subroutine getGs
!---------------------------------------------------------   
subroutine getExT(Nin,Ns)
use convert
use Lanc_span
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  integer*8 :: Ndim,i,jj
  integer :: Nd
  integer, allocatable :: v1(:)

  Ndim = dimnsnd(Nin,Ns)
  allocate(v1(Ns+1))
  
  exN = 0.0
  do i=1,Ndim
    v1 = Ind2Setnsnd(i,Nin,Ns)
    exN = exN + v1 * abs(x_of_t(i))**2.D0
  end do
  deallocate(v1)
  return
end subroutine getExT
!---------------------------------------------------------   
subroutine getRhoT(Nin,Ns)
use convert
use Lanc_span
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  integer*8 :: Ndim,ii,jj
  integer*8 :: ri,rj
  integer :: ni
  integer, allocatable :: v1(:),v2(:)

  Ndim = dimnsnd(Nin,Ns)
  allocate(v1(Ns+1),v2(Ns+1))
  
  rhot = 0.0
  do ri=1,2**Ns
    v1(1:Ns) = Ind2Setns(ri,Ns)
    do rj=1,2**Ns
      v2(1:Ns) = Ind2Setns(rj,Ns)
      do ni=0,Nin
        v1(Ns+1) = ni
        v2(Ns+1) = ni
        ii = Set2Indnsnd(Nin,Ns,v1)
        jj = Set2Indnsnd(Nin,Ns,v2)
        rhot(ri,rj) = rhot(ri,rj) + x_of_t(ii)*conjg(x_of_t(jj))
      end do
    end do  
  end do
  deallocate(v1,v2)
  return
end subroutine getRhoT
!--------------------------------------------------------- 
subroutine getRho(Nin,Ns)
use convert
use Lanc_span
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  integer*8 :: Ndim,ii,jj
  integer*8 :: ri,rj
  integer :: ni
  integer, allocatable :: v1(:),v2(:)

  Ndim = dimnsnd(Nin,Ns)
  allocate(v1(Ns+1),v2(Ns+1))
  
  rho = 0.0
  do ri=1,2**Ns
    v1(1:Ns) = Ind2Setns(ri,Ns)
    do rj=1,2**Ns
      v2(1:Ns) = Ind2Setns(rj,Ns)
      do ni=0,Nin
        v1(Ns+1) = ni
        v2(Ns+1) = ni
        ii = Set2Indnsnd(Nin,Ns,v1)
        jj = Set2Indnsnd(Nin,Ns,v2)
        rho(ri,rj) = rho(ri,rj) + hmat(ii,1)*hmat(jj,1)
      end do
    end do  
  end do
  deallocate(v1,v2)
  return
end subroutine getRho
!--------------------------------------------------------- 
subroutine init_coherent(Nin,Ns,alpha)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  real*8, intent(in) :: alpha
  integer*8 :: Ndim,jj
  integer :: i,Nd
  integer, allocatable :: v1(:)
  real*8 :: C0

  Ndim = dimnsnd(Nin,Ns)
  C0 = exp(-alpha**2/2.D0)
  allocate(v1(Ns+1))
  
  v1 = 0
  do i=0,Nin
    v1(Ns+1) = i
    jj = Set2Indnsnd(Nin,Ns,v1)
    psi0(jj) = C0 * alpha**i / sqrt(gamma(real(i)+1.D0))
  end do
  return
end subroutine init_coherent
!--------------------------------------------------------- 
subroutine init_S0coha(Nin,Ns,alpha)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  real*8, intent(in) :: alpha
  integer*8 :: Ndim,jj
  integer :: i,si,Nd,S0a(4,6)
  integer, allocatable :: v1(:)
  real*8 :: C0,Csa(6)

  Ndim = dimnsnd(Nin,Ns)
  C0 = exp(-alpha**2/2.D0)
  ! s part: 1/sqrt(12) [|1100> -2|1010> +|1001> +|0110> -2|0101> +|0011>]
  S0a = 0
  S0a(1,1) = 1
  S0a(2,1) = 1
  S0a(1,2) = 1
  S0a(3,2) = 1
  S0a(1,3) = 1
  S0a(4,3) = 1
  S0a(2,4) = 1
  S0a(3,4) = 1
  S0a(2,5) = 1
  S0a(4,5) = 1
  S0a(3,6) = 1
  S0a(4,6) = 1
  Csa(1) = 1.D0
  Csa(2) = -2.D0
  Csa(3) = 1.D0
  Csa(4) = 1.D0
  Csa(5) = -2.D0
  Csa(6) = 1.D0
  Csa = Csa / sqrt(12.D0)
  allocate(v1(Ns+1))
  
  v1 = 0
  do i=0,Nin
    v1(Ns+1) = i
    do si=1,6
      v1(1:Ns) = S0a(:,si)
      jj = Set2Indnsnd(Nin,Ns,v1)
      psi0(jj) = Csa(si) * C0 * alpha**i / sqrt(gamma(real(i)+1.D0))
    end do  
  end do
  return
end subroutine init_S0coha
!--------------------------------------------------------- 
subroutine init_S0ia(Nin,Ns,i_init)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns,i_init
  integer*8 :: Ndim,jj
  integer :: i,si,Nd,S0a(4,6)
  integer, allocatable :: v1(:)
  real*8 :: C0,Csa(6)

  Ndim = dimnsnd(Nin,Ns)
  ! s part: 1/sqrt(12) [|1100> -2|1010> +|1001> +|0110> -2|0101> +|0011>]
  S0a = 0
  S0a(1,1) = 1
  S0a(2,1) = 1
  S0a(1,2) = 1
  S0a(3,2) = 1
  S0a(1,3) = 1
  S0a(4,3) = 1
  S0a(2,4) = 1
  S0a(3,4) = 1
  S0a(2,5) = 1
  S0a(4,5) = 1
  S0a(3,6) = 1
  S0a(4,6) = 1
  Csa(1) = 1.D0
  Csa(2) = -2.D0
  Csa(3) = 1.D0
  Csa(4) = 1.D0
  Csa(5) = -2.D0
  Csa(6) = 1.D0
  Csa = Csa / sqrt(12.D0)
  allocate(v1(Ns+1))
  
  v1 = 0

    v1(Ns+1) = i_init
    do si=1,6
      v1(1:Ns) = S0a(:,si)
      jj = Set2Indnsnd(Nin,Ns,v1)
      psi0(jj) = Csa(si) 
    end do  

  return
end subroutine init_S0ia
!--------------------------------------------------------- 
subroutine init_S0cohb(Nin,Ns,alpha)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  real*8, intent(in) :: alpha
  integer*8 :: Ndim,jj
  integer :: i,si,Nd,S0b(4,4)
  integer, allocatable :: v1(:)
  real*8 :: C0,Csb(4)

  Ndim = dimnsnd(Nin,Ns)
  C0 = exp(-alpha**2/2.D0)
  ! s part:  1/2 [|1100> -|1001> -1/2|0110> +1/2|0011>]
  S0b = 0
  S0b(1,1) = 1
  S0b(2,1) = 1
  S0b(1,2) = 1
  S0b(4,2) = 1
  S0b(2,3) = 1
  S0b(3,3) = 1
  S0b(3,4) = 1
  S0b(4,4) = 1
  Csb(1) = 0.5D0
  Csb(2) = -0.5D0
  Csb(3) = -0.5D0
  Csb(4) = 0.5D0
  allocate(v1(Ns+1))
  
  v1 = 0
  do i=0,Nin
    v1(Ns+1) = i
    do si=1,4
      v1(1:Ns) = S0b(:,si)
      jj = Set2Indnsnd(Nin,Ns,v1)
      psi0(jj) = Csb(si) * C0 * alpha**i / sqrt(gamma(real(i)+1.D0))
    end do  
  end do
  return
end subroutine init_S0cohb
!---------------------------------------------------------
subroutine getPSI1(Nin,Ns,alpha)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  real*8, intent(in) :: alpha
  integer*8 :: Ndim,jj
  integer :: i,Nd
  integer, allocatable :: v1(:)
  real*8 :: C0
  ! |psi0> ~ |alpha> + |-alpha> 
  !        = 1/sqrt(cosh(alpha^2)) alpha^n / sqrt(n!) , n=0,2,...
  Ndim = dimnsnd(Nin,Ns)
  C0 = 1/sqrt(cosh(alpha**2))
  allocate(v1(Ns+1))
  
  v1 = 0
  do i=0,Nin,2
    v1(Ns+1) = i
    jj = Set2Indnsnd(Nin,Ns,v1)
    psi0(jj) = C0 * alpha**i / sqrt(gamma(real(i)+1.D0))
  end do
  return
end subroutine getPSI1
!---------------------------------------------------------
subroutine getPSI2(Nin,Ns,alpha)
use convert
use NAGmodule
use exNmodule
  integer, intent(in) :: Nin,Ns
  real*8, intent(in) :: alpha
  integer*8 :: Ndim,jj
  integer :: i,Nd
  integer, allocatable :: v1(:)
  real*8 :: C0
  ! |psi0> ~ |alpha> - |-alpha> 
  !        = 1/sqrt(sinh(alpha^2)) alpha^n / sqrt(n!) , n=0,2,...
  Ndim = dimnsnd(Nin,Ns)
  C0 = 1/sqrt(sinh(alpha**2))
  allocate(v1(Ns+1))
  
  v1 = 0
  do i=1,Nin,2
    v1(Ns+1) = i
    jj = Set2Indnsnd(Nin,Ns,v1)
    psi0(jj) = C0 * alpha**i / sqrt(gamma(real(i)+1.D0))
  end do
  return
end subroutine getPSI2
!--------------------------------------------------------- 
end module getIndV
