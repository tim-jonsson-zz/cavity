module convert
  implicit none
contains
integer*8 function Set2Indnsnd(Nin,Ns,vi)
  integer, intent(in) :: Nin,Ns
  integer, intent(in), allocatable :: vi(:)  ! (s1,s2,...,sn,nin)
  integer :: i,Nd
  integer*8 :: Ind

  ! Ind = s1 + 2*s2 + ... + 2^(Ns-1)*sn + iin*2^Ns 
  !          + id1*2^n*(Nin+1) 
  Ind = 1 
  do i=1,Ns
    Ind = Ind + 2**(i-1)*vi(i)
  end do
  Ind = Ind + vi(Ns+1)*2**Ns 
  Set2Indnsnd = Ind
  return
end function Set2Indnsnd 

integer*8 function dimnsnd(Nin,Ns)
  integer, intent(in) :: Nin,Ns

  dimnsnd = (Nin+1) * 2**Ns
  return
end function dimnsnd  

function Ind2Setnsnd(Ind,Nin,Ns) result(output)
  integer*8, intent(in) :: Ind
  integer, intent(in) :: Nin,Ns
  integer :: i,Nd
  integer*8 :: temp
  integer, allocatable :: output(:)
! Ind = s1 + 2*s2 + ... + 2^(Ns-1)*sn + iin*2^Ns 
!          + id1*2^n*(Nin+1) + id2*2^n*(Nin+1)*(Nd1+1) + ... + idn*2^n*(Nin+1)*(Nd1+1)*...*
  allocate(output(Ns+1))
  temp = Ind - 1
  do i =1,Ns
    output(i) = mod(temp,2)
    temp = temp / 2
  end do  
  output(Ns+1) = mod(temp,Nin+1)
end function Ind2Setnsnd

integer*8 function Set2Indns(Ns,vi)
  integer, intent(in) :: Ns
  integer, intent(in), allocatable :: vi(:)  ! (s1,s2,...,sn)
  integer :: i
  integer*8 :: Ind

  ! Ind = s1 + 2*s2 + ... + 2^(Ns-1)*sn 
  Ind = 1 
  do i=1,Ns
    Ind = Ind + 2**(i-1)*vi(i)
  end do
  Set2Indns = Ind
  return
end function Set2Indns

function Ind2Setns(Ind,Ns) result(output)
  integer*8, intent(in) :: Ind
  integer, intent(in) :: Ns
  integer :: i,Nd
  integer*8 :: temp
  integer, allocatable :: output(:)
! Ind = s1 + 2*s2 + ... + 2^(Ns-1)*sn 
  allocate(output(Ns))
  temp = Ind - 1
  do i =1,Ns
    output(i) = mod(temp,2)
    temp = temp / 2
  end do  
end function Ind2Setns

end module convert
