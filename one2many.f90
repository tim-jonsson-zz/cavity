module one2many
implicit none
contains
function wi2set(wi,Nwl) result(output)
  integer, intent(in) :: wi
  integer, intent(in), allocatable :: Nwl(:)
  integer :: i,Nd
  integer :: temp
  integer, allocatable :: output(:)
  
  Nd = size(Nwl)
  allocate(output(Nd))
  temp = wi-1
  do i=1,Nd
    output(i) = mod(temp,Nwl(i))
    temp = temp / Nwl(i)
  end do
end function  
end module one2many
