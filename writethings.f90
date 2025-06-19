module writethings
contains
      subroutine dirnames(Nin,Ns,Nout,para,alpha)
      use path_names
      implicit none
      integer, intent(in) :: Nin,Ns
      integer, intent(in), allocatable :: Nout(:)
      real*8, intent(in), allocatable :: para(:)    ! win, Delta, woutl, gin
      real*8, intent(in) :: alpha
      
      rpath = "/export/zz/cavity/a--a/"   ! v3-n2
      if (Nin>99 .and. Nout(1) >9) then
      write(paraname, &
      "('Na',I3,'Ns',I1,'Nb',I2,'_wa',F3.1,'_D',F3.1, '_alpha',F3.1)") Nin,Ns,Nout(1),para(1),para(2),alpha
      else if (Nin>9 .and. Nout(1) >9) then
      write(paraname, &
      "('Na',I2,'Ns',I1,'Nb',I2,'_wa',F3.1,'_D',F3.1, '_alpha',F3.1)") Nin,Ns,Nout(1),para(1),para(2),alpha
      else if (Nin<10 .and. Nout(1) <10) then
      write(paraname, &
      "('Na',I1,'Ns',I1,'Nb',I1,'_wa',F3.1,'_D',F3.1, '_alpha',F3.1)") Nin,Ns,Nout(1),para(1),para(2),alpha
      end if
      
      write(fname,"('/woa',F5.3,'.txt')") para(3)
      
      return
      end subroutine
end module writethings    
