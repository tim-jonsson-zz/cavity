program cavity_Tevo
use getIndV
use matind
use convert
use NAGmodule
use exNmodule
use sparseED
use KrylovH
use Lanc_span
use timevolution 
use decay
use things
use path_names
use writethings
implicit none
integer*8 :: i,Ind,Ndim,pi,pj,ii,jj
integer :: iin,iout,Nin,Ns,Nd,Ntstep,wi
integer, allocatable :: wl(:)
real*8, allocatable :: para(:)              ! win, Deltal, gin, g2; 1,Ns,1, 1
real*8 :: alpha,Gamma_out,qqqq
real*8 :: tt,t_start,t_finish,tw0,tw1
complex*16 :: trace_rho

call cpu_time(t_start)
Nin = 50
Ns = 4
allocate(para(Ns+3))

write(*,*) "Nin="
write(*,*) Nin
write(*,*) "Ns="
write(*,*) Ns


! win, Deltal, gin, g2
para(1) = 1.D0
para(2) = 1.D0
para(3) = 1.D0
para(4) = 1.D0
para(5) = 1.D0
para(Ns+2) = 0.1D0
para(Ns+3) = 0.D0

write(*,*) "win="
write(*,*) para(1)
write(*,*) "Delta="
write(*,*) para(2:(Ns+1))
write(*,*) "gin="
write(*,*) para(Ns+2)
write(*,*) "g2="
write(*,*) para(Ns+3)
! 
alpha = 1.0D0
Gamma_out = 0.02D0

! 
klimit = 10

Ndim = dimnsnd(Nin,Ns)
print '("H size",I10)', Ndim

allocate(psi0(Ndim))
allocate(exN(Ns+1),rhot(2**Ns,2**Ns))
allocate(diagv(Ndim))
allocate(bInd1(Ndim*Ns),bInd2(Ndim*Ns),bv(Ndim*Ns))
allocate(dInd1(Ndim),dInd2(Ndim),dv(Ndim))

call para2IndV(Nin,Ns,para)
call para2V(Nin,Ns,para)

call cpu_time(tw0)
print '("length=",I10)', i_of_b-1
print '("length=",I10)', i_of_d-1
print '("  Non-diag term time = ",f10.3," seconds.")',tw0-t_start

!! start with number state + S0
call init_S0ia(Nin,Ns,0)
write(*,*) "<Psi0|Psi0>:"
write(*,*) dot_product(psi0,psi0)

ndiag = Ndim
call make_alloc
  x_0 = psi0
  x_of_t = x_0
  call getExT(Nin,Ns)
  write(*,*) "exN(t0):"
  write(*,*) exN
  call getRhoT(Nin,Ns)
  write(*,*) "rho(t0):"
  open(1,file="data/t0_rho_ED_real.txt")
  open(2,file="data/t0_rho_ED_imag.txt")
  do i=1,2**Ns
    write(1,*) real(rhot(i,:))
    write(2,*) aimag(rhot(i,:))
    write(*,*) rhot(i,i)
    trace_rho = trace_rho + rhot(i,i)
  end do
  close(1)
  close(2)
  write(*,*) "Tr rho:"
  write(*,*) trace_rho
  !open(3,file="/export/zz/cavity/JC/rho_real.txt")
  !open(4,file="/export/zz/cavity/JC/rho_imag.txt")
  open(3,file="JC/rho_real.txt")
  open(4,file="JC/rho_imag.txt")
  write(3,*) real(rhot)
  write(4,*) aimag(rhot)
!! ground state done.

ht = 0.05D0
Ntstep = 2500
!ft = decay_factor(ht,Gamma_out)



!!!! T evo start
!para(2) = 0.5D0   ! Delta_1 -> 0.5

  call cpu_time(tw0)
 
  open(1,file="data/Nin_Nout_lanc.txt")
  do i=0,Ntstep-1
    tt = ht*(i+0.5D0)
    ft = 1.D0
    para(Ns+3) = 2.D0 * switchon_factor(tt,100.D0)
    call para2IndV(Nin,Ns,para)
    call mylanczos(ht,ft)
    call getExT(Nin,Ns)
    write(1,*) exN
    if (mod(i+1,10) == 0) then
      call getRhoT(Nin,Ns)
      write(3,*) real(rhot)
      write(4,*) aimag(rhot)
    end if
    x_0 = x_of_t
  end do
  close(1)
  
  ! after evolution
  write(*,*) "exN(tf):"
  write(*,*) exN
  call getRhoT(Nin,Ns)
  trace_rho = (0.D0,0.D0)
  write(*,*) "rho(tf):"
  open(1,file="data/tf_rho_ED_real.txt")
  open(2,file="data/tf_rho_ED_imag.txt")
  do i=1,2**Ns
    write(1,*) real(rhot(i,:))
    write(2,*) aimag(rhot(i,:))
    write(*,*) rhot(i,i)
    trace_rho = trace_rho + rhot(i,i)
  end do
  close(1)
  close(2)
  write(*,*) "Tr rho:"
  write(*,*) trace_rho
  call cpu_time(tw1)
  write(*,*) "One evolution done. Parameters:"
  write(*,*) para
  print '("  Time = ",f10.3," seconds.")',tw1-tw0
  write(*,*) ""

call cpu_time(t_finish)
print '("Total time = ",f20.3," seconds.")',t_finish-t_start
deallocate(psi0)
deallocate(diagv)
deallocate(bInd1,bInd2,bv,dInd1,dInd2,dv)
call de_alloc
end program cavity_Tevo


!-------------------------------------------------------------------------------   
      subroutine make_alloc
      use KrylovH
      use Lanc_span
      use timevolution          
      use nagmodule
      implicit none

      
      !deallocate(hmat,eval)
      allocate(x_0(ndiag))
      allocate(Hx_k(ndiag),x_km1(ndiag),x_k(ndiag),x_kp1(ndiag),y_temp(ndiag))
      allocate(tridia(klimit,klimit),autoval(klimit)) 
      allocate(c_of_t(klimit),x_of_t(ndiag),A_base(ndiag,klimit),phases(klimit)) 
      return
      end 
!-------------------------------------------------------------------------------  
      subroutine de_alloc
      use KrylovH
      use Lanc_span
      use timevolution          
      use nagmodule
      implicit none

      
      !deallocate(hmat,eval)
      deallocate(x_0)
      deallocate(Hx_k,x_km1,x_k,x_kp1,y_temp)
      deallocate(tridia,autoval) 
      deallocate(c_of_t,x_of_t,A_base,phases) 
      return
      end 
!---------------------------------------------------------   
      SUBROUTINE NAGdiag(ndiag)
      use NAGmodule
      
      implicit none
      integer*8 ndiag,info,lda,liwork,lwork,nmax
      
      real*8,  allocatable :: WORK(:)
      integer*8, allocatable :: IWORK(:)  
                   
      CHARACTER JOB, UPLO
      EXTERNAL dsyevd
      NMAX=ndiag
      LDA=NMAX
      LWORK=4*NMAX*NMAX+100
      LIWORK=5*NMAX
      LIWORK=10*NMAX      
      ALLOCATE(WORK(LWORK),IWORK(LIWORK))
      
      JOB='V'    
      UPLO='L' 
      
      CALL dsyevd(JOB,UPLO,ndiag,hmat,LDA,eval,WORK,LWORK,IWORK,LIWORK,INFO)
      
      IF (INFO.GT.0) THEN
      WRITE (*,*) 'Failure to converge.'
      stop
      endif

      deALLOCATE(WORK,IWORK)

      return
      end

!-------------------------------------------------------------------------------      

      SUBROUTINE mylanczos(htime,ftime)
      use KrylovH
      use Lanc_span
      use matind
      use NAGmodule
      use timevolution
      implicit none

      integer*8 i,j,l,ilamb,ix,iy,klimit0,m23,mm,nn
      real*8 htime,ftime,qqqq,a_k,b_k,b_kp1,dum,rr
      real*8, allocatable :: tridia0(:,:),autoval0(:)
      complex*16, allocatable :: phases0(:),c_of_t0(:),A_base0(:,:)
      complex*16 uc

      M23=ndiag
      uc=(0.,1.)
      tridia=0.

      x_k=x_0
      x_km1=0.
      b_k=0.
      A_base=0.

      klimit0=0
dolanc: do mm= 1,klimit

      	  if(mm>1) then
	    do nn=1,mm-1
	    dum=dot_product(A_base(:,nn),x_k)
	    x_k(:)=x_k(:)-dum*A_base(:,nn)
	    enddo
	  endif
	  dum=sqrt(real(dot_product(x_k,x_k)))
	  x_k=x_k/dum


          A_base(:,mm)=x_k(:)

	  Hx_k = diagv * x_k
	  
          DO L=1,i_of_b-1
            I=bInd1(L)
            J=bInd2(L)
            qqqq=bv(L)
            Hx_k(I)=Hx_k(I)+qqqq*x_k(J)
            Hx_k(J)=Hx_k(J)+qqqq*x_k(I)
          ENDDO
          
          DO L=1,i_of_d-1
            I=dInd1(L)
            J=dInd2(L)
            qqqq=dv(L)
            Hx_k(I)=Hx_k(I)+qqqq*x_k(J)
            Hx_k(J)=Hx_k(J)+qqqq*x_k(I)
          ENDDO                   

          a_k= dot_product(x_k,Hx_k)
          y_temp=Hx_k-a_k*x_k
          rr= dot_product(y_temp,y_temp)-b_k*b_k
	  if(rr<0.d0.or.abs(rr)<1.d-16)exit dolanc
	  klimit0=klimit0+1

          b_kp1=sqrt(rr)
          x_kp1=(y_temp-b_k*x_km1)/b_kp1
          tridia(mm,mm)=a_k
          if(mm < klimit)then
          tridia(mm,mm+1)=b_kp1
          tridia(mm+1,mm)=b_kp1
          endif
          x_km1 = x_k
          x_k   = x_kp1
          b_k=b_kp1
      enddo dolanc

      if(klimit0==0)then
      x_of_t=x_0*exp(-uc*htime*a_k)
      return
      endif

      allocate(hmat(klimit0,klimit0),eval(klimit0),tridia0(klimit0,klimit0),autoval0(klimit0))
      allocate(phases0(klimit0),c_of_t0(klimit0),A_base0(M23,klimit0))
      hmat=0.d0; eval=0.d0 

      do ix=1,klimit0
      do iy=1,klimit0
      hmat(ix,iy)=tridia(ix,iy)
      enddo
      enddo

      do ix=1,klimit0
      A_base0(:,ix)=A_base(:,ix)
      enddo

               		! expedient for name clash
      call NAGdiag(klimit0)
      autoval0=eval; tridia0=hmat		! expedient for name clash

      phases0(:)=exp(-uc*htime*autoval0(:))

      c_of_t0=0.0
      do mm=1,klimit0
         do ilamb=1,klimit0
           c_of_t0(mm)=c_of_t0(mm)+tridia0(1,ilamb)*tridia0(mm,ilamb)*phases0(ilamb)
         enddo
      enddo

      x_of_t=matmul(A_base0,c_of_t0)

      deallocate(hmat,eval,tridia0,autoval0,phases0,c_of_t0,A_base0)

      return
      end
