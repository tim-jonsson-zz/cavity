module fullED
contains
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
end module fullED
