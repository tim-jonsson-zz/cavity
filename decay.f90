module decay
  implicit none
 contains
 real*8 function decay_factor(t,Gamma_out)
   real*8, intent(in) :: t,Gamma_out
   decay_factor = exp(-t*Gamma_out)
   return
 end function decay_factor
 
 real*8 function switchon_factor(t,T_sw)
   real*8, intent(in) :: t,T_sw
   real*8 :: PI
   PI=4.D0*DATAN(1.D0)
   if (t >= T_sw) then
     switchon_factor = 1.D0
   else  
     switchon_factor = 0.5D0 - cos(PI*t/T_sw)/2.D0
   end if  
   return
 end function switchon_factor
 
 real*8 function switchoff_factor(t,T_sw1,T_sw2,tget)
   real*8, intent(in) :: t,T_sw1,T_sw2,tget
   real*8 :: PI
   PI=4.D0*DATAN(1.D0)
   if (t <= T_sw1) then
     switchoff_factor = 1.D0
   else if (t <= T_sw2) then
     switchoff_factor = 0.5D0*(1+tget) + (1-tget)*cos(PI*(t-T_sw1)/(T_sw2-T_sw1))/2.D0
   else
     switchoff_factor = tget
   end if  
   return
 end function switchoff_factor
 
 real*8 function osci_factor(t,w0)
   real*8, intent(in) :: t,w0
   real*8 :: PI
   PI=4.D0*DATAN(1.D0)
   
     osci_factor = cos(w0*t)
   
   return
 end function osci_factor
 
 real*8 function superG_factor(t,tc,w0,tau,n)
   real*8, intent(in) :: t,tc,w0,tau
   integer*8, intent(in) :: n
   
   superG_factor = exp(-dlog(2.D0)/2.D0*(2.D0*(t-tc)/tau)**(2*n))*sin(w0*(t-tc))
   
   return
 end function superG_factor
 
 real*8 function superG_factor2(t,tc,w0,tau,n)
   real*8, intent(in) :: t,tc,w0,tau
   integer*8, intent(in) :: n
   
   superG_factor2 = exp(-dlog(2.D0)/2.D0*(2.D0*(t-tc)/tau)**(2*n))*cos(w0*(t-tc))
   
   return
 end function superG_factor2
 end module decay
