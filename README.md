# cavity
$$H=\omega b^\dagger b + \sum_n\frac{\Delta_n(t)}{2}\sigma_n^z+g(t)\sum_n\sigma_n^x(b^\dagger+b)+g'(t)(b^\dagger+b)$$
## useful files
### getIndV.f90
contains subroutines to 
1. construct the H matrix indices & values
+   para2IndV(Nin,Ns,para): off-diagonal terms
+   para2V(Nin,Ns,para): diagonal terms
2. construct the initial state
+  init_coherent(Nin,Ns,alpha): |0000>|alpha>
+  init_S0coha(Nin,Ns,alpha): |S0a>|alpha>
+  init_S0ia(Nin,Ns,i_init): |S0a>|i_init>
   ...
3. calculate time-dependent expectation values and density matrix
+  getExT(Nin,Ns): calculate $\langle n_1\rangle$, $\langle n_2\rangle$, $\langle n_3\rangle$, $\langle n_4\rangle$, $\langle n_b\rangle$
+  getRhoT(Nin,Ns): calculate reduced density matrix (16x16)

### decay.f90
contains functions to switch on/off parameters or to send pulse

### lanc_pureS0_pulse.f90
The main file. 
1. Ns: number of two-level sites, fixed to be 4.
2. Nin: uplimit of photon numbers.
3. para: array of parameters. $\[\omega, \Delta_1, \Delta_2, \Delta_3, \Delta_4, g, g'\]$
