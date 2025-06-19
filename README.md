# cavity

## useful files
### getIndV.f90
contains subroutines to 
1. construct the basis and the H matrix indices & values
+   para2IndV(Nin,Ns,para): off-diagonal terms
+   para2V(Nin,Ns,para): diagonal terms
2. construct the initial state
+  init_coherent(Nin,Ns,alpha): |0000>|alpha>
+  init_S0coha(Nin,Ns,alpha): |S0a>|alpha>
+  init_S0ia(Nin,Ns,i_init): |S0a>|i_init>
   ...
3. calculate time-dependent expectation values and density matrix
+  getExT(Nin,Ns): calculate $\langle n1\rangle$,<n2>,<n3>,<n4>,<nb>
+  getRhoT(Nin,Ns): calculate reduced density matrix (16x16)
