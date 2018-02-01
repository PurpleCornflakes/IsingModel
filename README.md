Plese compile using c++11 standard.

```
g++ -o ss_gen -std=c++11 ss_generator.cpp
```

If the g++ is not g++7 version(e.g HPC and NSCC), just use
```
gcc -o eq_test -lstdc++ -std=c++0x equibm_test.cpp
```
Note '-lstdc++' is to link c++ libs.

To give Ising Lattice Configuration,
```
./ss_gen 2.27 fBC 5 5 5
```
T = 2.27 (which is Tc)
fBC/pBC: free/periodic Boundary Condition
5 5 5: valid spin lattice is of dimension (5,5,5)


#### E_T, C_T
32X32 spins
5000 iterations
![](imgs/e_T.png)
![](imgs/ee_T.png)
![](imgs/c_T.png)

#### M_T, chi_T
32X32 spins
5000 iterations
![](imgs/m_T.png)
![](imgs/chi_T.png)

#### M distribution
100X100 spins
5000 iterations
![](imgs/M_distribution.png)

