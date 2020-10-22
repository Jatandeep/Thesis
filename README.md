# Thesis

This is a finite element code based in deal.II to simulate crack propagation in brittle materials using phase field model. 
The above project has been developed for a master thesis in FAU Erlangen-Nuremebrg.
Please refer to Algorithm-1 in thesis-report in order to understand the project completely.

# Parameter file
In order to replicate the results in thesis-report, various parameter files have been added to folder named "parameter_files". 
We will take an example and explain the meaning of that parameter file. 
"**M_I_Tension_l=0.015_v=0_Fig4.6_uxbfxd-uxtfree.prm**":
1. **M_I** -> It denotes the pre-existing crack modeling strategy. (**M_I** | **M_Id** | **P_I**)
2. **Tension** -> It shows the type of test begin performed on the specimen. (**Tension** | **Shear** | **lefm**)
3. **l=0.015** -> Regularization length (**0.015** | **0.0075**)
4. **v=0** -> Viscosity (**0** | **1e-6** | **1e-5** | **0.5e-4**)
5. **Fig4.6** -> Corresponding figure number in master thesis report.
6. **uxbfxd-uxtfree** -> **Optional** parameter showing boundary conditions for tension test. If not mentioned, **uxbfxd-uxtfree** is assumed. (**uxbfxd-uxtfxd** | **uxbfxd-uxtfree** | **uxbfree-uxtfree**)

# How to run

You need to install deal.II (see http://www.dealii.org) from the official site or using spack (https://github.com/dealii/dealii/wiki/deal.II-in-Spack). 
Download the code and configure with:
```
  spack load dealii
```
```
  cmake .
```
Compile with:
```
  make
```
and finally run with:
```
  ./crack parameter.prm
```

# Helpful material
For getting an overview of block solvers, multithreading, input parameter handling and linear newton system from non-linear equations, 
please refer to Step-44 of deal.II tutorials (https://dealii.org/9.0.0/doxygen/deal.II/step_44.html).

# Notes
The code uses the deal.II 9.0.0 version. If proejct is run with newer versions, certain version specific modifications
might need to be done.
