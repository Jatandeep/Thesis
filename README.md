# Thesis

This is a finite element code based in deal.II to simulate crack propagation in brittle materials using phase field model. 
The above project has been developed for a master thesis in FAU Erlangen-Nuremebrg.
Please refer to Algorithm-1 in thesis-report in order to understand the project completely.

# Mesh file
The mesh name nomenclature used in parameter files name is described below:
1. **mesh01** -> This mesh is used for M_I and M_Id type crack in tension test and contains pre-refined mesh where crack is expected to grow with element size h = 0.001 mm (approx).(Tension_left_1.inp)
2. **mesh02** -> This mesh is used for P_I type crack for tension test and contains pre-refined mesh where crack is to be prescribed as well as where the crack is expected to grow with element size h = 0.001 mm (approx).(Tension_left_10.inp)
3. **mesh03** -> This mesh is used for M_I type crack for shear test and contains pre-refined mesh where crack is expected to grow with element size h = 0.002 mm (approx).(Shear_left_2.inp)
4. **mesh04** -> This mesh is used for M_I and M_Id type crack for lefm scenarios and subsequent global and local refinement will be implemented in the program giving us a element size of h = 0.001 mm(approx).(Tension_left_lefm_1.inp)
5. **mesh05** -> This mesh is used for P_I type crack for lefm scenarios and contains a pre-refined mesh where crack is to be prescribed and expected to grow with element size of h = 0.001 mm(approx).(Tension_left_lefm_2.inp)
6. **mesh06_01/02** -> This mesh is used for P_I type crack for tension test using method of single row of nodes as a crack  and is refined appropriately with local and global refinement to get element sizes of h = 0.001 mm(**mesh06-01**)and h = 0.0007 mm (**mesh06-02**) respectively.(Tension_left_7.inp)  

# Parameter file
In order to replicate the results in thesis-report, various parameter files have been added to folder named "parameter_files". 
We will take an example and explain the meaning of that parameter file. 
"**M_I_Tension_l-0.015_v-0_mesh01_uxbfxd-uxtfree.prm**":
1. **M_I** -> It denotes the pre-existing crack modeling strategy. (**M_I** | **M_Id** | **P_I**)
2. **Tension** -> It shows the type of test begin performed on the specimen. (**Tension** | **Shear** | **lefm**)
3. **l-0.015** -> Regularization length (**0.015** | **0.0075**)
4. **v-0** -> Viscosity (**0** | **1e-6** | **1e-5** | **0.5e-4**)
5. **mesh01** -> Corresponding mesh file name.
6. **uxbfxd-uxtfree** -> **Optional** parameter showing boundary conditions for tension test. If not mentioned, **uxbfxd-uxtfree** is assumed. (**uxbfxd-uxtfxd** | **uxbfxd-uxtfree** | **uxbfree-uxtfree**)

- You can find parameter files relation with thesis report figures [here](Figures.md)
- Detailed explanation of variables to be given in parameter file is given [here](Parameter.md) 


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
  ./crack /path/to/parameter/file
```

# Helpful material
For getting an overview of block solvers, multithreading, input parameter handling and linear newton system from non-linear equations, 
please refer to Step-44 of deal.II tutorials (https://dealii.org/9.0.0/doxygen/deal.II/step_44.html).

# Notes
The code uses the deal.II 9.0.0 version. If proejct is run with newer versions, certain version specific modifications
might need to be done.
