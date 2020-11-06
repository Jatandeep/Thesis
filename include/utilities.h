#ifndef UTILITIES_H
#define UTILITIES_H

#include "../include/PhasefieldSMP.h"

/*Overloaded Functions for printing a tensor*/
/*Print 4th Order Symmetric Tensor*/
template<int dim>
void print_tensor(const dealii::SymmetricTensor<4,dim>&C);

/*Print 4th Order Tensor*/
template<int dim>
void print_tensor(const dealii::Tensor<4,dim>&C);

/*Print 2nd Order Symmetric Tensor*/
template<int dim>
void print_tensor(const dealii::SymmetricTensor<2,dim>&C);

/*Print 2nd Order Tensor*/
template<int dim>
void print_tensor(const dealii::Tensor<2,dim>&C);

/*Returns a biaxial strain*/
template <int dim>
dealii::SymmetricTensor<2,dim> biaxial();

/*Returns a uniaxial strain*/
template <int dim>
dealii::SymmetricTensor<2,dim> uniaxial();

/*To define the definition of sign funciton*//*To be used in constitutive.cpp in 
delsig_dellmbda_minus/plus and delsig_dellmbda_b_minus/plus*/
double get_sign(double x);

#endif