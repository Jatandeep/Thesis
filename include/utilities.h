#pragma once
#include "../include/PhasefieldSMP.h"

/*Overloaded Functions for printing a tensor*/
template<int dim>
void print_tensor(const dealii::SymmetricTensor<4,dim>&C);

template<int dim>
void print_tensor(const dealii::Tensor<4,dim>&C);

template<int dim>
void print_tensor(const dealii::SymmetricTensor<2,dim>&C);

template<int dim>
void print_tensor(const dealii::Tensor<2,dim>&C);

template <int dim>
dealii::SymmetricTensor<2,dim> biaxial();

template <int dim>
dealii::SymmetricTensor<2,dim> uniaxial();
/*To define the definition of sign funciton*/
double get_sign(double x);

