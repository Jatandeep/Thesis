#pragma once

#include "../include/ElasticProblem.h"

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

