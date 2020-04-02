#pragma once

#include "../include/ElasticProblem.h"

template<int dim>
void print_tensor(dealii::SymmetricTensor<4,dim>&C);

template<int dim>
void print_tensor(dealii::Tensor<4,dim>&C);

template<int dim>
void print_tensor(dealii::SymmetricTensor<2,dim>&C);

template<int dim>
void print_tensor(dealii::Tensor<2,dim>&C);
