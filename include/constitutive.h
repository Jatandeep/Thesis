#pragma once

#include "../include/ElasticProblem.h"

template <int dim>
dealii::SymmetricTensor<4, dim> get_const_BigC(const double lambda
                                                 	,const double mu);

template <int dim>
dealii::SymmetricTensor<4,dim> get_BigC(const double lambda
					,const double mu
					,dealii::SymmetricTensor<2,dim>& epsilon);

template <int dim>
dealii::SymmetricTensor<2,dim> get_stress(const double lambda
					 ,const double mu
					 ,dealii::SymmetricTensor<2,dim> &epsilon);


template <int dim>
void print_C(dealii::SymmetricTensor<4,dim> &C);

template <int dim>
void print_eps(dealii::SymmetricTensor<2,dim> &eps);

template <int dim>
dealii::SymmetricTensor<2,dim> biaxial();

template <int dim>
dealii::SymmetricTensor<2,dim> uniaxial();
