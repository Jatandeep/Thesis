#pragma once

#include "../include/ElasticProblem.h"

template <int dim>
dealii::SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda
                                                 	,const double mu);

template <int dim>
dealii::SymmetricTensor<4,dim> get_BigC(const double lambda
					,const double mu
					,dealii::SymmetricTensor<2,dim>& epsilon);

template <int dim>
dealii::SymmetricTensor<2,dim> get_stress(const double lambda
					 ,const double mu
					 ,dealii::SymmetricTensor<2,dim> &epsilon);

