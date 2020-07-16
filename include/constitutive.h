#pragma once
#include "../include/utilities.h"


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
dealii::SymmetricTensor<2,dim> get_stress_plus(const double lambda
					      ,const double mu
					      ,dealii::SymmetricTensor<2,dim> &epsilon);
template <int dim>
dealii::SymmetricTensor<2,dim> get_stress_minus(const double lambda
					      ,const double mu
					      ,dealii::SymmetricTensor<2,dim> &epsilon);


template <int dim>
dealii::SymmetricTensor<4,dim> get_BigC_plus(const double lambda
					    ,const double mu
					    ,dealii::SymmetricTensor<2,dim>& epsilon);


template <int dim>
dealii::SymmetricTensor<4,dim> get_BigC_minus(const double lambda
					     ,const double mu
					     ,dealii::SymmetricTensor<2,dim>& epsilon);


