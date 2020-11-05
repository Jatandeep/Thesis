#ifndef CONSTITUTIVE_H
#define CONSTITUTIVE_H

#include "../include/utilities.h"


/*To get tensile/positive component of stress to be used in assembly of u*//*Thesis_report:Equation 2.38*/
template <int dim>
dealii::SymmetricTensor<2,dim> get_stress_plus(const double lambda
					      ,const double mu
					      ,dealii::SymmetricTensor<2,dim> &epsilon);

/*To get compressive/negative component of stress to be used in assembly of u*//*Thesis_report:Equation 2.39*/
template <int dim>
dealii::SymmetricTensor<2,dim> get_stress_minus(const double lambda
					      ,const double mu
					      ,dealii::SymmetricTensor<2,dim> &epsilon);

/*To get positive component of BigC to be used in assembly of u*//*Thesis_report:Equation 3.43*/
template <int dim>
dealii::SymmetricTensor<4,dim> get_BigC_plus(const double lambda
					    ,const double mu
					    ,dealii::SymmetricTensor<2,dim>& epsilon);

/*To get negative component of BigC to be used in assembly of u*//*Thesis_report:Equation 3.48*/
template <int dim>
dealii::SymmetricTensor<4,dim> get_BigC_minus(const double lambda
					     ,const double mu
					     ,dealii::SymmetricTensor<2,dim>& epsilon);

#endif
