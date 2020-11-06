#ifndef OTHERS_H
#define OTHERS_H

#include <iostream>
#include "../include/utilities.h"
#include "../include/constitutive.h"
namespace thesis
{
    /*A dealii type function for Tension test for prescribing load on top boundary*//*Thesis_report:Section4.1.1*/
    /*It takes 3 arguments. Number of newton iterations(itr_)
      as input as inhomogeneous bc needs to be applied in only first iteration. 
      "load_ratio_" decides the per step displacement load to be applied. "u_total_" is the total
      displacement that will be applied during the whole simulation.*/
    template <int dim>
    class BoundaryTension:public dealii::Function<dim>{
    public:
        BoundaryTension(unsigned int itr,const double load_ratio,const double u_total):dealii::Function<dim>(dim+1)
	{
	itr_ = itr;
	load_ratio_ = load_ratio;
	u_total_ = u_total;	
	}

    virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double load_ratio_,u_total_;
	unsigned int itr_;
    };

    /*A dealii type function for shear test for prescribing load on top boundary*//*Thesis_report:Section4.1.2*/
    /*It takes 3 arguments. Number of newton iterations(itr_)
      as input as inhomogeneous bc needs to be applied in only first iteration. 
      "load_ratio_" decides the per step displacement load to be applied. "u_total_" is the total
      displacement that will be applied during the whole simulation.*/
    template <int dim>
    class BoundaryShear:public dealii::Function<dim>{
    public:
        BoundaryShear(unsigned int itr,const double load_ratio,const double u_total):dealii::Function<dim>(dim+1)
	{
	itr_ = itr;
	load_ratio_ = load_ratio;
	u_total_ = u_total;
	}
	
    virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double load_ratio_,u_total_;
	unsigned int itr_;
    };

    /*For M_Id, this dealii type function assigns d=1 on crack boundary*//*Thesis_report:Section4.3.2*/
    /*Newton Iteration(itr_) required as we only need to provide Dirichlet condition
      of d=1 in first iteration.*/
    template <int dim>
    class InitialCrack:public dealii::Function<dim>{
    public:
        InitialCrack(unsigned int itr):dealii::Function<dim>(dim+1)
	{
    itr_ = itr;
	}
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

    virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
    unsigned int itr_;
    };

    /*For LEFM, this dealii type function prescribe the loading conditions on boundaries*//*Thesis_report:Section4.4*/
    /*Number of newton iterations(itr_) as input as inhomogeneous bc needs to be applied in only first iteration. "steps_" are required
      to get an per step increment of K_I."g_c,lamba,mu" required for K_IC calculation.*/
    template <int dim>
    class Reference_solution:public dealii::Function<dim>{
    public:
        Reference_solution(unsigned int itr,const double steps,const double g_c
                          ,const double lambda
                          ,const double mu):dealii::Function<dim>(dim+1)
	{
	itr_ = itr;
	steps_ = steps;
    g_c_ = g_c;
	lambda_ = lambda;
    mu_ = mu;
	}
	
    virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double steps_,lambda_,mu_,g_c_;
	unsigned int itr_;
    };    

}

/*Generate energy density plus values for calculating total Elastic energy(for statistics file)*//*Thesis_report:Equation 2.32*/ 
/*"lambda,mu" are lame constant. eps is strain. All are required to calculate +ve energy density.
Return type is scalar +ve elastic energy value*/
template <int dim>
double get_energy_density_plus(const double lambda
			      ,const double mu
			      ,dealii::SymmetricTensor<2,dim> &eps);
/*Generate energy density minus values for calculating total Elastic energy(for statistics file)*//*Thesis_report:Equation 2.33*/
/*"lambda,mu" are lame constant. eps is strain. All are required to calculate -ve energy density.
Return type is scalar -ve elastic energy value*/
template <int dim>
double get_energy_density_minus(const double lambda
			      ,const double mu
			      ,dealii::SymmetricTensor<2,dim> &eps);

/*Define the degradation function to be used in formulation*//*Thesis_report:Equation 2.36*/
/*"d" is phase field value. Return type is value of degradation function.*/
double get_deg_func(const double d);

/*Generate Youngs modulus and poisson ratio as a std::pair by taking lame parameters"lambda,mu".*/
std::pair<double,double> get_youngsM_poissonR(const double lambda,const double mu);

#endif