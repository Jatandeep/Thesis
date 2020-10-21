#pragma once
#include <iostream>
#include "../include/utilities.h"
#include "../include/constitutive.h"
namespace thesis
{

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


template <int dim>
double get_energy_density_plus(const double lambda
			      ,const double mu
			      ,dealii::SymmetricTensor<2,dim> &eps);
template <int dim>
double get_energy_density_minus(const double lambda
			      ,const double mu
			      ,dealii::SymmetricTensor<2,dim> &eps);

double get_deg_func(const double d);

std::pair<double,double> get_youngsM_poissonR(const double lambda,const double mu);

template<int dim>
double get_critical_load(dealii::Tensor<2,dim> stress);
