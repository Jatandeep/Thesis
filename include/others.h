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
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

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
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

    virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double load_ratio_,u_total_;
	unsigned int itr_;
    };

 
    template <int dim>
    class InitialCrack:public dealii::Function<dim>{
    public:
        InitialCrack(const double min_cell_diameter):dealii::Function<dim>(dim+1)
	{
	_min_cell_diameter = min_cell_diameter;
	}
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

    virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double _min_cell_diameter;
    };

   
    template <int dim>
    class ElasticBodyForce:public dealii::Function<dim>{
    public:
        ElasticBodyForce():dealii::Function<dim>(1){}

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;

    };

    template <int dim>
    class Reference_solution : public dealii::Function<dim>
    {
    public:
      Reference_solution(): dealii::Function<dim>(dim+1)
      {}
      virtual double value(const dealii::Point<dim> & p,
                       const unsigned int component = 0) const override;
    };

    template <int dim>
    void comparison(const double lambda,const double mu,dealii::SymmetricTensor<2,dim> &dummy);

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

double get_stress_intensity_const(const std::string test);

std::pair<double,double> get_youngsM_poissonR(const double lambda,const double mu);

template<int dim>
double get_critical_load(dealii::Tensor<2,dim> stress);
