#pragma once
#include <iostream>
#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
//#include "../include/ElasticProblem.h"
#include "../include/Phasefield.h"
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

namespace thesis
{

    template <int dim>
    class BoundaryForce:public dealii::Function<dim>{
    public:
        BoundaryForce():dealii::Function<dim>(dim+1)
	{}
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;

    };

    template <int dim>
    class BoundaryTension:public dealii::Function<dim>{
    public:
        BoundaryTension(const double load_ratio):dealii::Function<dim>(dim+1)
	{
	load_ratio_ = load_ratio;
	}
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double load_ratio_;
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

}


template <int dim>
double get_energy_density_plus(const double lambda
			      ,const double mu
			      ,dealii::SymmetricTensor<2,dim> &eps);
template <int dim>
double get_energy_density_minus(const double lambda
			      ,const double mu
			      ,dealii::SymmetricTensor<2,dim> &eps);

