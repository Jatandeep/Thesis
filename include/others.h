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
        BoundaryForce();
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;

    };

    template <int dim>
    class BoundaryTension:public dealii::Function<dim>{
    public:
        BoundaryTension(const double time):dealii::Function<dim>(dim+1)
	{
	time_ = time;
	}
	
	virtual double value (const dealii::Point<dim> &p,
                    const unsigned int component = 0) const;

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;
    private:
	double time_;
    };

    template <int dim>
    class ElasticBodyForce:public dealii::Function<dim>{
    public:
        ElasticBodyForce():dealii::Function<dim>(1){}

        virtual void vector_value (const dealii::Point<dim> &p,
                             dealii::Vector<double>  &value)const;

    };
}

/*
template <int dim>
double get_history(const double lambda
		,const double mu
		,dealii::BlockVector<double> &solution);
*/
