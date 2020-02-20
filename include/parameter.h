#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria.h>

namespace thesis
{

    class Geometrymodel{
    public:
        std::string name;
        std::string meshfile;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };

    class Materialmodel{
    public:
        double elastic_mod;
        double poisson_ratio;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };


    class FESys
    {
    public:
        unsigned int fe_degree;
        unsigned int quad_order;
        double gl_ref;
        unsigned int steps;
        double cycles;
        double tol,relax_prm;
        double lambda,mu;
        double act_ref, act_cors;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);

        friend class AllParameters;

    };


    class AllParameters:public FESys, public Geometrymodel, public Materialmodel{
    public:
        AllParameters(const std::string &filename);

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);

    };



}
