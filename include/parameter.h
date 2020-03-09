#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria.h>

namespace thesis
{
    namespace parameters
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
        double tol,relax_prm;
        double lambda,mu;
        double act_ref, act_cors;
        unsigned int n_time_steps;
        unsigned int max_new_ite;
        double res_tol;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);



    };


    class AllParameters{
    public:
        AllParameters(const std::string &filename);

        FESys fesys;
        Geometrymodel geometrymodel;
        Materialmodel materialmodel;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);

    };
    }

}
