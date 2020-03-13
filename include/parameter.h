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
        unsigned int gl_ref;
        unsigned int cg_steps;
        double cg_tol,relax_prm;
        double lambda,mu;
        double act_ref, act_cors;
        unsigned int start_time, end_time;
        double delta_t, time_tol;
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
