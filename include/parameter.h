#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria.h>

namespace thesis
{
    namespace parameters
    {

    class GeometryModel{
    public:
        std::string meshfile;
	double act_ref, act_cors,grid_scale;
	unsigned int gl_ref,lc_ref;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };

    class MaterialModel{
    public:
        double elastic_mod;
        double poisson_ratio;
        double lambda,mu;
	double viscosity;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };


    class FESys
    {
    public:
        unsigned int fe_degree;
        unsigned int quad_order;
        
        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };


    class NewtonRaphson{
    public:
        unsigned int max_new_ite;
        double res_tol_u,res_tol_d,nu_tol_u,nu_tol_d;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };


    class LinearSolver{
    public:
        double cg_tol,relax_prm;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };


    class Time{
    public:
        double start_time, end_time;
        double delta_t, time_tol,delta_t_f,time_change_interval;
	unsigned int op_freq;
        double alpha,beta;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };
	
    class PhaseFieldMethod{
    public:
        double g_c,l,k,u_total;
       
        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };
	
    class AllParameters{
    public:
        AllParameters(const std::string &filename);

        FESys fesys;
        GeometryModel geometrymodel;
        MaterialModel materialmodel;
	NewtonRaphson newtonraphson;
	LinearSolver linearsolver;
	Time time;	
	PhaseFieldMethod pf;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);

    };
    }

}
