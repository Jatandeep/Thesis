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
	double grid_scale,b,a,h,x;
	unsigned int gl_ref,lc_ref;

    static void declare_param(dealii::ParameterHandler& prm);
    void parse_param(dealii::ParameterHandler &prm);
    };

    class MaterialModel{
    public:
    double lambda,mu;

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
        double res_tol_u,res_tol_d,nu_tol_u;

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
        std::string time_adap;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };
	
    class PhaseFieldMethod{
    public:
        double g_c,l,k,u_total,st_tol,viscosity;
       
        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };

    class TestCase{
    public:
        std::string test;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };

    class BoundaryConditions{
    public:
        std::string uxb,uxt;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);
    };

    class ModelingStrategy{
    public:
        std::string strategy,comp_strategy;
        double fac_ft,steps_ft;

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
        TestCase test_case;
        BoundaryConditions bc;
        ModelingStrategy mod_strategy;

        static void declare_param(dealii::ParameterHandler& prm);
        void parse_param(dealii::ParameterHandler &prm);

    };
    }

}
