#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria.h>


using namespace dealii;

struct Geometrymodel{
    std::string name;
    std::string meshfile;

    static void declare_param(ParameterHandler& prm);
    void parse_param(ParameterHandler &prm);
};

struct Materialmodel{
    double elastic_mod;
    double poisson_ratio;

    static void declare_param(ParameterHandler& prm);
    void parse_param(ParameterHandler &prm);
};


struct FESys
{
    unsigned int fe_degree;
    unsigned int quad_order;
    double gl_ref;
    unsigned int steps;
    double cycles;
    double tol,relax_prm;
    double lambda,mu;

    static void declare_param(ParameterHandler& prm);
    void parse_param(ParameterHandler &prm);
};


struct AllParameters:public FESys, public Geometrymodel, public Materialmodel{
    AllParameters(const std::string &filename);

    static void declare_param(ParameterHandler& prm);
    void parse_param(ParameterHandler &prm);

};


