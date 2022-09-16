#pragma once

#include <array>
#include <memory>

#include "const.h"
#include "domain.h"

#define N_PML_XN 20;
#define N_PML_XP 20;
#define N_PML_YN 20;
#define N_PML_YP 20;

struct Solver {
    Domain &domain;
    unsigned max_solver_it;  // maximum number of solver iterations
    double tolerance;        // solver tolerance

    // FDTD parameter
    double m_pml = 3;
    double alpha_max = 0.5;
    double sigma_factor = 1;
    double kappa_max = 7;

    double sigma_factor_xn = sigma_factor;
    double kappa_max_xn = kappa_max;
    double sigma_factor_xp = sigma_factor;
    double kappa_max_xp = kappa_max;
    double sigma_factor_yn = sigma_factor;
    double kappa_max_yn = kappa_max;
    double sigma_factor_yp = sigma_factor;
    double kappa_max_yp = kappa_max;

    // remember to delete!!!!!
    double *rho_e_xn = new double[domain.geo.n_pml_xn];
    double *rho_e_xp = new double[domain.geo.n_pml_xp];
    double *rho_m_xn = new double[domain.geo.n_pml_xn];
    double *rho_m_xp = new double[domain.geo.n_pml_xp];

    double *rho_e_yn = new double[domain.geo.n_pml_yn];
    double *rho_e_yp = new double[domain.geo.n_pml_yp];
    double *rho_m_yn = new double[domain.geo.n_pml_yn];
    double *rho_m_yp = new double[domain.geo.n_pml_yp];

    double eps_r_pml_xn = 1;
    double mu_r_pml_xn = 1;
    double eps_r_pml_xp = 1;
    double mu_r_pml_xp = 1;
    double eps_r_pml_yn = 1;
    double mu_r_pml_yn = 1;
    double eps_r_pml_yp = 1;
    double mu_r_pml_yp = 1;

    // sigma
    double sigma_max_xn;
    double sigma_max_xp;
    double *sigma_e_xn = new double[domain.geo.n_pml_xn];
    double *sigma_m_xn = new double[domain.geo.n_pml_xn];
    double *sigma_e_xp = new double[domain.geo.n_pml_xp];
    double *sigma_m_xp = new double[domain.geo.n_pml_xp];

    double sigma_max_yn;
    double sigma_max_yp;
    double *sigma_e_yn = new double[domain.geo.n_pml_yn];
    double *sigma_m_yn = new double[domain.geo.n_pml_yn];
    double *sigma_e_yp = new double[domain.geo.n_pml_yp];
    double *sigma_m_yp = new double[domain.geo.n_pml_yp];

    // kappa
    double *kappa_e_xn = new double[domain.geo.n_pml_xn];
    double *kappa_e_xp = new double[domain.geo.n_pml_xp];
    double *kappa_m_xn = new double[domain.geo.n_pml_xn];
    double *kappa_m_xp = new double[domain.geo.n_pml_xp];

    double *kappa_e_yn = new double[domain.geo.n_pml_yn];
    double *kappa_e_yp = new double[domain.geo.n_pml_yp];
    double *kappa_m_yn = new double[domain.geo.n_pml_yn];
    double *kappa_m_yp = new double[domain.geo.n_pml_yp];

    // alpha
    double *alpha_e_xn = new double[domain.geo.n_pml_xn];
    double *alpha_e_xp = new double[domain.geo.n_pml_xp];
    double *alpha_m_xn = new double[domain.geo.n_pml_xn];
    double *alpha_m_xp = new double[domain.geo.n_pml_xp];

    double *alpha_e_yn = new double[domain.geo.n_pml_yn];
    double *alpha_e_yp = new double[domain.geo.n_pml_yp];
    double *alpha_m_yn = new double[domain.geo.n_pml_yn];
    double *alpha_m_yp = new double[domain.geo.n_pml_yp];

    // b
    double *b_e_xn_v = new double[domain.geo.n_pml_xn];
    double *b_e_xp_v = new double[domain.geo.n_pml_xp];
    double *b_m_xn_v = new double[domain.geo.n_pml_xn];
    double *b_m_xp_v = new double[domain.geo.n_pml_xp];

    double *b_e_yn_v = new double[domain.geo.n_pml_yn];
    double *b_e_yp_v = new double[domain.geo.n_pml_yp];
    double *b_m_yn_v = new double[domain.geo.n_pml_yn];
    double *b_m_yp_v = new double[domain.geo.n_pml_yp];

    // a
    double *a_e_xn_v = new double[domain.geo.n_pml_xn];
    double *a_e_xp_v = new double[domain.geo.n_pml_xp];
    double *a_m_xn_v = new double[domain.geo.n_pml_xn];
    double *a_m_xp_v = new double[domain.geo.n_pml_xp];

    double *a_e_yn_v = new double[domain.geo.n_pml_yn];
    double *a_e_yp_v = new double[domain.geo.n_pml_yp];
    double *a_m_yn_v = new double[domain.geo.n_pml_yn];
    double *a_m_yp_v = new double[domain.geo.n_pml_yp];

    /*constructor, sets world*/
    Solver(Domain &domain, int max_it, double tol);

    /*solves potential using Gauss-Seidel*/
    bool solve();

    /*computes electric field = -gradient(phi)*/
    void computeEF();

    void UpdateHz2d();
    void UpdateHx2d();
    void UpdateHy2d();
    void UpdateEz2d();
    void UpdateEx2d();
    void UpdateEy2d();

    // update the magnetic and electric field under 2d using FDTD method
    void UpdateElectromagnetic();
};
