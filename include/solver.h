#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "const.h"
#include "domain.h"

using vector = Eigen::ArrayXd;
using matrix = Eigen::MatrixXd;

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

    vector rho_e_xn;
    vector rho_e_xp;
    vector rho_m_xn;
    vector rho_m_xp;

    vector rho_e_yn;
    vector rho_e_yp;
    vector rho_m_yn;
    vector rho_m_yp;

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
    vector sigma_e_xn;
    vector sigma_m_xn;
    vector sigma_e_xp;
    vector sigma_m_xp;

    double sigma_max_yn;
    double sigma_max_yp;
    vector sigma_e_yn;
    vector sigma_m_yn;
    vector sigma_e_yp;
    vector sigma_m_yp;

    // kappa
    vector kappa_e_xn;
    vector kappa_e_xp;
    vector kappa_m_xn;
    vector kappa_m_xp;

    vector kappa_e_yn;
    vector kappa_e_yp;
    vector kappa_m_yn;
    vector kappa_m_yp;

    // alpha
    vector alpha_e_xn;
    vector alpha_e_xp;
    vector alpha_m_xn;
    vector alpha_m_xp;

    vector alpha_e_yn;
    vector alpha_e_yp;
    vector alpha_m_yn;
    vector alpha_m_yp;

    // b
    vector b_e_xn_v;
    vector b_e_xp_v;
    vector b_m_xn_v;
    vector b_m_xp_v;

    vector b_e_yn_v;
    vector b_e_yp_v;
    vector b_m_yn_v;
    vector b_m_yp_v;

    // a
    vector a_e_xn_v;
    vector a_e_xp_v;
    vector a_m_xn_v;
    vector a_m_xp_v;

    vector a_e_yn_v;
    vector a_e_yp_v;
    vector a_m_yn_v;
    vector a_m_yp_v;

    // 场量更新公式系数
    double CD_dt;
    double CD_dt_dx;
    double CD_dt_dy;

    double CB_dt;
    double CB_dt_dx;
    double CB_dt_dy;

    vector C_Dx_dyn_v;
    vector C_Dx_dyp_v;

    matrix Phi_ex_yn;
    matrix b_ex_yn;
    matrix a_ex_yn;
    matrix C_Dx_dyn;

    matrix Phi_ex_yp;
    matrix b_ex_yp;
    matrix a_ex_yp;
    matrix C_Dx_dyp;

    vector temp_dx;

    vector C_Dy_dxn_v;
    vector C_Dy_dxp_v;

    matrix Phi_ey_xn;
    matrix b_ey_xn;
    matrix a_ey_xn;
    matrix C_Dy_dxn;

    matrix Phi_ey_xp;
    matrix b_ey_xp;
    matrix a_ey_xp;
    matrix C_Dy_dxp;

    vector temp_dy;

    vector C_Bz_dxn_v;
    vector C_Bz_dxp_v;

    vector C_Bz_dyn_v;
    vector C_Bz_dyp_v;

    matrix Phi_mz_xn;
    matrix b_mz_xn;
    matrix a_mz_xn;
    matrix C_Bz_dxn;

    matrix Phi_mz_xp;
    matrix b_mz_xp;
    matrix a_mz_xp;
    matrix C_Bz_dxp;

    vector temp_bz;

    matrix Phi_mz_yn;
    matrix b_mz_yn;
    matrix a_mz_yn;
    matrix C_Bz_dyn;

    matrix Phi_mz_yp;
    matrix b_mz_yp;
    matrix a_mz_yp;
    matrix C_Bz_dyp;

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

    void Allocate(Domain &domain);
    void CalculFdtdPara(Domain &domain);
    void CalculFdtdCoeff(Domain &domain);
};
