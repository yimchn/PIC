#pragma once

#include <spdlog/spdlog.h>

#include <Eigen/Core>
#include <iostream>

#include "const.h"
#include "domain.h"

using vector = Eigen::ArrayXd;
using matrix = Eigen::ArrayXXd;


struct Solver {
    Domain &dm;
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

    vector rho_e_xn = vector::Zero(dm.geo.n_pml_xn);
    vector rho_e_xp = vector::Zero(dm.geo.n_pml_xp);
    vector rho_m_xn = vector::Zero(dm.geo.n_pml_xn);
    vector rho_m_xp = vector::Zero(dm.geo.n_pml_xp);

    vector rho_e_yn = vector::Zero(dm.geo.n_pml_yn);
    vector rho_e_yp = vector::Zero(dm.geo.n_pml_yp);
    vector rho_m_yn = vector::Zero(dm.geo.n_pml_yn);
    vector rho_m_yp = vector::Zero(dm.geo.n_pml_yp);

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
    vector sigma_e_xn = vector::Zero(dm.geo.n_pml_xn);
    vector sigma_m_xn = vector::Zero(dm.geo.n_pml_xn);
    vector sigma_e_xp = vector::Zero(dm.geo.n_pml_xp);
    vector sigma_m_xp = vector::Zero(dm.geo.n_pml_xp);

    double sigma_max_yn;
    double sigma_max_yp;
    vector sigma_e_yn = vector::Zero(dm.geo.n_pml_yn);
    vector sigma_m_yn = vector::Zero(dm.geo.n_pml_yn);
    vector sigma_e_yp = vector::Zero(dm.geo.n_pml_yp);
    vector sigma_m_yp = vector::Zero(dm.geo.n_pml_yp);

    // kappa
    vector kappa_e_xn = vector::Zero(dm.geo.n_pml_xn);
    vector kappa_e_xp = vector::Zero(dm.geo.n_pml_xp);
    vector kappa_m_xn = vector::Zero(dm.geo.n_pml_xn);
    vector kappa_m_xp = vector::Zero(dm.geo.n_pml_xp);

    vector kappa_e_yn = vector::Zero(dm.geo.n_pml_yn);
    vector kappa_e_yp = vector::Zero(dm.geo.n_pml_yp);
    vector kappa_m_yn = vector::Zero(dm.geo.n_pml_yn);
    vector kappa_m_yp = vector::Zero(dm.geo.n_pml_yp);

    // alpha
    vector alpha_e_xn = vector::Zero(dm.geo.n_pml_xn);
    vector alpha_e_xp = vector::Zero(dm.geo.n_pml_xp);
    vector alpha_m_xn = vector::Zero(dm.geo.n_pml_xn);
    vector alpha_m_xp = vector::Zero(dm.geo.n_pml_xp);

    vector alpha_e_yn = vector::Zero(dm.geo.n_pml_yn);
    vector alpha_e_yp = vector::Zero(dm.geo.n_pml_yp);
    vector alpha_m_yn = vector::Zero(dm.geo.n_pml_yn);
    vector alpha_m_yp = vector::Zero(dm.geo.n_pml_yp);

    // b
    vector b_e_xn_v = vector::Zero(dm.geo.n_pml_xn);
    vector b_e_xp_v = vector::Zero(dm.geo.n_pml_xp);
    vector b_m_xn_v = vector::Zero(dm.geo.n_pml_xn);
    vector b_m_xp_v = vector::Zero(dm.geo.n_pml_xp);

    vector b_e_yn_v = vector::Zero(dm.geo.n_pml_yn);
    vector b_e_yp_v = vector::Zero(dm.geo.n_pml_yp);
    vector b_m_yn_v = vector::Zero(dm.geo.n_pml_yn);
    vector b_m_yp_v = vector::Zero(dm.geo.n_pml_yp);

    // a
    vector a_e_xn_v = vector::Zero(dm.geo.n_pml_xn);
    vector a_e_xp_v = vector::Zero(dm.geo.n_pml_xp);
    vector a_m_xn_v = vector::Zero(dm.geo.n_pml_xn);
    vector a_m_xp_v = vector::Zero(dm.geo.n_pml_xp);

    vector a_e_yn_v = vector::Zero(dm.geo.n_pml_yn);
    vector a_e_yp_v = vector::Zero(dm.geo.n_pml_yp);
    vector a_m_yn_v = vector::Zero(dm.geo.n_pml_yn);
    vector a_m_yp_v = vector::Zero(dm.geo.n_pml_yp);

    // 场量更新公式系数
    double CD_dt;
    double CD_dt_dx;
    double CD_dt_dy;

    double CB_dt;
    double CB_dt_dx;
    double CB_dt_dy;

    // Dx
    vector C_Dx_dyn_v = vector::Zero(dm.geo.n_pml_yn);
    vector C_Dx_dyp_v = vector::Zero(dm.geo.n_pml_yp);

    matrix Phi_ex_yn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);
    matrix b_ex_yn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);
    matrix a_ex_yn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);
    matrix C_Dx_dyn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);

    matrix Phi_ex_yp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yp);
    matrix b_ex_yp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yp);
    matrix a_ex_yp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yp);
    matrix C_Dx_dyp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yp);

    // Dy
    vector C_Dy_dxn_v = vector::Zero(dm.geo.n_pml_xn);
    vector C_Dy_dxp_v = vector::Zero(dm.geo.n_pml_xp);

    matrix Phi_ey_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);
    matrix b_ey_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);
    matrix a_ey_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);
    matrix C_Dy_dxn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);

    matrix Phi_ey_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);
    matrix b_ey_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);
    matrix a_ey_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);
    matrix C_Dy_dxp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);

    // Dz
    vector C_Dz_dxn_v = vector::Zero(dm.geo.n_pml_xn);
    vector C_Dz_dxp_v = vector::Zero(dm.geo.n_pml_xp);

    vector C_Dz_dyn_v = vector::Zero(dm.geo.n_pml_yn);
    vector C_Dz_dyp_v = vector::Zero(dm.geo.n_pml_yp);

    matrix Phi_ez_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);
    matrix b_ez_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);
    matrix a_ez_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);
    matrix C_Dz_dxn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);

    matrix Phi_ez_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);
    matrix b_ez_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);
    matrix a_ez_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);
    matrix C_Dz_dxp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);

    matrix Phi_ez_yn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);
    matrix b_ez_yn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);
    matrix a_ez_yn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);
    matrix C_Dz_dyn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);

    matrix Phi_ez_yp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);
    matrix b_ez_yp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);
    matrix a_ez_yp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);
    matrix C_Dz_dyp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);

    // Bx
    vector C_Bx_dyn_v = vector::Zero(dm.geo.n_pml_yn);
    vector C_Bx_dyp_v = vector::Zero(dm.geo.n_pml_yp);

    matrix Phi_mx_yn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);
    matrix b_mx_yn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);
    matrix a_mx_yn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);
    matrix C_Bx_dyn = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yn);

    matrix Phi_mx_yp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);
    matrix b_mx_yp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);
    matrix a_mx_yp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);
    matrix C_Bx_dyp = matrix::Zero(dm.geo.nj, dm.geo.n_pml_yp);

    // By
    vector C_By_dxn_v = vector::Zero(dm.geo.n_pml_xn);
    vector C_By_dxp_v = vector::Zero(dm.geo.n_pml_xp);

    matrix Phi_my_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);
    matrix b_my_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);
    matrix a_my_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);
    matrix C_By_dxn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.ni);

    matrix Phi_my_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);
    matrix b_my_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);
    matrix a_my_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);
    matrix C_By_dxp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.ni);

    // Bz
    vector C_Bz_dxn_v = vector::Zero(dm.geo.n_pml_xn);
    vector C_Bz_dxp_v = vector::Zero(dm.geo.n_pml_xp);

    vector C_Bz_dyn_v = vector::Zero(dm.geo.n_pml_yn);
    vector C_Bz_dyp_v = vector::Zero(dm.geo.n_pml_yp);

    matrix Phi_mz_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);
    matrix b_mz_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);
    matrix a_mz_xn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);
    matrix C_Bz_dxn = matrix::Zero(dm.geo.n_pml_xn, dm.geo.nj);

    matrix Phi_mz_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);
    matrix b_mz_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);
    matrix a_mz_xp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);
    matrix C_Bz_dxp = matrix::Zero(dm.geo.n_pml_xp, dm.geo.nj);

    matrix Phi_mz_yn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);
    matrix b_mz_yn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);
    matrix a_mz_yn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);
    matrix C_Bz_dyn = matrix::Zero(dm.geo.ni, dm.geo.n_pml_yn);

    matrix Phi_mz_yp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_xn);
    matrix b_mz_yp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_xn);
    matrix a_mz_yp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_xn);
    matrix C_Bz_dyp = matrix::Zero(dm.geo.ni, dm.geo.n_pml_xn);

    /*constructor, sets world*/
    Solver(Domain &domain, int max_it, double tol);

    /*solves potential using Gauss-Seidel*/
    bool solve();

    /*computes electric field = -gradient(phi)*/
    void computeEF();

    void InitFdtdPara();
    void CalculFdtdCoeff();
    // update the magnetic and electric field under 2d using FDTD method
    void StepForward();
    matrix &UpdateSource();
    matrix &UpdateDz();
    matrix &UpdateBx();
    matrix &UpdateBy();
};