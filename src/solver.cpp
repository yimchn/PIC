#include "solver.h"

#include <math.h>

#include <iostream>

#include "const.h"
#include "field.h"

Solver::Solver(Domain &domain, int max_it, double tol)
    : domain(domain), max_solver_it(max_it), tolerance(tol) {
    Allocate(domain);
    CalculFdtdPara(domain);
    CalculFdtdCoeff(domain);
}
void Solver::CalculFdtdCoeff(Domain &domain) {
    CD_dt = domain.getDt();
    CD_dt_dx = domain.getDt() / domain.geo.dh[0];
    CD_dt_dy = domain.getDt() / domain.geo.dh[1];

    CB_dt = domain.getDt();
    CB_dt_dx = domain.getDt() / domain.geo.dh[0];
    CB_dt_dy = domain.getDt() / domain.geo.dh[1];

    C_Dx_dyn_v = CD_dt_dy / kappa_e_yn;
    C_Dx_dyp_v = CD_dt_dy / kappa_max_yp;

    // PML区域内Dx的更新公式系数
    temp_dx.setConstant(1);
    for (int i = 0; i < domain.geo.n_pml_yn; ++i) {
        b_ex_yn.col(i) = b_e_yn_v(i) * temp_dx;
        a_ex_yn.col(i) = a_e_yn_v(i) * temp_dx;
        C_Dx_dyn.col(i) = C_Dx_dyn_v(i) * temp_dx;
    }
    for (int i = 0; i < domain.geo.n_pml_yp; ++i) {
        b_ex_yp.col(i) = b_e_yp_v(i) * temp_dx;
        a_ex_yp.col(i) = a_e_yp_v(i) * temp_dx;
        C_Dx_dyp.col(i) = C_Dx_dyp_v(i) * temp_dx;
    }

    // PML区域内Dy的公式更新系数
    C_Dy_dxn_v = CD_dt_dx / kappa_e_xn;
    C_Dy_dxp_v = CD_dt_dx / kappa_e_xp;

    temp_dy.setConstant(1);
    for (int i = 0; i < domain.geo.n_pml_xn; ++i) {
        b_ey_xn.row(i) = b_e_xn_v(i) * temp_dy;
        a_ey_xn.row(i) = a_e_xn_v(i) * temp_dy;
        C_Dy_dxn.row(i) = C_Dy_dxn_v(i) * temp_dy;
    }
    for (int i = 0; i < domain.geo.n_pml_xp; ++i) {
        b_ey_xp.row(i) = b_e_xp_v(i) * temp_dy;
        a_ey_xp.row(i) = a_e_xp_v(i) * temp_dy;
        C_Dy_dxp.row(i) = C_Dy_dxp_v(i) * temp_dy;
    }

    // PML区域内By的更新公式系数
    temp_bz.setConstant(1);
    for (int i = 0; i < domain.geo.n_pml_xn; ++i) {
        b_mz_xn.row(i) = b_m_xn_v(i) * temp_bz;
        a_mz_xn.row(i) = a_m_xn_v(i) * temp_bz;
        C_Bz_dxn.row(i) = C_Bz_dxn_v(i) * temp_bz;
    }
    for (int i = 0; i < domain.geo.n_pml_xp; ++i) {
        b_mz_xp.row(i) = b_m_xp_v(i) * temp_bz;
        a_mz_xp.row(i) = a_m_xp_v(i) * temp_bz;
        C_Bz_dxp.row(i) = C_Bz_dxp_v(i) * temp_bz;
    }

    temp_bz.resize(domain.geo.ni);
    temp_bz.setConstant(1);
    for (int i = 0; i < domain.geo.n_pml_yn; ++i) {
        b_mz_yn.col(i) = b_m_yn_v(i) * temp_bz;
        a_mz_yn.col(i) = a_m_yn_v(i) * temp_bz;
        C_Bz_dyn.col(i) = C_Bz_dyn_v(i) * temp_bz;
    }
    for (int i = 0; i < domain.geo.n_pml_xp; ++i) {
        b_mz_xp.col(i) = b_m_yp_v(i) * temp_bz;
        a_mz_xp.col(i) = a_m_yp_v(i) * temp_bz;
        C_Bz_dxp.col(i) = C_Bz_dyp_v(i) * temp_bz;
    }
}

void Solver::CalculFdtdPara(Domain &domain) {
    // 计算电磁场量节点离开PML内边界的距离
    // x轴方向
    for (double i = domain.geo.n_pml_xn - 0.75; i > 0.25; i--)
        rho_e_xn[static_cast<int>(i)] = i / domain.geo.n_pml_xn;

    for (double i = 0.25; i < domain.geo.n_pml_xp - 0.75; i++)
        rho_e_xp[static_cast<int>(i)] = i / domain.geo.n_pml_xp;

    for (double i = domain.geo.n_pml_xn - 0.25; i > 0.75; i--)
        rho_m_xn[static_cast<int>(i)] = i / domain.geo.n_pml_xn;

    for (double i = 0.75; i < domain.geo.n_pml_xp - 0.25; i++)
        rho_m_xp[static_cast<int>(i)] = i / domain.geo.n_pml_xp;

    // y轴方向
    for (double i = domain.geo.n_pml_yn - 0.75; i > 0.25; i--)
        rho_e_yn[static_cast<int>(i)] = i / domain.geo.n_pml_yn;

    for (double i = 0.25; i < domain.geo.n_pml_yp - 0.75; i++)
        rho_e_yp[static_cast<int>(i)] = i / domain.geo.n_pml_yp;

    for (double i = domain.geo.n_pml_yn - 0.25; i > 0.75; i--)
        rho_m_yn[static_cast<int>(i)] = i / domain.geo.n_pml_yn;

    for (double i = 0.75; i < domain.geo.n_pml_yp - 0.25; i++)
        rho_m_yp[static_cast<int>(i)] = i / domain.geo.n_pml_yp;

    // 设置sigma参数
    sigma_max_xn =
        sigma_factor_xn * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_xn * mu_r_pml_xn) * domain.geo.dh[0]);
    sigma_max_xp =
        sigma_factor_xp * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_xp * mu_r_pml_xp) * domain.geo.dh[0]);

    sigma_max_yn =
        sigma_factor_yn * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_yn * mu_r_pml_yn) * domain.geo.dh[1]);
    sigma_max_yp =
        sigma_factor_yp * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_yp * mu_r_pml_yp) * domain.geo.dh[1]);

    sigma_e_xn = sigma_max_xn * pow(rho_e_xn, m_pml);
    sigma_m_xn =
        Const::MU_0 / Const::EPS_0 * sigma_max_xn * pow(rho_m_xn, m_pml);
    //    for (int i = 0; i < domain.geo.n_pml_xn; i++) {
    //        sigma_e_xn[i] = sigma_max_xn * pow(rho_e_xn[i], m_pml);
    //        sigma_m_xn[i] =
    //            Const::MU_0 / Const::EPS_0 * sigma_max_xn * pow(rho_m_xn[i],
    //            m_pml);
    //    }

    sigma_e_xp = sigma_max_xp * pow(rho_e_xp, m_pml);
    sigma_m_xp =
        Const::MU_0 / Const::EPS_0 * sigma_max_xp * pow(rho_m_xp, m_pml);

    sigma_e_yn = sigma_max_yn * pow(rho_e_yn, m_pml);
    sigma_m_yn =
        Const::MU_0 / Const::EPS_0 * sigma_max_yn * pow(rho_m_yn, m_pml);

    sigma_e_yp = sigma_max_yp * pow(rho_e_yp, m_pml);
    sigma_m_yp =
        Const::MU_0 / Const::EPS_0 * sigma_max_yp * pow(rho_m_yp, m_pml);

    // 设置kappa参数
    kappa_e_xn = 1 + (kappa_max_xn - 1) * pow(rho_e_xn, m_pml);
    kappa_m_xn = 1 + (kappa_max_xn - 1) * pow(rho_m_xn, m_pml);

    kappa_e_xp = 1 + (kappa_max_xp - 1) * pow(rho_e_xp, m_pml);
    kappa_m_xp = 1 + (kappa_max_xp - 1) * pow(rho_m_xp, m_pml);

    kappa_e_yn = 1 + (kappa_max_yn - 1) * pow(rho_e_yn, m_pml);
    kappa_m_yn = 1 + (kappa_max_yn - 1) * pow(rho_m_yn, m_pml);

    kappa_e_yp = 1 + (kappa_max_yp - 1) * pow(rho_e_yp, m_pml);
    kappa_m_yp = 1 + (kappa_max_yp - 1) * pow(rho_m_yp, m_pml);

    // 设置alpha参数
    alpha_e_xn = alpha_max * (1 - rho_e_xn);
    alpha_m_xn = Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_xn);

    alpha_e_xp = alpha_max * (1 - rho_e_xp);
    alpha_m_xp = Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_xp);

    alpha_e_yn = alpha_max * (1 - rho_e_yn);
    alpha_m_yn = Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_yn);

    alpha_e_yp = alpha_max * (1 - rho_e_yp);
    alpha_m_yp = Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_yp);

    // b
    b_e_xn_v = exp(-(sigma_e_xn / kappa_e_xn + alpha_e_xn) * domain.getDt() /
                   Const::EPS_0);
    b_m_xn_v = exp(-(sigma_m_xn / kappa_m_xn + alpha_m_xn) * domain.getDt() /
                   Const::MU_0);

    b_e_xp_v = exp(-(sigma_e_xp / kappa_e_xp + alpha_e_xp) * domain.getDt() /
                   Const::EPS_0);
    b_m_xp_v = exp(-(sigma_m_xp / kappa_m_xp + alpha_m_xp) * domain.getDt() /
                   Const::MU_0);

    b_e_yn_v = exp(-(sigma_e_yn / kappa_e_yn + alpha_e_yn) * domain.getDt() /
                   Const::EPS_0);
    b_m_yn_v = exp(-(sigma_m_yn / kappa_m_yn + alpha_m_yn) * domain.getDt() /
                   Const::MU_0);

    b_e_yp_v = exp(-(sigma_e_yp / kappa_e_yp + alpha_e_yp) * domain.getDt() /
                   Const::EPS_0);
    b_m_yp_v = exp(-(sigma_m_yp / kappa_m_yp + alpha_m_yp) * domain.getDt() /
                   Const::MU_0);

    // a
    a_e_xn_v = sigma_e_xn * (b_e_xn_v - 1) /
               (domain.geo.dh[0] * kappa_e_xn *
                (sigma_e_xn + alpha_e_xn * kappa_e_xn));
    a_m_xn_v = sigma_m_xn * (b_m_xn_v - 1) /
               (domain.geo.dh[0] * kappa_m_xn *
                (sigma_m_xn + alpha_m_xn * kappa_m_xn));

    a_e_xp_v = sigma_e_xp * (b_e_xp_v - 1) /
               (domain.geo.dh[0] * kappa_e_xp *
                (sigma_e_xp + alpha_e_xp * kappa_e_xp));
    a_m_xp_v = sigma_m_xp * (b_m_xp_v - 1) /
               (domain.geo.dh[0] * kappa_m_xp *
                (sigma_m_xp + alpha_m_xp * kappa_m_xp));

    a_e_yn_v = sigma_e_yn * (b_e_yn_v - 1) /
               (domain.geo.dh[0] * kappa_e_yn *
                (sigma_e_yn + alpha_e_yn * kappa_e_yn));
    a_m_yn_v = sigma_m_yn * (b_m_yn_v - 1) /
               (domain.geo.dh[0] * kappa_m_yn *
                (sigma_m_yn + alpha_m_yn * kappa_m_yn));

    a_e_yp_v = sigma_e_yp * (b_e_yp_v - 1) /
               (domain.geo.dh[0] * kappa_e_yp *
                (sigma_e_yp + alpha_e_yp * kappa_e_yp));
    a_m_yp_v = sigma_m_yp * (b_m_yp_v - 1) /
               (domain.geo.dh[0] * kappa_m_yp *
                (sigma_m_yp + alpha_m_yp * kappa_m_yp));
}

void Solver::Allocate(Domain &domain) {
    rho_e_xn.resize(domain.geo.n_pml_xn);
    rho_e_xp.resize(domain.geo.n_pml_xp);
    rho_m_xn.resize(domain.geo.n_pml_xn);
    rho_m_xp.resize(domain.geo.n_pml_xp);

    rho_e_yn.resize(domain.geo.n_pml_yn);
    rho_e_yp.resize(domain.geo.n_pml_yp);
    rho_m_yn.resize(domain.geo.n_pml_yn);
    rho_m_yp.resize(domain.geo.n_pml_yp);

    // sigma
    sigma_e_xn.resize(domain.geo.n_pml_xn);
    sigma_m_xn.resize(domain.geo.n_pml_xn);
    sigma_e_xp.resize(domain.geo.n_pml_xp);
    sigma_m_xp.resize(domain.geo.n_pml_xp);

    sigma_e_yn.resize(domain.geo.n_pml_yn);
    sigma_m_yn.resize(domain.geo.n_pml_yn);
    sigma_e_yp.resize(domain.geo.n_pml_yp);
    sigma_m_yp.resize(domain.geo.n_pml_yp);

    // kappa
    kappa_e_xn.resize(domain.geo.n_pml_xn);
    kappa_e_xp.resize(domain.geo.n_pml_xp);
    kappa_m_xn.resize(domain.geo.n_pml_xn);
    kappa_m_xp.resize(domain.geo.n_pml_xp);

    kappa_e_yn.resize(domain.geo.n_pml_yn);
    kappa_e_yp.resize(domain.geo.n_pml_yp);
    kappa_m_yn.resize(domain.geo.n_pml_yn);
    kappa_m_yp.resize(domain.geo.n_pml_yp);

    // alpha
    alpha_e_xn.resize(domain.geo.n_pml_xn);
    alpha_e_xp.resize(domain.geo.n_pml_xp);
    alpha_m_xn.resize(domain.geo.n_pml_xn);
    alpha_m_xp.resize(domain.geo.n_pml_xp);

    alpha_e_yn.resize(domain.geo.n_pml_yn);
    alpha_e_yp.resize(domain.geo.n_pml_yp);
    alpha_m_yn.resize(domain.geo.n_pml_yn);
    alpha_m_yp.resize(domain.geo.n_pml_yp);

    // b
    b_e_xn_v.resize(domain.geo.n_pml_xn);
    b_e_xp_v.resize(domain.geo.n_pml_xp);
    b_m_xn_v.resize(domain.geo.n_pml_xn);
    b_m_xp_v.resize(domain.geo.n_pml_xp);

    b_e_yn_v.resize(domain.geo.n_pml_yn);
    b_e_yp_v.resize(domain.geo.n_pml_yp);
    b_m_yn_v.resize(domain.geo.n_pml_yn);
    b_m_yp_v.resize(domain.geo.n_pml_yp);

    // a
    a_e_xn_v.resize(domain.geo.n_pml_xn);
    a_e_xp_v.resize(domain.geo.n_pml_xp);
    a_m_xn_v.resize(domain.geo.n_pml_xn);
    a_m_xp_v.resize(domain.geo.n_pml_xp);

    a_e_yn_v.resize(domain.geo.n_pml_yn);
    a_e_yp_v.resize(domain.geo.n_pml_yp);
    a_m_yn_v.resize(domain.geo.n_pml_yn);
    a_m_yp_v.resize(domain.geo.n_pml_yp);

    // Coeff
    C_Dx_dyn_v.resize(domain.geo.n_pml_yn);
    C_Dx_dyp_v.resize(domain.geo.n_pml_yp);

    Phi_ex_yn.resize(domain.geo.ni, domain.geo.n_pml_yn);
    b_ex_yn.resize(domain.geo.ni, domain.geo.n_pml_yn);
    a_ex_yn.resize(domain.geo.ni, domain.geo.n_pml_yn);
    C_Dx_dyn.resize(domain.geo.ni, domain.geo.n_pml_yn);

    Phi_ex_yp.resize(domain.geo.ni, domain.geo.n_pml_yp);
    b_ex_yp.resize(domain.geo.ni, domain.geo.n_pml_yp);
    a_ex_yp.resize(domain.geo.ni, domain.geo.n_pml_yp);
    C_Dx_dyp.resize(domain.geo.ni, domain.geo.n_pml_yp);

    temp_dx.resize(domain.geo.ni);

    C_Dy_dxn_v.resize(domain.geo.n_pml_xn);
    C_Dy_dxp_v.resize(domain.geo.n_pml_xp);

    Phi_ey_xn.resize(domain.geo.n_pml_xn, domain.geo.nj);
    b_ey_xn.resize(domain.geo.n_pml_xn, domain.geo.nj);
    a_ey_xn.resize(domain.geo.n_pml_xn, domain.geo.nj);
    C_Dy_dxn.resize(domain.geo.n_pml_xn, domain.geo.nj);

    Phi_ey_xp.resize(domain.geo.n_pml_xp, domain.geo.nj);
    b_ey_xp.resize(domain.geo.n_pml_xp, domain.geo.nj);
    a_ey_xp.resize(domain.geo.n_pml_xp, domain.geo.nj);
    C_Dy_dxp.resize(domain.geo.n_pml_xp, domain.geo.nj);

    temp_dy.resize(domain.geo.nj);

    C_Bz_dxn_v.resize(domain.geo.n_pml_xn);
    C_Bz_dxp_v.resize(domain.geo.n_pml_xp);

    C_Bz_dyn_v.resize(domain.geo.n_pml_yn);
    C_Bz_dyp_v.resize(domain.geo.n_pml_yp);

    Phi_mz_xn.resize(domain.geo.n_pml_xn, domain.geo.nj);
    b_mz_xn.resize(domain.geo.n_pml_xn, domain.geo.nj);
    b_mz_xn.resize(domain.geo.n_pml_xn, domain.geo.nj);
    C_Bz_dxn.resize(domain.geo.n_pml_xn, domain.geo.nj);

    Phi_mz_xp.resize(domain.geo.n_pml_xp, domain.geo.nj);
    b_mz_xp.resize(domain.geo.n_pml_xp, domain.geo.nj);
    b_mz_xp.resize(domain.geo.n_pml_xp, domain.geo.nj);
    C_Bz_dxp.resize(domain.geo.n_pml_xp, domain.geo.nj);

    temp_bz.resize(domain.geo.nj);

    Phi_mz_yn.resize(domain.geo.ni, domain.geo.n_pml_yn);
    b_mz_yn.resize(domain.geo.ni, domain.geo.n_pml_yn);
    b_mz_yn.resize(domain.geo.ni, domain.geo.n_pml_yn);
    C_Bz_dyn.resize(domain.geo.ni, domain.geo.n_pml_yn);

    Phi_mz_xp.resize(domain.geo.ni, domain.geo.n_pml_yp);
    b_mz_xp.resize(domain.geo.ni, domain.geo.n_pml_yp);
    b_mz_xp.resize(domain.geo.ni, domain.geo.n_pml_yp);
    C_Bz_dxp.resize(domain.geo.ni, domain.geo.n_pml_yp);
}

/*solves Poisson equation using Gauss-Seidel*/
bool Solver::solve() {
    // references to avoid having to write world.phi
    Field<double> &phi = domain.phi;
    Field<double> &rho = domain.rho;

    // precompute 1/(dx^2)
    Vec2d dh = domain.geo.dh;
    double idx2 = 1.0 / (dh[0] * dh[0]);
    double idy2 = 1.0 / (dh[1] * dh[1]);

    double L2 = 0;  // norm
    bool converged = false;

    /*solve potential*/
    for (unsigned it = 0; it < max_solver_it; it++) {
        for (int i = 1; i < domain.geo.ni - 1; i++)
            for (int j = 1; j < domain.geo.nj - 1; j++) {
                // standard internal open node
                double phi_new = (rho[i][j] / Const::EPS_0 +
                                  idx2 * (phi[i - 1][j] + phi[i + 1][j]) +
                                  idy2 * (phi[i][j - 1] + phi[i][j + 1])) /
                                 (2 * idx2 + 2 * idy2);

                /*SOR*/
                phi[i][j] = phi[i][j] + 1.4 * (phi_new - phi[i][j]);
            }

        /*check for convergence*/
        if (it % 25 == 0) {
            double sum = 0;
            for (int i = 1; i < domain.geo.ni - 1; i++)
                for (int j = 1; j < domain.geo.nj - 1; j++) {
                    double R = -phi[i][j] * (2 * idx2 + 2 * idy2) +
                               rho[i][j] / Const::EPS_0 +
                               idx2 * (phi[i - 1][j] + phi[i + 1][j]) +
                               idy2 * (phi[i][j - 1] + phi[i][j + 1]);

                    sum += R * R;
                }

            L2 = sqrt(sum / (domain.geo.ni * domain.geo.nj));
            if (L2 < tolerance) {
                converged = true;
                break;
            }
        }
    }

    if (!converged) {
        std::cerr << "GS failed to converge, L2=" << L2 << std::endl;
    }
    return converged;
}

/*computes electric field = -gradient(phi) using 2nd order differencing*/
void Solver::computeEF() {
    // reference to phi to avoid needing to write world.phi
    Field<double> &phi = domain.phi;

    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];

    for (int i = 0; i < domain.geo.ni; i++)
        for (int j = 0; j < domain.geo.nj; j++) {
            Vec3d &E = domain.E[i][j];  // reference to (i,j,k) E
                                        // vec3

            /*x component*/
            if (i == 0) /*forward*/
                E[0] = -(-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) /
                       (2 * dx);
            else if (i == domain.geo.ni - 1) /*backward*/
                E[0] = -(phi[i - 2][j] - 4 * phi[i - 1][j] + 3 * phi[i][j]) /
                       (2 * dx);
            else /*central*/
                E[0] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dx);

            /*y component*/
            if (j == 0)
                E[1] = -(-3 * phi[i][j] + 4 * phi[i][j + 1] - phi[i][j + 2]) /
                       (2 * dy);
            else if (j == domain.geo.nj - 1)
                E[1] = -(phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) /
                       (2 * dy);
            else
                E[1] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dy);
        }
}

void Solver::UpdateHx2d() {
    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = domain.getDt();

#pragma omp parallel for
    for (int i = 0; i < domain.geo.ni; i++) {
#pragma omp parallel for
        for (int j = 0; j < domain.geo.nj - 1; j++) {
            Field<Vec3d> &E = domain.E;
            Field<Vec3d> &H = domain.H;
            Field<Vec3d> &J = domain.J;

            H[i][j][0] = H[i][j][0] - dt * (E[i][j + 1][2] - E[i][j][2]) /
                                          (Const::MU_0 * dy);
        }
    }
}

void Solver::UpdateHy2d() {
    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = domain.getDt();

#pragma omp parallel for
    for (int i = 0; i < domain.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 0; j < domain.geo.nj; j++) {
            Field<Vec3d> &E = domain.E;
            Field<Vec3d> &H = domain.H;
            Field<Vec3d> &J = domain.J;

            H[i][j][1] = H[i][j][1] + dt * (E[i + 1][j][2] - E[i][j][2]) /
                                          (Const::MU_0 * dx);
        }
    }
}

void Solver::UpdateHz2d() {
    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = domain.getDt();

#pragma omp parallel for
    for (int i = 0; i < domain.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 0; j < domain.geo.nj - 1; j++) {
            Field<Vec3d> &E = domain.E;
            Field<Vec3d> &H = domain.H;
            Field<Vec3d> &J = domain.J;

            H[i][j][2] =
                H[i][j][2] -
                dt * (E[i + 1][j][1] - E[i][j][1]) / (Const::MU_0 * dx) +
                dt * (E[i][j + 1][0] - E[i][j][0]) / (Const::MU_0 * dy);
        }
    }
}

void Solver::UpdateEx2d() {
    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = domain.getDt();

#pragma omp parallel for
    for (int i = 1; i < domain.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 1; j < domain.geo.nj - 1; j++) {
            Field<Vec3d> &E = domain.E;
            Field<Vec3d> &H = domain.H;
            Field<Vec3d> &J = domain.J;

            E[i][j][0] =
                E[i][j][0] +
                dt * (H[i][j][2] - H[i][j - 1][2]) / (Const::EPS_0 * dy) -
                J[i][j][0] / Const::EPS_0;
        }
    }
}

void Solver::UpdateEy2d() {
    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = domain.getDt();

#pragma omp parallel for
    for (int i = 1; i < domain.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 1; j < domain.geo.nj - 1; j++) {
            Field<Vec3d> &E = domain.E;
            Field<Vec3d> &H = domain.H;
            Field<Vec3d> &J = domain.J;

            E[i][j][1] =
                E[i][j][1] -
                dt * (H[i][j][2] - H[i - 1][j][2]) / (Const::EPS_0 * dx) -
                J[i][j][1] / Const::EPS_0;
        }
    }
}

void Solver::UpdateEz2d() {
    Vec2d dh = domain.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = domain.getDt();

#pragma omp parallel for
    for (int i = 1; i < domain.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 1; j < domain.geo.nj - 1; j++) {
            Field<Vec3d> &E = domain.E;
            Field<Vec3d> &H = domain.H;
            Field<Vec3d> &J = domain.J;

            E[i][j][2] =
                E[i][j][2] +
                dt * (H[i][j][1] - H[i - 1][j][1]) / (Const::EPS_0 * dx) -
                dt * (H[i][j][0] - H[i][j - 1][0]) / (Const::EPS_0 * dy) -
                dt * J[i][j][2] / Const::EPS_0;
        }
    }
}

void Solver::UpdateElectromagnetic() {
    UpdateEz2d();
    // UpdateEx2d();
    // UpdateEy2d();
    // UpdateHz2d();
    UpdateHx2d();
    UpdateHy2d();
}
