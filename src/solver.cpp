#include "solver.h"

#include <Eigen/src/Core/util/IndexedViewHelper.h>
#include <math.h>

#include <iostream>

#include "const.h"
#include "field.h"

Solver::Solver(Domain &dm, int max_it, double tol)
    : dm(dm), max_solver_it(max_it), tolerance(tol) {
    Allocate(dm);
    CalculFdtdPara(dm);
    CalculFdtdCoeff(dm);
}
void Solver::CalculFdtdCoeff(Domain &dm) {
    CD_dt = dm.getDt();
    CD_dt_dx = CD_dt / dm.geo.dh[0];
    CD_dt_dy = CD_dt / dm.geo.dh[1];

    CB_dt = dm.getDt();
    CB_dt_dx = CB_dt / dm.geo.dh[0];
    CB_dt_dy = CB_dt / dm.geo.dh[1];

    C_Dx_dyn_v = CD_dt_dy / kappa_e_yn;
    C_Dx_dyp_v = CD_dt_dy / kappa_e_yp;

    // PML区域内Dx的更新公式系数
    for (int i = 0; i < dm.geo.n_pml_yn; ++i) {
        b_ex_yn.col(i).array() = b_e_yn_v(i);
        b_ex_yn.matrix();
        a_ex_yn.col(i) = a_e_yn_v(i);
        C_Dx_dyn.col(i) = C_Dx_dyn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_yp; ++i) {
        b_ex_yp.col(i) = b_e_yp_v(i);
        a_ex_yp.col(i) = a_e_yp_v(i);
        C_Dx_dyp.col(i) = C_Dx_dyp_v(i);
    }
    spdlog::debug("The coefficient of Dx has been calculated");

    // PML区域内Dy的公式更新系数
    C_Dy_dxn_v = CD_dt_dx / kappa_e_xn;
    C_Dy_dxp_v = CD_dt_dx / kappa_e_xp;

    for (int i = 0; i < dm.geo.n_pml_xn; ++i) {
        b_ey_xn.row(i) = b_e_xn_v(i);
        a_ey_xn.row(i) = a_e_xn_v(i);
        C_Dy_dxn.row(i) = C_Dy_dxn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_xp; ++i) {
        b_ey_xp.row(i) = b_e_xp_v(i);
        a_ey_xp.row(i) = a_e_xp_v(i);
        C_Dy_dxp.row(i) = C_Dy_dxp_v(i);
    }
    spdlog::debug("The coefficient of Dy has been calculated");

    // PML区域内Bz的更新公式系数
    C_Bz_dxn_v = CB_dt_dx / kappa_m_xn;
    C_Bz_dxp_v = CB_dt_dx / kappa_m_xp;

    C_Bz_dyn_v = CB_dt_dy / kappa_m_yn;
    C_Bz_dyp_v = CB_dt_dy / kappa_m_yp;

    for (int i = 0; i < dm.geo.n_pml_xn; i++) {
        b_mz_xn.row(i) = b_m_xn_v(i);
        a_mz_xn.row(i) = a_m_xn_v(i);
        C_Bz_dxn.row(i) = C_Bz_dxn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_xp; ++i) {
        b_mz_xp.row(i) = b_m_xp_v(i);
        a_mz_xp.row(i) = a_m_xp_v(i);
        C_Bz_dxp.row(i) = C_Bz_dxp_v(i);
    }

    for (int i = 0; i < dm.geo.n_pml_yn; ++i) {
        b_mz_yn.col(i) = b_m_yn_v(i);
        a_mz_yn.col(i) = a_m_yn_v(i);
        C_Bz_dyn.col(i) = C_Bz_dyn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_yp; ++i) {
        b_mz_xp.col(i) = b_m_yp_v(i);
        a_mz_xp.col(i) = a_m_yp_v(i);
        C_Bz_dxp.col(i) = C_Bz_dyp_v(i);
    }
    spdlog::debug("The coefficient of Bz has been calculated");

    spdlog::debug("All of the coefficients have been calculated");
}

void Solver::CalculFdtdPara(Domain &dm) {
    // 计算电磁场量节点离开PML内边界的距离
    // x轴方向
    for (double i = dm.geo.n_pml_xn - 0.75; i > 0.25; i--)
        rho_e_xn[static_cast<int>(i)] = i / dm.geo.n_pml_xn;

    for (double i = 0.25; i < dm.geo.n_pml_xp - 0.75; i++)
        rho_e_xp[static_cast<int>(i)] = i / dm.geo.n_pml_xp;

    for (double i = dm.geo.n_pml_xn - 0.25; i > 0.75; i--)
        rho_m_xn[static_cast<int>(i)] = i / dm.geo.n_pml_xn;

    for (double i = 0.75; i < dm.geo.n_pml_xp - 0.25; i++)
        rho_m_xp[static_cast<int>(i)] = i / dm.geo.n_pml_xp;

    // y轴方向
    for (double i = dm.geo.n_pml_yn - 0.75; i > 0.25; i--)
        rho_e_yn[static_cast<int>(i)] = i / dm.geo.n_pml_yn;

    for (double i = 0.25; i < dm.geo.n_pml_yp - 0.75; i++)
        rho_e_yp[static_cast<int>(i)] = i / dm.geo.n_pml_yp;

    for (double i = dm.geo.n_pml_yn - 0.25; i > 0.75; i--)
        rho_m_yn[static_cast<int>(i)] = i / dm.geo.n_pml_yn;

    for (double i = 0.75; i < dm.geo.n_pml_yp - 0.25; i++)
        rho_m_yp[static_cast<int>(i)] = i / dm.geo.n_pml_yp;

    // 设置sigma参数
    sigma_max_xn =
        sigma_factor_xn * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_xn * mu_r_pml_xn) * dm.geo.dh[0]);
    sigma_max_xp =
        sigma_factor_xp * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_xp * mu_r_pml_xp) * dm.geo.dh[0]);

    sigma_max_yn =
        sigma_factor_yn * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_yn * mu_r_pml_yn) * dm.geo.dh[1]);
    sigma_max_yp =
        sigma_factor_yp * (m_pml + 1) /
        (150 * Const::PI * sqrt(eps_r_pml_yp * mu_r_pml_yp) * dm.geo.dh[1]);

    sigma_e_xn = sigma_max_xn * pow(rho_e_xn, m_pml);
    sigma_m_xn =
        Const::MU_0 / Const::EPS_0 * sigma_max_xn * pow(rho_m_xn, m_pml);

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
    b_e_xn_v = exp(-(sigma_e_xn / kappa_e_xn + alpha_e_xn) * dm.getDt() /
                   Const::EPS_0);
    b_m_xn_v =
        exp(-(sigma_m_xn / kappa_m_xn + alpha_m_xn) * dm.getDt() / Const::MU_0);

    b_e_xp_v = exp(-(sigma_e_xp / kappa_e_xp + alpha_e_xp) * dm.getDt() /
                   Const::EPS_0);
    b_m_xp_v =
        exp(-(sigma_m_xp / kappa_m_xp + alpha_m_xp) * dm.getDt() / Const::MU_0);

    b_e_yn_v = exp(-(sigma_e_yn / kappa_e_yn + alpha_e_yn) * dm.getDt() /
                   Const::EPS_0);
    b_m_yn_v =
        exp(-(sigma_m_yn / kappa_m_yn + alpha_m_yn) * dm.getDt() / Const::MU_0);

    b_e_yp_v = exp(-(sigma_e_yp / kappa_e_yp + alpha_e_yp) * dm.getDt() /
                   Const::EPS_0);
    b_m_yp_v =
        exp(-(sigma_m_yp / kappa_m_yp + alpha_m_yp) * dm.getDt() / Const::MU_0);

    // a
    a_e_xn_v =
        sigma_e_xn * (b_e_xn_v - 1) /
        (dm.geo.dh[0] * kappa_e_xn * (sigma_e_xn + alpha_e_xn * kappa_e_xn));
    a_m_xn_v =
        sigma_m_xn * (b_m_xn_v - 1) /
        (dm.geo.dh[0] * kappa_m_xn * (sigma_m_xn + alpha_m_xn * kappa_m_xn));

    a_e_xp_v =
        sigma_e_xp * (b_e_xp_v - 1) /
        (dm.geo.dh[0] * kappa_e_xp * (sigma_e_xp + alpha_e_xp * kappa_e_xp));
    a_m_xp_v =
        sigma_m_xp * (b_m_xp_v - 1) /
        (dm.geo.dh[0] * kappa_m_xp * (sigma_m_xp + alpha_m_xp * kappa_m_xp));

    a_e_yn_v =
        sigma_e_yn * (b_e_yn_v - 1) /
        (dm.geo.dh[0] * kappa_e_yn * (sigma_e_yn + alpha_e_yn * kappa_e_yn));
    a_m_yn_v =
        sigma_m_yn * (b_m_yn_v - 1) /
        (dm.geo.dh[0] * kappa_m_yn * (sigma_m_yn + alpha_m_yn * kappa_m_yn));

    a_e_yp_v =
        sigma_e_yp * (b_e_yp_v - 1) /
        (dm.geo.dh[0] * kappa_e_yp * (sigma_e_yp + alpha_e_yp * kappa_e_yp));
    a_m_yp_v =
        sigma_m_yp * (b_m_yp_v - 1) /
        (dm.geo.dh[0] * kappa_m_yp * (sigma_m_yp + alpha_m_yp * kappa_m_yp));

    spdlog::debug("All of the parameters have been calculated");
}

void Solver::Allocate(Domain &dm) {
    rho_e_xn.resize(dm.geo.n_pml_xn);
    rho_e_xp.resize(dm.geo.n_pml_xp);
    rho_m_xn.resize(dm.geo.n_pml_xn);
    rho_m_xp.resize(dm.geo.n_pml_xp);

    rho_e_yn.resize(dm.geo.n_pml_yn);
    rho_e_yp.resize(dm.geo.n_pml_yp);
    rho_m_yn.resize(dm.geo.n_pml_yn);
    rho_m_yp.resize(dm.geo.n_pml_yp);

    // sigma
    sigma_e_xn.resize(dm.geo.n_pml_xn);
    sigma_m_xn.resize(dm.geo.n_pml_xn);
    sigma_e_xp.resize(dm.geo.n_pml_xp);
    sigma_m_xp.resize(dm.geo.n_pml_xp);

    sigma_e_yn.resize(dm.geo.n_pml_yn);
    sigma_m_yn.resize(dm.geo.n_pml_yn);
    sigma_e_yp.resize(dm.geo.n_pml_yp);
    sigma_m_yp.resize(dm.geo.n_pml_yp);

    // kappa
    kappa_e_xn.resize(dm.geo.n_pml_xn);
    kappa_e_xp.resize(dm.geo.n_pml_xp);
    kappa_m_xn.resize(dm.geo.n_pml_xn);
    kappa_m_xp.resize(dm.geo.n_pml_xp);

    kappa_e_yn.resize(dm.geo.n_pml_yn);
    kappa_e_yp.resize(dm.geo.n_pml_yp);
    kappa_m_yn.resize(dm.geo.n_pml_yn);
    kappa_m_yp.resize(dm.geo.n_pml_yp);

    // alpha
    alpha_e_xn.resize(dm.geo.n_pml_xn);
    alpha_e_xp.resize(dm.geo.n_pml_xp);
    alpha_m_xn.resize(dm.geo.n_pml_xn);
    alpha_m_xp.resize(dm.geo.n_pml_xp);

    alpha_e_yn.resize(dm.geo.n_pml_yn);
    alpha_e_yp.resize(dm.geo.n_pml_yp);
    alpha_m_yn.resize(dm.geo.n_pml_yn);
    alpha_m_yp.resize(dm.geo.n_pml_yp);

    // b
    b_e_xn_v.resize(dm.geo.n_pml_xn);
    b_e_xp_v.resize(dm.geo.n_pml_xp);
    b_m_xn_v.resize(dm.geo.n_pml_xn);
    b_m_xp_v.resize(dm.geo.n_pml_xp);

    b_e_yn_v.resize(dm.geo.n_pml_yn);
    b_e_yp_v.resize(dm.geo.n_pml_yp);
    b_m_yn_v.resize(dm.geo.n_pml_yn);
    b_m_yp_v.resize(dm.geo.n_pml_yp);

    // a
    a_e_xn_v.resize(dm.geo.n_pml_xn);
    a_e_xp_v.resize(dm.geo.n_pml_xp);
    a_m_xn_v.resize(dm.geo.n_pml_xn);
    a_m_xp_v.resize(dm.geo.n_pml_xp);

    a_e_yn_v.resize(dm.geo.n_pml_yn);
    a_e_yp_v.resize(dm.geo.n_pml_yp);
    a_m_yn_v.resize(dm.geo.n_pml_yn);
    a_m_yp_v.resize(dm.geo.n_pml_yp);

    // Coeff
    C_Dx_dyn_v.resize(dm.geo.n_pml_yn);
    C_Dx_dyp_v.resize(dm.geo.n_pml_yp);

    Phi_ex_yn.resize(dm.geo.ni, dm.geo.n_pml_yn);
    b_ex_yn.resize(dm.geo.ni, dm.geo.n_pml_yn);
    a_ex_yn.resize(dm.geo.ni, dm.geo.n_pml_yn);
    C_Dx_dyn.resize(dm.geo.ni, dm.geo.n_pml_yn);

    Phi_ex_yp.resize(dm.geo.ni, dm.geo.n_pml_yp);
    b_ex_yp.resize(dm.geo.ni, dm.geo.n_pml_yp);
    a_ex_yp.resize(dm.geo.ni, dm.geo.n_pml_yp);
    C_Dx_dyp.resize(dm.geo.ni, dm.geo.n_pml_yp);

    C_Dy_dxn_v.resize(dm.geo.n_pml_xn);
    C_Dy_dxp_v.resize(dm.geo.n_pml_xp);

    Phi_ey_xn.resize(dm.geo.n_pml_xn, dm.geo.nj);
    b_ey_xn.resize(dm.geo.n_pml_xn, dm.geo.nj);
    a_ey_xn.resize(dm.geo.n_pml_xn, dm.geo.nj);
    C_Dy_dxn.resize(dm.geo.n_pml_xn, dm.geo.nj);

    Phi_ey_xp.resize(dm.geo.n_pml_xp, dm.geo.nj);
    b_ey_xp.resize(dm.geo.n_pml_xp, dm.geo.nj);
    a_ey_xp.resize(dm.geo.n_pml_xp, dm.geo.nj);
    C_Dy_dxp.resize(dm.geo.n_pml_xp, dm.geo.nj);

    C_Bz_dxn_v.resize(dm.geo.n_pml_xn);
    C_Bz_dxp_v.resize(dm.geo.n_pml_xp);

    C_Bz_dyn_v.resize(dm.geo.n_pml_yn);
    C_Bz_dyp_v.resize(dm.geo.n_pml_yp);

    Phi_mz_xn.resize(dm.geo.n_pml_xn, dm.geo.nj);
    b_mz_xn.resize(dm.geo.n_pml_xn, dm.geo.nj);
    a_mz_xn.resize(dm.geo.n_pml_xn, dm.geo.nj);
    C_Bz_dxn.resize(dm.geo.n_pml_xn, dm.geo.nj);

    Phi_mz_xp.resize(dm.geo.n_pml_xp, dm.geo.nj);
    b_mz_xp.resize(dm.geo.n_pml_xp, dm.geo.nj);
    a_mz_xp.resize(dm.geo.n_pml_xp, dm.geo.nj);
    C_Bz_dxp.resize(dm.geo.n_pml_xp, dm.geo.nj);

    Phi_mz_yn.resize(dm.geo.ni, dm.geo.n_pml_yn);
    b_mz_yn.resize(dm.geo.ni, dm.geo.n_pml_yn);
    a_mz_yn.resize(dm.geo.ni, dm.geo.n_pml_yn);
    C_Bz_dyn.resize(dm.geo.ni, dm.geo.n_pml_yn);

    Phi_mz_xp.resize(dm.geo.ni, dm.geo.n_pml_yp);
    b_mz_xp.resize(dm.geo.ni, dm.geo.n_pml_yp);
    a_mz_xp.resize(dm.geo.ni, dm.geo.n_pml_yp);
    C_Bz_dxp.resize(dm.geo.ni, dm.geo.n_pml_yp);

    spdlog::debug("All of the parameters have been allocated");
}

/*solves Poisson equation using Gauss-Seidel*/
bool Solver::solve() {
    // references to avoid having to write world.phi
    Field<double> &phi = dm.phi;
    Field<double> &rho = dm.rho;

    // precompute 1/(dx^2)
    Vec2d dh = dm.geo.dh;
    double idx2 = 1.0 / (dh[0] * dh[0]);
    double idy2 = 1.0 / (dh[1] * dh[1]);

    double L2 = 0;  // norm
    bool converged = false;

    /*solve potential*/
    for (unsigned it = 0; it < max_solver_it; it++) {
        for (int i = 1; i < dm.geo.ni - 1; i++)
            for (int j = 1; j < dm.geo.nj - 1; j++) {
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
            for (int i = 1; i < dm.geo.ni - 1; i++)
                for (int j = 1; j < dm.geo.nj - 1; j++) {
                    double R = -phi[i][j] * (2 * idx2 + 2 * idy2) +
                               rho[i][j] / Const::EPS_0 +
                               idx2 * (phi[i - 1][j] + phi[i + 1][j]) +
                               idy2 * (phi[i][j - 1] + phi[i][j + 1]);

                    sum += R * R;
                }

            L2 = sqrt(sum / (dm.geo.ni * dm.geo.nj));
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
    Field<double> &phi = dm.phi;

    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];

    for (int i = 0; i < dm.geo.ni; i++)
        for (int j = 0; j < dm.geo.nj; j++) {
            Vec3d &E = dm.E[i][j];  // reference to (i,j,k) E
                                    // vec3

            /*x component*/
            if (i == 0) /*forward*/
                E[0] = -(-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) /
                       (2 * dx);
            else if (i == dm.geo.ni - 1) /*backward*/
                E[0] = -(phi[i - 2][j] - 4 * phi[i - 1][j] + 3 * phi[i][j]) /
                       (2 * dx);
            else /*central*/
                E[0] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dx);

            /*y component*/
            if (j == 0)
                E[1] = -(-3 * phi[i][j] + 4 * phi[i][j + 1] - phi[i][j + 2]) /
                       (2 * dy);
            else if (j == dm.geo.nj - 1)
                E[1] = -(phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) /
                       (2 * dy);
            else
                E[1] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dy);
        }
}

void Solver::UpdateHx2d() {
    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = dm.getDt();

#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni; i++) {
#pragma omp parallel for
        for (int j = 0; j < dm.geo.nj - 1; j++) {
            Field<Vec3d> &E = dm.E;
            Field<Vec3d> &H = dm.H;
            Field<Vec3d> &J = dm.J;

            H[i][j][0] = H[i][j][0] - dt * (E[i][j + 1][2] - E[i][j][2]) /
                                          (Const::MU_0 * dy);
        }
    }
}

void Solver::UpdateHy2d() {
    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = dm.getDt();

#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 0; j < dm.geo.nj; j++) {
            Field<Vec3d> &E = dm.E;
            Field<Vec3d> &H = dm.H;
            Field<Vec3d> &J = dm.J;

            H[i][j][1] = H[i][j][1] + dt * (E[i + 1][j][2] - E[i][j][2]) /
                                          (Const::MU_0 * dx);
        }
    }
}

void Solver::UpdateHz2d() {
    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = dm.getDt();

#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 0; j < dm.geo.nj - 1; j++) {
            Field<Vec3d> &E = dm.E;
            Field<Vec3d> &H = dm.H;
            Field<Vec3d> &J = dm.J;

            H[i][j][2] =
                H[i][j][2] -
                dt * (E[i + 1][j][1] - E[i][j][1]) / (Const::MU_0 * dx) +
                dt * (E[i][j + 1][0] - E[i][j][0]) / (Const::MU_0 * dy);
        }
    }
}

void Solver::UpdateEx2d() {
    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = dm.getDt();

#pragma omp parallel for
    for (int i = 1; i < dm.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 1; j < dm.geo.nj - 1; j++) {
            Field<Vec3d> &E = dm.E;
            Field<Vec3d> &H = dm.H;
            Field<Vec3d> &J = dm.J;

            E[i][j][0] =
                E[i][j][0] +
                dt * (H[i][j][2] - H[i][j - 1][2]) / (Const::EPS_0 * dy) -
                J[i][j][0] / Const::EPS_0;
        }
    }
}

void Solver::UpdateEy2d() {
    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = dm.getDt();

#pragma omp parallel for
    for (int i = 1; i < dm.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 1; j < dm.geo.nj - 1; j++) {
            Field<Vec3d> &E = dm.E;
            Field<Vec3d> &H = dm.H;
            Field<Vec3d> &J = dm.J;

            E[i][j][1] =
                E[i][j][1] -
                dt * (H[i][j][2] - H[i - 1][j][2]) / (Const::EPS_0 * dx) -
                J[i][j][1] / Const::EPS_0;
        }
    }
}

void Solver::UpdateEz2d() {
    Vec2d dh = dm.geo.dh;
    double dx = dh[0];
    double dy = dh[1];
    double dt = dm.getDt();

#pragma omp parallel for
    for (int i = 1; i < dm.geo.ni - 1; i++) {
#pragma omp parallel for
        for (int j = 1; j < dm.geo.nj - 1; j++) {
            Field<Vec3d> &E = dm.E;
            Field<Vec3d> &H = dm.H;
            Field<Vec3d> &J = dm.J;

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

void Solver::UpdateBoundary(Domain &dm, double I, double f) {
    double t = dm.getTime();

// 设置电流边界
#pragma omp parallel for
    for (int i = 6; i < 20; i++) {
        for (int j = 6; j < 20; j++) {
            //            dm.J[5][j][2] = -I * sin(2 * Const::PI * f * t);
            //            dm.J[20][j][2] = I * sin(2 * Const::PI * f * t);
            //            dm.J[i][5][2] = -I * cos(2 * Const::PI * f * t);
            //            dm.J[i][20][2] = I * cos(2 * Const::PI * f * t);
            dm.Jz(5, j) = -I * sin(2 * Const::PI * f * t);
            dm.Jz(20, j) = I * sin(2 * Const::PI * f * t);
            dm.Jz(i, 5) = -I * cos(2 * Const::PI * f * t);
            dm.Jz(i, 20) = I * cos(2 * Const::PI * f * t);
        }
    }
    spdlog::debug("Current has been updated");

    // Dx更新主循环程序
    // 负y轴方向CMPL区域中的Dx更新程序
    auto r = seq(0, last);
    auto knD = seq(2, dm.geo.n_pml_yn + 1);
    auto knD_pre = seq(1, dm.geo.n_pml_yn);

    Phi_ex_yn = b_ex_yn * Phi_ex_yn + a_ex_yn * (dm.Hy(seq(0, last), knD) -
                                                 dm.Hy(seq(0, last), knD_pre));
    std::cout << Phi_ex_yn << std::endl;
    dm.Dx(seq(0, last), knD) =
        dm.Dx(seq(0, last), knD) -
        C_Dx_dyn * (dm.Hy(seq(0, last), knD) - dm.Hy(seq(0, last), knD_pre)) -
        CD_dt * Phi_ex_yn;

    spdlog::debug("Dx PML y- update successful");

    // 正y轴方向的CPML区域中的Dx更新程序
    auto knP = seq(dm.geo.nj - dm.geo.n_pml_yp + 1, dm.geo.nj);
    auto knP_pre = seq(dm.geo.nj - dm.geo.n_pml_yp, dm.geo.nj - 1);
#pragma omp parallel for
    for (int i = dm.geo.nj - dm.geo.n_pml_yp + 1; i < dm.geo.nj; i++) {
        Phi_ex_yp =
            b_ex_yp * Phi_ex_yp + a_ex_yp * (dm.Hy.col(i) - dm.Hy.col(i - 1));
        dm.Dx.col(i) = dm.Dx.col(i) -
                       C_Dx_dyp * (dm.Hy.col(i) - dm.Hy.col(i - 1)) -
                       CD_dt * Phi_ex_yp;
    }
    spdlog::debug("Dx PML y+ update successful");

    // 非PML域更新Dx
    auto kD = seq(dm.geo.n_pml_yn + 2, dm.geo.nj - dm.geo.n_pml_yn);
    auto kD_pre = seq(dm.geo.n_pml_yn + 1, dm.geo.nj - dm.geo.n_pml_yn - 1);
#pragma omp parallel for
    for (int i = dm.geo.n_pml_yn + 2; i < dm.geo.nj - dm.geo.n_pml_yp; i++) {
        dm.Dx.col(i) =
            dm.Dx.col(i) - CD_dt_dy * (dm.Hy.col(i) - dm.Hy.col(i - 1));
    }
    spdlog::debug("Dx  update successful");

    // Dy更新主循环程序
    // 负x轴方向PML
#pragma omp parallel for
    for (int i = 2; i < dm.geo.n_pml_xn + 1; ++i) {
        Phi_ey_xn =
            b_ey_xn * Phi_ey_xn + a_ey_xn * (dm.Hy.row(i) - dm.Hy.row(i - 1));
        dm.Dy.row(i) = dm.Dy.row(i) +
                       C_Dy_dxn * (dm.Hy.row(i) - dm.Hy.row(i - 1)) +
                       CD_dt * Phi_ey_xn;
    }
    spdlog::debug("Dy PML x- update successful");

    // 正y轴方向PML
#pragma omp parallel for
    for (int i = dm.geo.ni - dm.geo.n_pml_xp + 1; i < dm.geo.ni; i++) {
        Phi_ey_xp =
            b_ey_xp * Phi_ey_xp + a_ey_xp * (dm.Hy.row(i) - dm.Hy.row(i - 1));
        dm.Dy.row(i) = dm.Dy.row(i) +
                       C_Dy_dxp * (dm.Hy.row(i) - dm.Hy.row(i - 1)) +
                       CD_dt * Phi_ey_xp;
    }
    spdlog::debug("Dy PML x+ update successful");

    // 非PML
#pragma omp parallel for
    for (int i = dm.geo.n_pml_xn + 2; i < dm.geo.ni - dm.geo.n_pml_xp; i++) {
        dm.Dy.row(i) = dm.Dy.row(i) + CD_dt_dx * (dm.Hy(i) - dm.Hy(i - 1));
    }
    spdlog::debug("Dx update successful");

    // Bz更新
    // x
    // 负PML
#pragma omp parallel for
    for (int i = 1; i < dm.geo.n_pml_xn; i++) {
        Phi_mz_xn =
            b_mz_xn * Phi_mz_xn + a_mz_xn * (dm.Dx.row(i + 1) - dm.Dx.row(i));
        dm.Hz.row(i) = dm.Hz.row(i) +
                       C_Bz_dxn * (dm.Dx.row(i + 1) - dm.Dx.row(i)) +
                       CB_dt * Phi_mz_xn;
    }
    spdlog::debug("Bz PML x- update successful");

    // 正PML
#pragma omp parallel for
    for (int i = dm.geo.ni - dm.geo.n_pml_xp + 1; i < dm.geo.ni; i++) {
        Phi_mz_xp =
            b_mz_xp * Phi_mz_xp + a_mz_xp * (dm.Dx.row(i + 1) - dm.Dx.row(i));
        dm.Hz.row(i) =
            dm.Hz.row(i) + CB_dt_dx * (dm.Dx.row(i + 1) - dm.Dx.row(i));
    }
    spdlog::debug("Bz PML x+ update successful");

    // 非PML
#pragma omp parallel for
    for (int i = dm.geo.n_pml_xn + 1; i < dm.geo.ni - dm.geo.n_pml_xp; ++i) {
        dm.Hz.row(i) =
            dm.Hz.row(i) + CB_dt_dx * (dm.Dx.row(i + 1) - dm.Dx.row(i));
    }
    spdlog::debug("Bz y update successful");

    // y
    // 负PML
#pragma omp parallel for
    for (int i = 1; i < dm.geo.n_pml_yn; i++) {
        Phi_mz_yn =
            b_mz_yn * Phi_mz_yn + a_mz_yn * (dm.Dy.col(i + 1) - dm.Dy.col(i));
        dm.Hz.col(i) = dm.Hz.col(i) -
                       C_Bz_dyn * (dm.Dy.col(i + 1) - dm.Dy.col(i)) -
                       CB_dt * Phi_mz_yn;
    }
    spdlog::debug("Bz PML y- update successful");

    // 正PML
#pragma omp parallel for
    for (int i = dm.geo.nj - dm.geo.n_pml_yp + 1; i < dm.geo.nj; i++) {
        Phi_mz_yp =
            b_mz_yp * Phi_mz_yp + a_mz_yp * (dm.Dy.col(i + 1) - dm.Dy.col(i));
        dm.Hz.col(i) = dm.Hz.col(i) -
                       C_Bz_dyp * (dm.Dy.col(i + 1) - dm.Dy.col(i)) -
                       CB_dt * Phi_mz_yp;
    }
    spdlog::debug("Bz PML y+ update successful");

    // 非PML
#pragma omp parallel for
    for (int i = dm.geo.n_pml_yn + 1; i < dm.geo.nj - dm.geo.n_pml_yp; ++i) {
        dm.Hx.col(i) =
            dm.Hx.col(i) - CB_dt_dy * (dm.Dy.col(i + 1) - dm.Dy.col(i));
    }
    spdlog::debug("Bz y update successful");
}
