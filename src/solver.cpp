#include "solver.h"

Solver::Solver(Domain &dm, int max_it, double tol)
    : dm(dm), max_solver_it(max_it), tolerance(tol) {
    spdlog::debug("All of the parameters have been allocated");
    InitFdtdPara();
    CalculFdtdCoeff();
}

void Solver::CalculFdtdCoeff() {
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
        b_ex_yn.col(i) = b_e_yn_v(i);
        a_ex_yn.col(i) = a_e_yn_v(i);
        C_Dx_dyn.col(i) = C_Dx_dyn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_yp; ++i) {
        b_ex_yp.col(i) = b_e_yp_v(i);
        a_ex_yp.col(i) = a_e_yp_v(i);
        C_Dx_dyp.col(i) = C_Dx_dyp_v(i);
    }

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

    // Dz
    C_Dz_dxn_v = CD_dt_dx / kappa_m_xn;
    C_Dz_dxp_v = CD_dt_dx / kappa_m_xp;

    C_Dz_dyn_v = CD_dt_dy / kappa_m_yn;
    C_Dz_dyp_v = CD_dt_dy / kappa_m_yp;

    for (int i = 0; i < dm.geo.n_pml_xn; i++) {
        b_ez_xn.row(i) = b_e_xn_v(i);
        a_ez_xn.row(i) = a_e_xn_v(i);
        C_Dz_dxn.row(i) = C_Dz_dxn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_xp; i++) {
        b_ez_xp.row(i) = b_e_xp_v(i);
        a_ez_xp.row(i) = a_e_xp_v(i);
        C_Dz_dxp.row(i) = C_Dz_dxp_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_yn; ++i) {
        b_ez_yn.col(i) = b_e_yn_v(i);
        a_ez_yn.col(i) = a_e_yn_v(i);
        C_Dz_dyn.col(i) = C_Dz_dyn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_yp; ++i) {
        b_ez_yp.col(i) = b_e_yp_v(i);
        a_ez_yp.col(i) = a_e_yp_v(i);
        C_Dz_dyp.col(i) = C_Dz_dyp_v(i);
    }

    C_Bx_dyn_v = CB_dt_dy / kappa_e_yn;
    C_Bx_dyp_v = CB_dt_dy / kappa_e_yp;

    // PML区域内Dx的更新公式系数
    for (int i = 0; i < dm.geo.n_pml_yn; ++i) {
        b_mx_yn.col(i) = b_m_yn_v(i);
        a_mx_yn.col(i) = a_m_yn_v(i);
        C_Bx_dyn.col(i) = C_Bx_dyn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_yp; ++i) {
        b_mx_yp.col(i) = b_m_yp_v(i);
        a_mx_yp.col(i) = a_m_yp_v(i);
        C_Bx_dyp.col(i) = C_Bx_dyp_v(i);
    }

    // PML区域内Dy的公式更新系数
    C_By_dxn_v = CB_dt_dx / kappa_e_xn;
    C_By_dxp_v = CB_dt_dx / kappa_e_xp;

    for (int i = 0; i < dm.geo.n_pml_xn; ++i) {
        b_my_xn.row(i) = b_m_xn_v(i);
        a_my_xn.row(i) = a_m_xn_v(i);
        C_By_dxn.row(i) = C_By_dxn_v(i);
    }
    for (int i = 0; i < dm.geo.n_pml_xp; ++i) {
        b_my_xp.row(i) = b_m_xp_v(i);
        a_my_xp.row(i) = a_m_xp_v(i);
        C_By_dxp.row(i) = C_By_dxp_v(i);
    }

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
    for (int i = 0; i < dm.geo.n_pml_xp; i++) {
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
        b_mz_yp.col(i) = b_m_yp_v(i);
        a_mz_yp.col(i) = a_m_yp_v(i);
        C_Bz_dyp.col(i) = C_Bz_dyp_v(i);
    }

    spdlog::debug("All of the coefficients have been calculated");
}

void Solver::InitFdtdPara() {
    // 计算电磁场量节点离开PML内边界的距离
    // x轴方向
    for (double i = dm.geo.n_pml_xn - 0.75; i >= 0.25; i--) {
        rho_e_xn(dm.geo.n_pml_xn - static_cast<int>(i) - 1) =
            i / dm.geo.n_pml_xn;
    }

    for (double i = 0.25; i <= dm.geo.n_pml_xp - 0.75; i++) {
        rho_e_xp(static_cast<int>(i)) = i / dm.geo.n_pml_xp;
    }

    for (double i = dm.geo.n_pml_xn - 0.25; i >= 0.75; i--) {
        rho_m_xn(dm.geo.n_pml_xn - static_cast<int>(i) - 1) =
            i / dm.geo.n_pml_xn;
    }

    for (double i = 0.75; i <= dm.geo.n_pml_xp - 0.25; i++) {
        rho_m_xp(static_cast<int>(i)) = i / dm.geo.n_pml_xp;
    }

    // y轴方向
    for (double i = dm.geo.n_pml_yn - 0.75; i >= 0.25; i--) {
        rho_e_yn(dm.geo.n_pml_yn - static_cast<int>(i) - 1) =
            i / dm.geo.n_pml_yn;
    }

    for (double i = 0.25; i <= dm.geo.n_pml_yp - 0.75; i++) {
        rho_e_yp(static_cast<int>(i)) = i / dm.geo.n_pml_yp;
    }

    for (double i = dm.geo.n_pml_yn - 0.25; i >= 0.75; i--) {
        rho_m_yn(dm.geo.n_pml_yn - static_cast<int>(i) - 1) =
            i / dm.geo.n_pml_yn;
    }

    for (double i = 0.75; i <= dm.geo.n_pml_yp - 0.25; i++) {
        rho_m_yp(static_cast<int>(i)) = i / dm.geo.n_pml_yp;
    }

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

void Solver::StepForward() {
    // 更新源
    UpdateSource();

    // TMz模式计算
    TMz();

    // TEz模式计算
    TEz();
}

void Solver::CalculateEx() { dm.Ex = dm.Dx / Const::EPS_0; }

void Solver::CalculateEy() { dm.Ey = dm.Dy / Const::EPS_0; }

void Solver::CalculateEz() { dm.Ez = dm.Dz / Const::EPS_0; }

void Solver::CalculateHx() { dm.Hx = dm.Bx / Const::MU_0; }

void Solver::CalculateHy() { dm.Hy = dm.By / Const::MU_0; }

void Solver::CalculateHz() { dm.Hz = dm.Bz / Const::MU_0; }

void Solver::UpdateBx() {
// Bx
// -y
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 0; j < dm.geo.n_pml_yn - 1; ++j) {
            Phi_mx_yn(i, j) = b_mx_yn(i, j) * Phi_mx_yn(i, j) +
                              a_mx_yn(i, j) * (dm.Ez(i, j + 1) - dm.Ez(i, j));
            dm.Bx(i, j) -= C_Bx_dyn(i, j) * (dm.Ez(i, j + 1) - dm.Ez(i, j)) +
                           CB_dt * Phi_mx_yn(i, j);
        }
    }
    spdlog::debug("bx -y ok");

    // 正y轴方向的CPML区域中的Dx更新程序
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 0; j < dm.geo.n_pml_yp - 1; ++j) {
            int tmp = dm.geo.ni - dm.geo.n_pml_yp + j;

            Phi_mx_yp(i, j) =
                b_mx_yp(i, j) * Phi_mx_yp(i, j) +
                a_mx_yp(i, j) * (dm.Ez(i, tmp + 1) - dm.Ez(i, tmp));

            dm.Bx(i, j) -=
                C_Bx_dyp(i, j) * (dm.Ez(i, tmp + 1) - dm.Ez(i, tmp)) +
                CB_dt * Phi_mx_yp(i, j);
        }
    }
    spdlog::debug("bx +y ok");

    // 非PML域更新Dx
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = dm.geo.n_pml_yn; j < dm.geo.ni - dm.geo.n_pml_yp; ++j) {
            dm.Bx(i, j) -= CB_dt_dy * (dm.Ez(i, j + 1) - dm.Ez(i, j));
        }
    }
    spdlog::debug("bx y ok");
}

void Solver::UpdateBy() {
// By
// 负x轴方向PML
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.n_pml_xn - 1; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            Phi_my_xn(i, j) = b_my_xn(i, j) * Phi_my_xn(i, j) +
                              a_my_xn(i, j) * (dm.Ez(i + 1, j) - dm.Ez(i, j));

            dm.By(i, j) += C_By_dxn(i, j) * (dm.Ez(i + 1, j) - dm.Ez(i, j)) +
                           CB_dt * Phi_my_xn(i, j);
        }
    }
    spdlog::debug("by -x ok");

    // 正x轴方向PML
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.n_pml_xp - 1; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            int tmp = dm.geo.nj - dm.geo.n_pml_xp + i;

            Phi_my_xp(i, j) =
                b_my_xp(i, j) * Phi_my_xp(i, j) +
                a_my_xp(i, j) * (dm.Ez(tmp + 1, j) - dm.Ez(tmp, j));

            dm.By(i, j) +=
                C_By_dxp(i, j) * (dm.Ez(tmp + 1, j) - dm.Ez(tmp, j)) +
                CB_dt * Phi_my_xp(i, j);
        }
    }

    // 非PML
#pragma omp parallel for collapse(2)
    for (int i = dm.geo.n_pml_xn; i < dm.geo.nj - dm.geo.n_pml_xp; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            dm.By(i, j) += CB_dt_dx * (dm.Ez(i + 1, j) - dm.Ez(i, j));
        }
    }
}

void Solver::UpdateBz() {
// Bz更新公式中Δx部分负x轴方向(xn)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.n_pml_xn - 1; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            Phi_mz_xn(i, j) = b_mz_xn(i, j) * Phi_mz_xn(i, j) +
                              a_mz_xn(i, j) * (dm.Ex(i + 1, j) - dm.Ex(i, j));

            dm.Bz(i, j) += C_Bz_dxn(i, j) * (dm.Ex(i + 1, j) - dm.Hx(i, j)) +
                           CB_dt * (Phi_mz_xn(i, j) - dm.Mz(i, j));
        }
    }

// Bz更新公式中Δx部分正x轴方向(xp)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.n_pml_xp - 1; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            int tmp = dm.geo.nj - dm.geo.n_pml_xp + i;
            Phi_mz_xp(i, j) =
                b_mz_xp(i, j) * Phi_mz_xp(i, j) +
                a_mz_xp(i, j) * (dm.Ex(tmp + 1, j) - dm.Ex(tmp, j));

            dm.Bz(i, j) +=
                C_Bz_dxp(i, j) * (dm.Ex(tmp + 1, j) - dm.Ex(tmp, j)) -
                CB_dt * (Phi_mz_xp(i, j) - dm.Mz(i, j));
        }
    }

// Bz更新公式中Δx部分非PML计算域更新
#pragma omp parallel for collapse(2)
    for (int i = dm.geo.n_pml_xn; i < dm.geo.nj - dm.geo.n_pml_xp; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            dm.Bz(i, j) += CB_dt_dx * (dm.Ex(i + 1, j) - dm.Ex(i, j)) -
                           CB_dt * dm.Mz(i, j);
        }
    }

// Bz更新公式中Δy部分负y轴方向(yn)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 0; j < dm.geo.n_pml_yn - 1; ++j) {
            Phi_mz_yn(i, j) = b_mz_yn(i, j) * Phi_mz_yn(i, j) +
                              a_mz_yn(i, j) * (dm.Ey(i, j + 1) - dm.Ey(i, j));

            dm.Bz(i, j) -= C_Bz_dyn(i, j) * (dm.Ey(i, j + 1) - dm.Ey(i, j)) +
                           CB_dt * Phi_mz_yn(i, j);
        }
    }

// Bz更新公式中Δy部分正y轴方向(yp)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 0; j < dm.geo.n_pml_yp - 1; ++j) {
            int tmp = dm.geo.ni - dm.geo.n_pml_yp + j;

            Phi_mz_yp(i, j) =
                b_mz_yp(i, j) * Phi_mz_yp(i, j) +
                a_mz_yp(i, j) * (dm.Ey(i, tmp + 1) - dm.Ey(i, tmp));

            dm.Bz(i, j) -=
                C_Bz_dyp(i, j) * (dm.Ey(i, tmp + 1) - dm.Ey(i, tmp)) +
                CB_dt * Phi_mz_yp(i, j);
        }
    }

// Bz更新公式中Δy部分非PML计算域更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = dm.geo.n_pml_yn; j < dm.geo.ni - dm.geo.n_pml_yp; ++j) {
            dm.Bz(i, j) -= CB_dt_dy * (dm.Ey(i, j + 1) - dm.Ey(i, j));
        }
    }
}

void Solver::UpdateDx() {
// Dx更新公式中Δy部分负y轴(yn)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 1; j < dm.geo.n_pml_yn; ++j) {
            Phi_ex_yn(i, j) = b_ex_yn(i, j) * Phi_ex_yn(i, j) +
                              a_ex_yn(i, j) * (dm.Hz(i, j) - dm.Hz(i, j - 1));
            dm.Dx(i, j) -= C_Dx_dyn(i, j) * (dm.Hz(i, j + 1) - dm.Hz(i, j)) +
                           CD_dt * Phi_ex_yn(i, j);
        }
    }
    spdlog::debug("bx -y ok");

// Dx更新公式中Δy部分负正y轴(yp)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 1; j < dm.geo.n_pml_yp; ++j) {
            int tmp = dm.geo.ni - dm.geo.n_pml_yp + j;

            Phi_ex_yp(i, j) =
                b_ex_yp(i, j) * Phi_ex_yp(i, j) +
                a_ex_yp(i, j) * (dm.Hz(i, tmp) - dm.Hz(i, tmp - 1));

            dm.Dx(i, j) -=
                C_Dx_dyp(i, j) * (dm.Hz(i, tmp) - dm.Hz(i, tmp - 1)) +
                CD_dt * Phi_ex_yp(i, j);
        }
    }

// Dx更新公式中Δy部分非PML更新
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = dm.geo.n_pml_yn; j < dm.geo.ni - dm.geo.n_pml_yp; ++j) {
            dm.Dx(i, j) -= CD_dt_dy * (dm.Hz(i, j) - dm.Hz(i, j - 1));
        }
    }
}

void Solver::UpdateDy() {
// Dy更新公式中Δx部分负x轴(xn)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 1; i < dm.geo.n_pml_xn; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            Phi_ey_xn(i, j) = b_ey_xn(i, j) * Phi_ey_xn(i, j) +
                              a_ey_xn(i, j) * (dm.Hz(i, j) - dm.Hz(i - 1, j));

            dm.Dy(i, j) += C_Dy_dxn(i, j) * (dm.Hz(i, j) - dm.Hz(i - 1, j)) +
                           CD_dt * Phi_ey_xn(i, j);
        }
    }

// Dy更新公式中Δx部分正x轴(xp)PML更新
#pragma omp parallel for collapse(2)
    for (int i = 1; i < dm.geo.n_pml_xp; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            int tmp = dm.geo.nj - dm.geo.n_pml_xp + i;

            Phi_ey_xp(i, j) =
                b_ey_xp(i, j) * Phi_ey_xp(i, j) +
                a_ey_xp(i, j) * (dm.Hz(tmp, j) - dm.Hz(tmp - 1, j));

            dm.Dy(i, j) +=
                C_Dy_dxp(i, j) * (dm.Hz(tmp, j) - dm.Hz(tmp - 1, j)) +
                CD_dt * Phi_ey_xp(i, j);
        }
    }

// Dy更新公式中Δx部分非PML更新
#pragma omp parallel for collapse(2)
    for (int i = dm.geo.n_pml_xn; i < dm.geo.nj - dm.geo.n_pml_xp; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            dm.Dy(i, j) += CD_dt_dx * (dm.Hz(i, j) - dm.Hz(i - 1, j));
        }
    }
}

void Solver::UpdateDz() {
// Dz
// x
// 负PML
#pragma omp parallel for collapse(2)
    for (int i = 1; i < dm.geo.n_pml_xn; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            Phi_ez_xn(i, j) = b_ez_xn(i, j) * Phi_ez_xn(i, j) +
                              a_ez_xn(i, j) * (dm.Hy(i, j) - dm.Hy(i - 1, j));

            dm.Dz(i, j) += C_Dz_dxn(i, j) * (dm.Hy(i, j) - dm.Hy(i - 1, j)) +
                           CD_dt * (Phi_ez_xn(i, j) - dm.Jz(i, j));
        }
    }

    // 正PML
#pragma omp parallel for collapse(2)
    for (int i = 1; i < dm.geo.n_pml_xp; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            int tmp = dm.geo.nj - dm.geo.n_pml_xp + i;
            Phi_ez_xp(i, j) =
                b_ez_xp(i, j) * Phi_ez_xp(i, j) +
                a_ez_xp(i, j) * (dm.Hy(tmp, j) - dm.Hy(tmp - 1, j));

            dm.Dz(i, j) +=
                C_Dz_dxp(i, j) * (dm.Hy(tmp, j) - dm.Hy(tmp - 1, j)) -
                CD_dt * (Phi_ez_xp(i, j) - dm.Jz(i, j));
        }
    }

    // 非PML
#pragma omp parallel for collapse(2)
    for (int i = dm.geo.n_pml_xn; i < dm.geo.nj - dm.geo.n_pml_xp; ++i) {
        for (int j = 0; j < dm.geo.ni; ++j) {
            dm.Dz(i, j) += CD_dt_dx * (dm.Hy(i, j) - dm.Hy(i - 1, j)) -
                           CD_dt * dm.Jz(i, j);
        }
    }

    // y
    // 负PML
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 1; j < dm.geo.n_pml_yn; ++j) {
            Phi_ez_yn(i, j) = b_ez_yn(i, j) * Phi_ez_yn(i, j) +
                              a_ez_yn(i, j) * (dm.Hx(i, j) - dm.Hx(i, j - 1));

            dm.Dz(i, j) -= C_Dz_dyn(i, j) * (dm.Hx(i, j) - dm.Hx(i, j - 1)) +
                           CD_dt * Phi_ez_yn(i, j);
        }
    }

    // 正PML
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = 1; j < dm.geo.n_pml_yp; ++j) {
            int tmp = dm.geo.ni - dm.geo.n_pml_yp + j;

            Phi_ez_yp(i, j) =
                b_ez_yp(i, j) * Phi_ez_yp(i, j) +
                a_ez_yp(i, j) * (dm.Hx(i, tmp) - dm.Hx(i, tmp - 1));

            dm.Dz(i, j) -=
                C_Dz_dyp(i, j) * (dm.Hx(i, tmp) - dm.Hx(i, tmp - 1)) +
                CD_dt * Phi_ez_yp(i, j);
        }
    }

    // 非PML
#pragma omp parallel for collapse(2)
    for (int i = 0; i < dm.geo.nj; ++i) {
        for (int j = dm.geo.n_pml_yn; j < dm.geo.ni - dm.geo.n_pml_yp; ++j) {
            dm.Dz(i, j) -= CD_dt_dy * (dm.Hx(i, j) - dm.Hx(i, j - 1));
        }
    }
}

void Solver::UpdateSource() {
    // double sin_current = 1600;
    // double cos_current = 1600;
    double sin_current = dm.I * sin(2 * Const::PI * dm.f * dm.time);
    double cos_current = dm.I * cos(2 * Const::PI * dm.f * dm.time);

    double cen = static_cast<int>(((dm.geo.ni + 1) / 2) - 1);

    auto sr1 = Eigen::seq(cen + 23, cen + 25);
    auto sr2 = Eigen::seq(cen - 25, cen - 23);
    auto sc1 = Eigen::seq(cen + 18, cen + 20);
    auto sc2 = Eigen::seq(cen - 20, cen - 18);

    auto cr1 = Eigen::seq(cen + 18, cen + 20);
    auto cr2 = Eigen::seq(cen - 20, cen - 18);
    auto cc1 = Eigen::seq(cen + 23, cen + 25);
    auto cc2 = Eigen::seq(cen - 25, cen - 23);

    // +sin
    dm.Jz(sr1, sc1) = sin_current;
    dm.Jz(sr1, sc2) = sin_current;

    // -sin
    dm.Jz(sr2, sc1) = -sin_current;
    dm.Jz(sr2, sc2) = -sin_current;

    // +cos
    dm.Jz(cr1, cc2) = cos_current;
    dm.Jz(cr2, cc2) = cos_current;

    // -cos
    dm.Jz(cr1, cc1) = -cos_current;
    dm.Jz(cr2, cc1) = -cos_current;
}

void Solver::TMz() {
    UpdateDz();
    CalculateEz();
    UpdateBx();
    UpdateBy();
    CalculateHx();
    CalculateHy();
};

void Solver::TEz() {
    UpdateBz();
    CalculateHz();
    UpdateDx();
    UpdateDy();
    CalculateEx();
    CalculateEy();
}