#include "solver.h"

#include <math.h>

#include <iostream>

#include "const.h"
#include "field.h"

Solver::Solver(Domain &domain, int max_it, double tol)
    : domain(domain), max_solver_it(max_it), tolerance(tol) {
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

    for (int i = 0; i < domain.geo.n_pml_xn; i++) {
        sigma_e_xn[i] = sigma_max_xn * pow(rho_e_xn[i], m_pml);
        sigma_m_xn[i] =
            Const::MU_0 / Const::EPS_0 * sigma_max_xn * pow(rho_m_xn[i], m_pml);
    }

    for (int i = 0; i < domain.geo.n_pml_xp; i++) {
        sigma_e_xp[i] = sigma_max_xp * pow(rho_e_xp[i], m_pml);
        sigma_m_xp[i] =
            Const::MU_0 / Const::EPS_0 * sigma_max_xp * pow(rho_m_xp[i], m_pml);
    }

    for (int i = 0; i < domain.geo.n_pml_yn; i++) {
        sigma_e_yn[i] = sigma_max_yn * pow(rho_e_yn[i], m_pml);
        sigma_m_yn[i] =
            Const::MU_0 / Const::EPS_0 * sigma_max_yn * pow(rho_m_yn[i], m_pml);
    }

    for (int i = 0; i < domain.geo.n_pml_yp; i++) {
        sigma_e_yp[i] = sigma_max_yp * pow(rho_e_yp[i], m_pml);
        sigma_m_yp[i] =
            Const::MU_0 / Const::EPS_0 * sigma_max_yp * pow(rho_m_yp[i], m_pml);
    }

    // 设置kappa参数
    for (int i = 0; i < domain.geo.n_pml_xn; i++) {
        kappa_e_xn[i] = 1 + (kappa_max_xn - 1) * pow(rho_e_xn[i], m_pml);
        kappa_m_xn[i] = 1 + (kappa_max_xn - 1) * pow(rho_m_xn[i], m_pml);
    }

    for (int i = 0; i < domain.geo.n_pml_xp; i++) {
        kappa_e_xp[i] = 1 + (kappa_max_xp - 1) * pow(rho_e_xp[i], m_pml);
        kappa_m_xp[i] = 1 + (kappa_max_xp - 1) * pow(rho_m_xp[i], m_pml);
    }

    for (int i = 0; i < domain.geo.n_pml_yn; i++) {
        kappa_e_yn[i] = 1 + (kappa_max_yn - 1) * pow(rho_e_yn[i], m_pml);
        kappa_m_yn[i] = 1 + (kappa_max_yn - 1) * pow(rho_m_yn[i], m_pml);
    }

    for (int i = 0; i < domain.geo.n_pml_yp; i++) {
        kappa_e_yp[i] = 1 + (kappa_max_yp - 1) * pow(rho_e_yp[i], m_pml);
        kappa_m_yp[i] = 1 + (kappa_max_yp - 1) * pow(rho_m_yp[i], m_pml);
    }

    // 设置alpha参数
    for (int i = 0; i < domain.geo.n_pml_xn; i++) {
        alpha_e_xn[i] = alpha_max * (1 - rho_e_xn[i]);
        alpha_m_xn[i] =
            Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_xn[i]);
    }

    for (int i = 0; i < domain.geo.n_pml_xp; i++) {
        alpha_e_xp[i] = alpha_max * (1 - rho_e_xp[i]);
        alpha_m_xp[i] =
            Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_xp[i]);
    }

    for (int i = 0; i < domain.geo.n_pml_yn; i++) {
        alpha_e_yn[i] = alpha_max * (1 - rho_e_yn[i]);
        alpha_m_yn[i] =
            Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_yn[i]);
    }

    for (int i = 0; i < domain.geo.n_pml_yp; i++) {
        alpha_e_yp[i] = alpha_max * (1 - rho_e_yp[i]);
        alpha_m_yp[i] =
            Const::MU_0 / Const::EPS_0 * alpha_max * (1 - rho_m_yp[i]);
    }

    // b
    for (int i = 0; i < domain.geo.n_pml_xn; i++) {
        b_e_xn_v[i] = exp(-(sigma_e_xn[i] / kappa_e_xn[i] + alpha_e_xn[i]) *
                          domain.getDt() / Const::EPS_0);
        b_m_xn_v[i] = exp(-(sigma_m_xn[i] / kappa_m_xn[i] + alpha_m_xn[i]) *
                          domain.getDt() / Const::MU_0);
    }

    for (int i = 0; i < domain.geo.n_pml_xp; i++) {
        b_e_xp_v[i] = exp(-(sigma_e_xp[i] / kappa_e_xp[i] + alpha_e_xp[i]) *
                          domain.getDt() / Const::EPS_0);
        b_m_xp_v[i] = exp(-(sigma_m_xp[i] / kappa_m_xp[i] + alpha_m_xp[i]) *
                          domain.getDt() / Const::MU_0);
    }

    for (int i = 0; i < domain.geo.n_pml_yn; i++) {
        b_e_yn_v[i] = exp(-(sigma_e_yn[i] / kappa_e_yn[i] + alpha_e_yn[i]) *
                          domain.getDt() / Const::EPS_0);
        b_m_yn_v[i] = exp(-(sigma_m_yn[i] / kappa_m_yn[i] + alpha_m_yn[i]) *
                          domain.getDt() / Const::MU_0);
    }

    for (int i = 0; i < domain.geo.n_pml_yp; i++) {
        b_e_yp_v[i] = exp(-(sigma_e_yp[i] / kappa_e_yp[i] + alpha_e_yp[i]) *
                          domain.getDt() / Const::EPS_0);
        b_m_yp_v[i] = exp(-(sigma_m_yp[i] / kappa_m_yp[i] + alpha_m_yp[i]) *
                          domain.getDt() / Const::MU_0);
    }

    // a
    for (int i = 0; i < domain.geo.n_pml_xn; i++) {
        a_e_xn_v[i] = sigma_e_xn[i] * (b_e_xn_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_e_xn[i] *
                       (sigma_e_xn[i] + alpha_e_xn[i] * kappa_e_xn[i]));
        a_m_xn_v[i] = sigma_m_xn[i] * (b_m_xn_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_m_xn[i] *
                       (sigma_m_xn[i] + alpha_m_xn[i] * kappa_m_xn[i]));
    }

    for (int i = 0; i < domain.geo.n_pml_xp; i++) {
        a_e_xp_v[i] = sigma_e_xp[i] * (b_e_xp_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_e_xp[i] *
                       (sigma_e_xp[i] + alpha_e_xp[i] * kappa_e_xp[i]));
        a_m_xp_v[i] = sigma_m_xp[i] * (b_m_xp_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_m_xp[i] *
                       (sigma_m_xp[i] + alpha_m_xp[i] * kappa_m_xp[i]));
    }

    for (int i = 0; i < domain.geo.n_pml_yn; i++) {
        a_e_yn_v[i] = sigma_e_yn[i] * (b_e_yn_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_e_yn[i] *
                       (sigma_e_yn[i] + alpha_e_yn[i] * kappa_e_yn[i]));
        a_m_yn_v[i] = sigma_m_yn[i] * (b_m_yn_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_m_yn[i] *
                       (sigma_m_yn[i] + alpha_m_yn[i] * kappa_m_yn[i]));
    }

    for (int i = 0; i < domain.geo.n_pml_yp; i++) {
        a_e_yp_v[i] = sigma_e_yp[i] * (b_e_yp_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_e_yp[i] *
                       (sigma_e_yp[i] + alpha_e_yp[i] * kappa_e_yp[i]));
        a_m_yp_v[i] = sigma_m_yp[i] * (b_m_yp_v[i] - 1) /
                      (domain.geo.dh[0] * kappa_m_yp[i] *
                       (sigma_m_yp[i] + alpha_m_yp[i] * kappa_m_yp[i]));
    }
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