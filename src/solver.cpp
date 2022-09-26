#include "solver.h"

#include <math.h>

#include <iostream>

#include "const.h"
#include "field.h"

Solver::Solver(Domain &dm, int max_it, double tol)
        : dm(dm), max_solver_it(max_it), tolerance(tol) {
    Init();
    spdlog::debug("Initial success");
    UpdatePml();
    spdlog::debug("All of the parameters have been allocated");
}

void Solver::Init() {
    dt = dm.getDt();
    delta = dm.geo.dh[0];
    eps_r = 1;
    mu_r = 1;
    sigma = 0;
    gradient_k = 1;
    R0 = 1.0e-7;
    m = 3;
    eta = sqrt(Const::MU_0 / Const::EPS_0);

    double tmp1 = Const::EPS_0 * eps_r / dt;
    double tmp2 = sigma / 2;
    c1 = dt / (Const::EPS_0 * eps_r * delta);
    c2 = (tmp1 - tmp2) / (tmp1 + tmp2);
    c3 = 1 / (tmp1 + tmp2);
    c4 = tmp1;
    c5 = -tmp2;

    d1 = 1;
    d2 = dt / delta;
    d3 = 1;
    d4 = 1 / Const::MU_0;

    c1x = c1;
    c2x = c2;
    c3x = c3;
    c4x = c4;
    c5x = c5;

    c1y = c1;
    c2y = c2;
    c3y = c3;
    c4y = c4;
    c5y = c5;

    d1z = d1;
    d2z = d2;
    d3z = d3;
    d4z = d4;

    d1x = dt / (Const::MU_0 * mu_r * delta);
    d2x = c2;
    d3x = c3;
    d4x = c4;
    d5x = c5;

    d1y = -dt / (Const::MU_0 * mu_r * delta);
    d2x = c2;
    d3x = c3;
    d4x = c4;
    d5x = c5;

    c1z = d1;
    c2z = c2;
    c3z = c3;
    c4z = 1 / Const::EPS_0;
}

void Solver::UpdatePml() {
    boundary_width = static_cast<double>(dm.geo.npml * delta);
    sigma_max = -log(R0) * (m + 1) / (2 * eta * boundary_width);
    boundary_factor = sigma_max / (delta * (pow(boundary_width, m) * (m + 1)));
    gradient_k = 1;
    double gradient_conductivity = 0;
    // double x1 = 0;
    // double x2 = 0;

    // #pragma omp parallel for
    for (int i = 0; i < dm.geo.npml; i++) {
        double x = 0.0;
        double x1 = (x + 0.5) * delta;
        double x2 = (x - 0.5) * delta;

        if (i == 0) {
            gradient_conductivity = boundary_factor * (pow(x1, m + 1));
        } else {
            gradient_conductivity =
                    boundary_factor * (pow(x1, m + 1) - pow(x2, m + 1));
        }

        double tmp1 = gradient_k * Const::EPS_0 / dt;
        double tmp2 = gradient_conductivity / 2.0;

        gradientC2(i) = (tmp1 - tmp2) / (tmp1 + tmp2);
        gradientC3(i) = 1 / (tmp1 + tmp2);

        x1 = (x + 1) * delta;
        x2 = (x + 0) * delta;

        gradient_conductivity =
                boundary_factor * (pow(x1, m + 1) - pow(x2, m + 1));

        tmp1 = gradient_k * Const::EPS_0 / dt;
        tmp2 = gradient_conductivity / 2;

        gradientC4(i) = (tmp1 + tmp2);
        gradientC5(i) = (tmp1 - tmp2);

        tmp1 = gradient_k / dt;
        tmp2 = gradient_conductivity / (2 * Const::EPS_0);

        gradientD1(i) = (tmp1 - tmp2) / (tmp1 + tmp2);
        gradientD2(i) = (1 / delta) / (tmp1 + tmp2);

        tmp1 *= Const::MU_0;
        tmp2 *= Const::MU_0;

        gradientD3(i) = (tmp1 - tmp2) / (tmp1 + tmp2);
        gradientD4(i) = (1 / dt) / (tmp1 + tmp2);
        ++x;
    }

    // ex, ey, ez ------ front/back (y)
#pragma omp parallel for
    for (int j = 0; j < dm.geo.npml; j++) {
        for (int i = 0; i < dm.geo.ni; i++) {
            c2x(i, dm.geo.npml - j) = gradientC2(j);
            c3x(i, dm.geo.npml - j) = gradientC3(j);
            c2x(i, dm.geo.nj - dm.geo.npml + j) = gradientC2(j);
            c3x(i, dm.geo.nj - dm.geo.npml + j) = gradientC3(j);
        }

        c4y(dm.geo.npml - j - 1) = gradientC4(j);
        c5y(dm.geo.npml - j - 1) = gradientC5(j);
        c4y(dm.geo.nj - dm.geo.npml + j) = gradientC4(j);
        c5y(dm.geo.nj - dm.geo.npml + j) = gradientC5(j);

        d3z(dm.geo.npml - j - 1) = gradientD3(j);
        d4z(dm.geo.npml - j - 1) = gradientD4(j);
        d3z(dm.geo.nj - dm.geo.npml + j) = gradientD3(j);
        d4z(dm.geo.nj - dm.geo.npml + j) = gradientD4(j);
    }

    // ey, hz ----- left/right (x)
    for (int i = 0; i < dm.geo.npml; i++) {
        for (int j = 0; j < dm.geo.nj; j++) {
            c2y(dm.geo.npml - i, j) = gradientC2(i);
            c3y(dm.geo.npml - i, j) = gradientC3(i);
            c2y(dm.geo.ni - dm.geo.npml + i, j) = gradientC2(i);
            c3y(dm.geo.ni - dm.geo.npml + i, j) = gradientC3(i);
        }

        c4x(dm.geo.npml - i - 1) = gradientC4(i);
        c5x(dm.geo.npml - i - 1) = gradientC5(i);
        c4x(dm.geo.ni - dm.geo.npml + 1) = gradientC4(i);
        c5x(dm.geo.ni - dm.geo.npml + 1) = gradientC5(i);

        d1z(dm.geo.npml - i - 1) = gradientD1(i);
        d2z(dm.geo.npml - i - 1) = gradientD2(i);
        d1z(dm.geo.ni - dm.geo.npml + i) = gradientD1(i);
        d2z(dm.geo.ni - dm.geo.npml + i) = gradientD2(i);
    }
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

void Solver::EvaluateFdtd() {
    // Update magnetic fields (Hx and Hy)
#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni; ++i) {
        for (int j = 1; j < dm.geo.nj; ++j) {
            double tmp = dm.Bx(i, j);
            dm.Bx(i, j) += c1x * (dm.Ez(i, j) - dm.Ez(i, j - 1));
            dm.Hx(i, j) = c2x(i, j) * dm.Hx(i, j) +
                          c3x(i, j) * (c4x(i) * dm.Bx(i, j) + c5x(i) * tmp);
        }
    }
    spdlog::debug("node 1 ok");

#pragma omp parallel for
    for (int i = 1; i < dm.geo.ni; ++i) {
        for (int j = 0; j < dm.geo.nj; ++j) {
            double tmp = dm.By(i, j);
            dm.By(i, j) += c1y * (dm.Ez(i - 1, j) - dm.Ez(i, j));
            dm.Hy(i, j) = c2y(i, j) * dm.Hy(i, j) +
                          c3y(i, j) * (c4y(i) * dm.By(i, j) + c5y(i) * tmp);
        }
    }
    spdlog::debug("node 2 ok");

    // Update electric fields Ez
#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni; ++i) {
        for (int j = 0; j < dm.geo.nj; ++j) {
            double tmp = dm.Dz(i, j);
            double tmp2 = dm.Jz(i, j);
            dm.Dz(i, j) =
                    d1z(i) * dm.Dz(i, j) +
                    d2z(i) * (dm.Hx(i, j + 1) - dm.Hx(i, j) + dm.Hy(i, j) -
                              dm.Hy(i + 1, j)) - d4z(j);
            dm.Ez(i, j) = d3z(j) * dm.Ez(i, j) + d4z(j) * (dm.Dz(i, j) - tmp);
        }
    }
    spdlog::debug("node 3 ok");

    // Update electirc fields (Ex and Ey)
#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni; ++i) {
        for (int j = 1; j < dm.geo.nj; ++j) {
            double tmp = dm.Dx(i, j);
            dm.Dx(i, j) += c1x * (dm.Hz(i, j) - dm.Hz(i, j - 1));
            dm.Ex(i, j) = c2x(i, j) * dm.Ex(i, j) +
                          c3x(i, j) * (c4x(i) * dm.Dx(i, j) + c5x(i) * tmp);
        }
    }
    spdlog::debug("node 4 ok");

#pragma omp parallel for
    for (int i = 1; i < dm.geo.ni; ++i) {
        for (int j = 0; j < dm.geo.nj; ++j) {
            double tmp = dm.Dy(i, j);
            dm.Dy(i, j) += c1y * (dm.Hz(i - 1, j) - dm.Hz(i, j));
            dm.Ey(i, j) = c2y(i, j) * dm.Ey(i, j) +
                          c3y(i, j) * (c4y(i) * dm.Dy(i, j) + c5y(i) * tmp);
        }
    }
    spdlog::debug("node 5 ok");

    // Update magnetic fields Hz
#pragma omp parallel for
    for (int i = 0; i < dm.geo.ni; ++i) {
        for (int j = 0; j < dm.geo.nj; ++j) {
            double tmp = dm.Bz(i, j);
            dm.Bz(i, j) =
                    d1z(i) * dm.Bz(i, j) +
                    d2z(i) * (dm.Ex(i, j + 1) - dm.Ex(i, j) + dm.Ey(i, j) -
                              dm.Ey(i + 1, j));
            dm.Hz(i, j) = d3z(j) * dm.Hz(i, j) + d4z(j) * (dm.Bz(i, j) - tmp);
        }
    }
    spdlog::debug("node 6 ok");
}

void Solver::EvaluateSource(double I, double f, double t) {
#pragma omp parallel for
    for (int i = 6 + 3; i < 20 - 3; i++) {
        for (int j = 6 + 3; j < 20 - 3; j++) {
            dm.Jz(5 + 2, j) = -I * sin(2 * Const::PI * f * t);
            dm.Jz(20 - 2, j) = I * sin(2 * Const::PI * f * t);
            dm.Jz(i, 5 + 2) = -I * cos(2 * Const::PI * f * t);
            dm.Jz(i, 20 - 2) = I * cos(2 * Const::PI * f * t);
        }
    }
}
