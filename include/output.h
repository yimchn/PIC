#pragma once

#include <chrono>
#include <fstream>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

#include "domain.h"
#include "species.h"

template <typename Mixin>
class Output : public Mixin {
   public:
    Domain& dm;

    Output(Domain& dm, double it, double tol);
    // void fields(Domain& dm, std::vector<Species>& species);
    void OutputFields();
    void DiagB();
    void DiagE();
    void DiagJ();
    void OutputDiag(Domain& dm, std::vector<Species>& species);
    void Launch();
};

template <typename Mixin>
Output<Mixin>::Output(Domain& dm, double it, double tol)
    : Mixin(dm, it, tol), dm(dm) {}

template <typename Mixin>
void Output<Mixin>::OutputFields() {
    std::stringstream name;
    name << "../results/fields_" << std::setfill('0') << std::setw(10)
         << dm.getTs() << ".dat";

    /*open output file*/
    std::ofstream out(name.str());
    if (!out.is_open()) {
        std::cerr << "Could not open " << name.str() << std::endl;
        return;
    }

    out << "TITLE = Magnetic Field 2d\n";
    out << "VARIABLES = \"X\", "
           "\"Y\",\"Hx\",\"Hy\",\"Hz\",\"Ex\",\"Ey\",\"Ez\",\"Jx\",\"Jy\","
           "\"Jz\"\n";
    out << "ZONE i=" << dm.geo.ni << " j=" << dm.geo.nj
        << " SOLUTIONTIME=" << dm.time << "\n";
    for (int i = 0; i < dm.geo.ni; i++) {
        for (int j = 0; j < dm.geo.nj; j++) {
            out << i << " ";
            out << j << " ";
            out << dm.Hx(i, j) << " ";
            out << dm.Hy(i, j) << " ";
            out << dm.Hz(i, j) << " ";
            out << dm.Dx(i, j) << " ";
            out << dm.Dy(i, j) << " ";
            out << dm.Dz(i, j) << " ";
            out << dm.Jx(i, j) << " ";
            out << dm.Jy(i, j) << " ";
            out << dm.Jz(i, j) << "\n";
        }
    }

    out.close();
}

template <typename Mixin>
void Output<Mixin>::Launch() {
    using namespace indicators;

    // BlockProgressBar bar{
    //     option::BarWidth{80}, option::ForegroundColor{Color::white},
    //     option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
    //     option::MaxProgress{dm.num_ts}};
    ProgressBar bar{
        option::BarWidth{80},
        option::Start{"["},
        option::Fill{"="},
        option::Lead{">"},
        option::Remainder{" "},
        option::End{" ]"},
        option::MaxProgress{dm.num_ts},
        // option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
    };

    // 进度条上的文字
    std::cout << "Calculating..."
              << "\n";
    std::cout << "Total number of iterations: " << dm.num_ts << "\n";
    std::cout << "Duartaion of a time step:" << dm.dt << "\n";

    show_console_cursor(false);
    while (dm.advanceTime()) {
        Mixin::StepForward();

        if (dm.ts % 1000 == 0) {
            OutputFields();
        }

        // 在进度条最后显示迭代次数
        bar.set_option(option::PostfixText{std::to_string(dm.ts) + "/" +
                                           std::to_string(dm.num_ts)});

        // 更新进度条
        bar.tick();
    }

    // bar.mark_as_completed();
    // show_console_cursor(true);
}

/*save runtime diagnostics to a file*/
template <typename Mixin>
void Output<Mixin>::OutputDiag(Domain& domain, std::vector<Species>& species) {
    std::ofstream f_diag;

    // is the file open?
    if (!f_diag.is_open()) {
        f_diag.open("../diag/runtime_diags.csv");
        f_diag << "ts,time,wall_time";
        for (Species& sp : species)
            f_diag << ",mp_count." << sp.name << ",real_count." << sp.name
                   << ",px." << sp.name << ",py." << sp.name << ",pz."
                   << sp.name << ",KE." << sp.name;
        f_diag << ",PE,E_total" << std::endl;
    }

    f_diag << domain.getTs() << "," << domain.getTime();
    f_diag << "," << domain.getWallTime();

    double tot_KE = 0;
    for (Species& sp : species) {
        double KE = sp.getKE();  // species kinetic energy
        tot_KE += KE;            // increment total energy
        Vec2d mom = sp.getMomentum();

        f_diag << "," << sp.getNp() << "," << sp.getRealCount() << "," << mom[0]
               << "," << mom[1] << "," << mom[2] << "," << KE;
    }

    // write out system potential and total energy
    double PE = domain.getPE();
    f_diag << "," << PE << "," << (tot_KE + PE);

    f_diag << "\n";  // use \n to avoid flush to disc
    if (domain.getTs() % 25 == 0) f_diag.flush();
}

template <typename Mixin>
void Output<Mixin>::DiagB() {
    std::stringstream name;
    name << "../diag/B/B_" << std::setfill('0') << std::setw(10) << dm.getTs()
         << ".dat";

    /*open output file*/
    std::ofstream out(name.str());
    if (!out.is_open()) {
        std::cerr << "Could not open " << name.str() << std::endl;
        return;
    }

    out << sqrt(pow(dm.Hx, 2) + pow(dm.Hy, 2));

    out.close();
}

template <typename Mixin>
void Output<Mixin>::DiagE() {
    std::stringstream name;
    name << "../diag/E/E_" << std::setfill('0') << std::setw(10) << dm.getTs()
         << ".dat";

    /*open output file*/
    std::ofstream out(name.str());
    if (!out.is_open()) {
        std::cerr << "Could not open " << name.str() << std::endl;
        return;
    }

    out << dm.Dz;

    out.close();
}

template <typename Mixin>
void Output<Mixin>::DiagJ() {
    std::stringstream name;
    name << "../diag/J/J_" << std::setfill('0') << std::setw(10) << dm.getTs()
         << ".dat";

    /*open output file*/
    std::ofstream out(name.str());
    if (!out.is_open()) {
        std::cerr << "Could not open " << name.str() << std::endl;
        return;
    }

    out << dm.Jz;

    out.close();
}
