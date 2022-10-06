#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

#include "domain.h"
#include "species.h"

template <typename T>
class PostProcessing : public T {
   protected:
    // void fields(Domain& dm, std::vector<Species>& species);
    void OutputFields(Domain& dm);
    void DiagB(Domain& dm);
    void DiagE(Domain& dm);
    void DiagJ(Domain& dm);
    void OutputDiag(Domain& dm, std::vector<Species>& species);
    void Display(Domain dm);
};
// namespace Output {
// void fields(Domain &dm, std::vector<Species> &species);
// void fields(Domain &dm);
// void fields(Domain &dm, int step);
// void B(Domain &dm);
// void E(Domain &dm);
// void J(Domain &dm);
// void screenOutput(Domain &dm, std::vector<Species> &species);
// void diagOutput(Domain &dm, std::vector<Species> &species);
// void Display(Domain &dm, Field<Vec3d> field);
// void ProgressBar(Domain dm);
// }  // namespace Output
template <typename T>
void PostProcessing<T>::OutputFields(Domain& dm) {
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

template <typename T>
void PostProcessing<T>::Display(Domain dm) {
    //using namespace indicators;
    indicators::BlockProgressBar bar{
        option::BarWidth{80}, option::ForegroundColor{Color::white},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
        option::MaxProgress{dm.num_ts}};

    std::cout << "Total number of iterations: " << dm.num_ts << "\n";
    for (size_t i = 0; i < dm.num_ts; ++i) {
        // 执行相应的操作
        T::Solve();

        // 在进度条后面显示迭代进度
        bar.set_option(option::PostfixText {
            std::to_string(i) + "/" + std::to_string(dm.num_ts);
        });

        // 更新进度条
        bar.tick();
    }

    bar.mark_as_completed();

    indicators::show_console_cursor(true);

    // double progress =
    //	static_cast<double>(dm.ts) / static_cast<double>(dm.num_ts);
    // int bar_width = 70;

    // std::cout << "[";
    // int pos = bar_width * progress;
    // for (int i = 0; i < bar_width; ++i) {
    //	if (i < pos)
    //		std::cout << "=";
    //	else if (i == pos)
    //		std::cout << ">";
    //	else
    //		std::cout << " ";
    // }
    // std::cout << "] " << std::fixed << std::setprecision(2) << (progress *
    // 100)
    //	<< "%\r";
    // std::cout.flush();
}

/*save runtime diagnostics to a file*/
template <typename T>
void PostProcessing<T>::OutputDiag(Domain& domain,
                                   std::vector<Species>& species) {
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

template <typename T>
void PostProcessing<T>::DiagB(Domain& dm) {
    std::stringstream name;
    name << "../diag/B/B_" << std::setfill('0') << std::setw(10) << dm.getTs()
         << ".dat";

    /*open output file*/
    std::ofstream out(name.str());
    if (!out.is_open()) {
        std::cerr << "Could not open " << name.str() << std::endl;
        return;
    }

    out << std::sqrt(std::pow(dm.Hx, 2) + std::pow(dm.Hy, 2));

    out.close();
}

template <typename T>
void PostProcessing<T>::DiagE(Domain& dm) {
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

template <typename T>
void PostProcessing<T>::DiagJ(Domain& dm) {
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
