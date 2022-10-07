#pragma once

#include <chrono>
#include <fstream>
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
    void OutputDat();
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
void Output<Mixin>::OutputDat() {
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
           "\"Y\",\"Bx\",\"By\",\"Bz\",\"Ex\",\"Ey\",\"Ez\",\"Jx\",\"Jy\","
           "\"Jz\"\n";
    out << "ZONE i=" << dm.geo.ni << " j=" << dm.geo.nj
        << " SOLUTIONTIME=" << dm.time << "\n";
    for (int i = 0; i < dm.geo.ni; i++) {
        for (int j = 0; j < dm.geo.nj; j++) {
            out << i << " ";
            out << j << " ";
            out << dm.Bx(i, j) << " ";
            out << dm.By(i, j) << " ";
            out << dm.Bz(i, j) << " ";
            out << dm.Ex(i, j) << " ";
            out << dm.Ey(i, j) << " ";
            out << dm.Ez(i, j) << " ";
            out << dm.Jx(i, j) << " ";
            out << dm.Jy(i, j) << " ";
            out << dm.Jz(i, j) << "\n";
        }
    }

    // out << "<?xml version=\"1.0\"?>\n";
    // out << "<VTKFile type=\"ImageData\">\n";
    // Vec2d origin = dm.geo.x0;
    // Vec2d delta = dm.geo.dh;
    // out << "<ImageData WholeExtent=\"0 " << dm.geo.ni - 1 << " 0 "
    //     << dm.geo.nj - 1 << "\"\n";
    // out << "Origin=\"" << origin[0] << " " << origin[1] << "\" ";
    // out << "Spacing=\"" << delta[0] << " " << delta[1] << "\">\n";

    // /*output data stored on nodes (point data)*/
    // out << "<PointData>\n";

    // // // 单元面积
    // // out << "<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" "
    // //        "format=\"ascii\" type=\"Float64\">\n";
    // // out << dm.geo.node_area;
    // // out << "</DataArray>\n";

    // // 电场信息
    // out << "<DataArray Name=\"E\" NumberOfComponents=\"3\" format=\"ascii\" "
    //        "type=\"Float64\">\n";
    // for (int i = 0; i < dm.geo.ni; ++i) {
    //     for (int j = 0; j < dm.geo.nj; ++j) {
    //         out << dm.Dx(i, j) << " " << dm.Dy(i, j) << " " << dm.Dz(i, j)
    //             << "\n";
    //     }
    // }
    // out << "</DataArray>\n";

    // // 磁场信息
    // out << "<DataArray Name=\"B\" NumberOfComponents=\"3\" format=\"ascii\" "
    //        "type=\"Float64\">\n";
    // for (int i = 0; i < dm.geo.ni; ++i) {
    //     for (int j = 0; j < dm.geo.nj; ++j) {
    //         out << dm.Hx(i, j) << " " << dm.Hy(i, j) << " " << dm.Hz(i, j)
    //             << "\n";
    //     }
    // }
    // out << "</DataArray>\n";

    // // 电流密度信息
    // out << "<DataArray Name=\"J\" NumberOfComponents=\"1\" format=\"ascii\" "
    //        "type=\"Float64\">\n";
    // for (int i = 0; i < dm.geo.ni; ++i) {
    //     for (int j = 0; j < dm.geo.nj; ++j) {
    //         out << dm.Jx(i, j) << " " << dm.Jy(i, j) << " " << dm.Jz(i, j)
    //             << "\n";
    //     }
    // }
    // out << "</DataArray>\n";

    // /*close out tags*/
    // out << "</PointData>\n";
    // out << "</ImageData>\n";
    // out << "</VTKFile>\n";

    out.close();
}

template <typename Mixin>
void Output<Mixin>::Launch() {
    std::cout << "Calculating..."
              << "\n";
    std::cout << "Frequency: " << dm.f / 1000 << " [kHz]"
              << "\n";
    std::cout << "Current density: " << dm.I << " [A/m^2]"
              << "\n";
    std::cout << "Duartaion of a time step: " << dm.dt << " [s]"
              << "\n\n";

    while (dm.advanceTime()) {
        std::cout << "\rCurrent step / Total step: " << dm.ts << "/"
                  << dm.num_ts << std::flush;
        Mixin::StepForward();

        if (dm.ts % 1000 == 0) {
            OutputDat();
        }
    }

    std::cout << "\nCalculate complete\n";
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
