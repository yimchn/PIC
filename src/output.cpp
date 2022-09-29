#include "output.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "domain.h"
#include "species.h"

#define GET_VARIABLE_NAME(Variable) (#Variable)

/*saves fields in VTK format*/
// void Output::fields(Domain &dm, std::vector<Species> &species) {
//     std::stringstream name;
//     name << "../results/fields_" << std::setfill('0') << std::setw(5)
//          << dm.getTs() << ".vti";

//     /*open output file*/
//     std::ofstream out(name.str());
//     if (!out.is_open()) {
//         std::cerr << "Could not open " << name.str() << std::endl;
//         return;
//     }

//     /*ImageData is vtk format for structured Cartesian meshes*/
//     out << "<VTKFile type=\"ImageData\">\n";
//     Vec2d x0 = dm.getX0();
//     Vec2d dh = dm.getDh();
//     out << "<ImageData Origin=\"" << x0[0] << " " << x0[1] << "\" ";
//     out << "Spacing=\"" << dh[0] << " " << dh[1] << "\" ";
//     out << "WholeExtent=\"0 " << dm.ni - 1 << " 0 " << dm.nj - 1
//         << "\">\n";

//     /*output data stored on nodes (point data)*/
//     out << "<PointData>\n";

//     /*node volumes, scalar*/
//     out << "<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" "
//            "format=\"ascii\" type=\"Float64\">\n";
//     out << dm.node_area;
//     out << "</DataArray>\n";

//     /*potential, scalar*/
//     out << "<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\"
//     "
//            "type=\"Float64\">\n";
//     out << dm.phi;
//     out << "</DataArray>\n";

//     /*charge density, scalar*/
//     out << "<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\"
//     "
//            "type=\"Float64\">\n";
//     out << dm.rho;
//     out << "</DataArray>\n";

//     /*species number densities*/
//     for (Species &sp : species) {
//         out << "<DataArray Name=\"nd." << sp.name
//             << "\" NumberOfComponents=\"1\" format=\"ascii\" "
//                "type=\"Float64\">\n";
//         out << sp.den;
//         out << "</DataArray>\n";
//     }

//     /*electric field, 3 component vector*/
//     out << "<DataArray Name=\"ef\" NumberOfComponents=\"2\" format=\"ascii\"
//     "
//            "type=\"Float64\">\n";
//     out << dm.E;
//     out << "</DataArray>\n";

//     /*close out tags*/
//     out << "</PointData>\n";
//     out << "</ImageData>\n";
//     out << "</VTKFile>\n";
//     out.close();
// }
// void Output::fields(Domain &dm, std::vector<Species> &species) {
void Output::fields(Domain &domain) {
    std::stringstream name;
    name << "../results/fields_" << std::setfill('0') << std::setw(10)
         << domain.getTs() << ".dat";

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
    out << "ZONE i=" << domain.geo.ni << " j=" << domain.geo.nj
        << " SOLUTIONTIME=" << domain.time << "\n";
    for (int i = 0; i < domain.geo.ni; i++) {
        for (int j = 0; j < domain.geo.nj; j++) {
            out << i << " ";
            out << j << " ";
            out << domain.Hx(i, j) << " ";
            out << domain.Hy(i, j) << " ";
            //            out << sqrt(pow(domain.Hx(i, j), 2) + pow(domain.Hy(i,
            //            j), 2))
            //                << " ";
            out << domain.Hz(i, j) << " ";
            out << domain.Dx(i, j) << " ";
            out << domain.Dy(i, j) << " ";
            out << domain.Dz(i, j) << " ";
            out << domain.Jx(i, j) << " ";
            out << domain.Jy(i, j) << " ";
            out << domain.Jz(i, j) << "\n";
        }
    }

    // out << dm.E;

    out.close();
}

void Output::ProgressBar(double cur, double total) {
    double progress = static_cast<double>(cur) / static_cast<double>(total);
    int bar_width = 70;

    std::cout << "[";
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(2) << (progress * 100)
              << "%\r";
    std::cout.flush();
}

// writes information to the screen
void Output::screenOutput(Domain &domain, std::vector<Species> &species) {
    std::cout << "ts: " << domain.getTs();
    for (Species &sp : species)
        std::cout << std::setprecision(3) << "\t " << sp.name << ":"
                  << sp.getNp();
    std::cout << std::endl;
}

// file stream handle
namespace Output {
std::ofstream f_diag;
}

/*save runtime diagnostics to a file*/
void Output::diagOutput(Domain &domain, std::vector<Species> &species) {
    using namespace Output;  // to get access to f_diag

    // is the file open?
    if (!f_diag.is_open()) {
        f_diag.open("../runtime_diags.csv");
        f_diag << "ts,time,wall_time";
        for (Species &sp : species)
            f_diag << ",mp_count." << sp.name << ",real_count." << sp.name
                   << ",px." << sp.name << ",py." << sp.name << ",pz."
                   << sp.name << ",KE." << sp.name;
        f_diag << ",PE,E_total" << std::endl;
    }

    f_diag << domain.getTs() << "," << domain.getTime();
    f_diag << "," << domain.getWallTime();

    double tot_KE = 0;
    for (Species &sp : species) {
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

void Output::Display(Domain &domain, Field<Vec3d> field) {
    for (int i = 0; i < domain.geo.ni; i++) {
        for (int j = 0; j < domain.geo.nj; j++) {
            if (j == domain.geo.nj - 1)
                // std::cout << std::setprecision(2) << field[i][j] << "\n";
                std::cout << field[i][j] << "\n";
            else
                // std::cout << std::setprecision(2) << field[i][j] << "\t";
                std::cout << field[i][j] << "\t";
        }
    }
    std::cout << "----------------------------------------------\n";
}
