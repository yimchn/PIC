//#include "output.h"
//
//#include <fstream>
//#include <iomanip>
//#include <iostream>
//#include <sstream>
//
//#include "domain.h"
//#include "species.h"
//
//template <typename T>
//void Output::fields(Domain& dm) {
//	std::stringstream name;
//	name << "../results/fields_" << std::setfill('0') << std::setw(10)
//		<< dm.getTs() << ".dat";
//
//	/*open output file*/
//	std::ofstream out(name.str());
//	if (!out.is_open()) {
//		std::cerr << "Could not open " << name.str() << std::endl;
//		return;
//	}
//
//	out << "TITLE = Magnetic Field 2d\n";
//	out << "VARIABLES = \"X\", "
//		"\"Y\",\"Hx\",\"Hy\",\"Hz\",\"Ex\",\"Ey\",\"Ez\",\"Jx\",\"Jy\","
//		"\"Jz\"\n";
//	out << "ZONE i=" << dm.geo.ni << " j=" << dm.geo.nj
//		<< " SOLUTIONTIME=" << dm.time << "\n";
//	for (int i = 0; i < dm.geo.ni; i++) {
//		for (int j = 0; j < dm.geo.nj; j++) {
//			out << i << " ";
//			out << j << " ";
//			out << dm.Hx(i, j) << " ";
//			out << dm.Hy(i, j) << " ";
//			out << dm.Hz(i, j) << " ";
//			out << dm.Dx(i, j) << " ";
//			out << dm.Dy(i, j) << " ";
//			out << dm.Dz(i, j) << " ";
//			out << dm.Jx(i, j) << " ";
//			out << dm.Jy(i, j) << " ";
//			out << dm.Jz(i, j) << "\n";
//		}
//	}
//
//	out.close();
//}
//
//void Output::fields(Domain& domain, int step) {
//	std::stringstream name;
//	name << "../results/fields_" << std::setfill('0') << std::setw(10) << step
//		<< ".dat";
//
//	/*open output file*/
//	std::ofstream out(name.str());
//	if (!out.is_open()) {
//		std::cerr << "Could not open " << name.str() << std::endl;
//		return;
//	}
//
//	out << domain.Dz;
//
//	out.close();
//}
//
//void Output::ProgressBar(Domain dm) {
//	using namespace indicators;
//	BlockProgressBar{
//		option::BarWidth{80},
//		option::ForegroundColor{Color::white},
//		option::FontStyles{
//			std::vector<FontStyle>{FontStyle::bold}},
//			option::MaxProgress{dm.num_ts}
//	};
//
//	std::cout << "Total number of iterations: " << dm.num_ts << "\n";
//	for (size_t i = 0; i < dm.num_ts; ++i) {
//		// 执行相应的操作
//		result.emplace_back(i * i);
//
//		// 在进度条后面显示迭代进度
//		bar.set_option(option::PostfixText{
//			std::to_string(i) + "/" + std::to_string(dm.num_ts);
//			});
//
//		// 更新进度条
//		bar.tick();
//	}
//
//	bar.mark_as_completed();
//
//	show_console_cursor(true);
//
//	//double progress =
//	//	static_cast<double>(dm.ts) / static_cast<double>(dm.num_ts);
//	//int bar_width = 70;
//
//	//std::cout << "[";
//	//int pos = bar_width * progress;
//	//for (int i = 0; i < bar_width; ++i) {
//	//	if (i < pos)
//	//		std::cout << "=";
//	//	else if (i == pos)
//	//		std::cout << ">";
//	//	else
//	//		std::cout << " ";
//	//}
//	//std::cout << "] " << std::fixed << std::setprecision(2) << (progress * 100)
//	//	<< "%\r";
//	//std::cout.flush();
//}
//
//// writes information to the screen
//void Output::screenOutput(Domain& domain, std::vector<Species>& species) {
//	std::cout << "ts: " << domain.getTs();
//	for (Species& sp : species)
//		std::cout << std::setprecision(3) << "\t " << sp.name << ":"
//		<< sp.getNp();
//	std::cout << std::endl;
//}
//
//// file stream handle
//namespace Output {
//	std::ofstream f_diag;
//}
//
///*save runtime diagnostics to a file*/
//void Output::diagOutput(Domain& domain, std::vector<Species>& species) {
//	using namespace Output;  // to get access to f_diag
//
//	// is the file open?
//	if (!f_diag.is_open()) {
//		f_diag.open("../runtime_diags.csv");
//		f_diag << "ts,time,wall_time";
//		for (Species& sp : species)
//			f_diag << ",mp_count." << sp.name << ",real_count." << sp.name
//			<< ",px." << sp.name << ",py." << sp.name << ",pz."
//			<< sp.name << ",KE." << sp.name;
//		f_diag << ",PE,E_total" << std::endl;
//	}
//
//	f_diag << domain.getTs() << "," << domain.getTime();
//	f_diag << "," << domain.getWallTime();
//
//	double tot_KE = 0;
//	for (Species& sp : species) {
//		double KE = sp.getKE();  // species kinetic energy
//		tot_KE += KE;            // increment total energy
//		Vec2d mom = sp.getMomentum();
//
//		f_diag << "," << sp.getNp() << "," << sp.getRealCount() << "," << mom[0]
//			<< "," << mom[1] << "," << mom[2] << "," << KE;
//	}
//
//	// write out system potential and total energy
//	double PE = domain.getPE();
//	f_diag << "," << PE << "," << (tot_KE + PE);
//
//	f_diag << "\n";  // use \n to avoid flush to disc
//	if (domain.getTs() % 25 == 0) f_diag.flush();
//}
//
//void Output::Display(Domain& domain, Field<Vec3d> field) {
//	for (int i = 0; i < domain.geo.ni; i++) {
//		for (int j = 0; j < domain.geo.nj; j++) {
//			if (j == domain.geo.nj - 1)
//				// std::cout << std::setprecision(2) << field[i][j] << "\n";
//				std::cout << field[i][j] << "\n";
//			else
//				// std::cout << std::setprecision(2) << field[i][j] << "\t";
//				std::cout << field[i][j] << "\t";
//		}
//	}
//	std::cout << "----------------------------------------------\n";
//}
//
//void Output::B(Domain& dm) {
//	std::stringstream name;
//	name << "../results/B_" << std::setfill('0') << std::setw(10)
//		<< dm.getTs() << ".dat";
//
//	/*open output file*/
//	std::ofstream out(name.str());
//	if (!out.is_open()) {
//		std::cerr << "Could not open " << name.str() << std::endl;
//		return;
//	}
//
//	// out << "TITLE = Magnetic Field for FDTD 2d case\n";
//	// out << "VARIABLES = \"X\", "
//	//        "\"Y\",\"Bx\",\"By\",\"Bz\",\"B\"\n";
//	// out << " SOLUTIONTIME=" << dm.time << "\n";
//	// for (int i = 0; i < dm.geo.ni; i++) {
//	//     for (int j = 0; j < dm.geo.nj; j++) {
//	//         out << i << " ";
//	//         out << j << " ";
//	//         out << dm.Hx(i, j) << " ";
//	//         out << dm.Hy(i, j) << " ";
//	//         out << dm.Hz(i, j) << " ";
//	//         out << sqrt(pow(dm.Hx(i, j), 2) + pow(dm.Hy(i, j), 2)) << "\n";
//	//     }
//	// }
//	out << dm.Hy;
//
//	out.close();
//}
//
//void Output::E(Domain& dm) {
//	std::stringstream name;
//	name << "../results/E_" << std::setfill('0') << std::setw(10)
//		<< dm.getTs() << ".dat";
//
//	/*open output file*/
//	std::ofstream out(name.str());
//	if (!out.is_open()) {
//		std::cerr << "Could not open " << name.str() << std::endl;
//		return;
//	}
//
//	// out << "TITLE = Electric Field for FDTD 2d case\n";
//	// out << "VARIABLES = \"X\", "
//	//        "\"Y\",\"Ex\",\"Ey\",\"Ez\",\"E\"\n";
//	// out << "ZONE i=" << dm.geo.ni << " j=" << dm.geo.nj
//	//     << " SOLUTIONTIME=" << dm.time << "\n";
//	// for (int i = 0; i < dm.geo.ni; i++) {
//	//     for (int j = 0; j < dm.geo.nj; j++) {
//	//         out << i << " ";
//	//         out << j << " ";
//	//         out << dm.Dx(i, j) << " ";
//	//         out << dm.Dy(i, j) << " ";
//	//         out << dm.Dz(i, j) << " ";
//	//         out << sqrt(pow(dm.Dx(i, j), 2) + pow(dm.Dy(i, j), 2)) << "\n";
//	//     }
//	// }
//	out << dm.Dz;
//
//	out.close();
//}
//
//void Output::J(Domain& dm) {
//	std::stringstream name;
//	name << "../results/J/J_" << std::setfill('0') << std::setw(10)
//		<< dm.getTs() << ".dat";
//
//	/*open output file*/
//	std::ofstream out(name.str());
//	if (!out.is_open()) {
//		std::cerr << "Could not open " << name.str() << std::endl;
//		return;
//	}
//
//	out << "TITLE = Current density for FDTD 2d case\n";
//	out << "VARIABLES = \"X\", "
//		"\"Y\",\"Jx\",\"Jy\",\"Jz\",\"J\"\n";
//	out << "ZONE i=" << dm.geo.ni << " j=" << dm.geo.nj
//		<< " SOLUTIONTIME=" << dm.time << "\n";
//	for (int i = 0; i < dm.geo.ni; i++) {
//		for (int j = 0; j < dm.geo.nj; j++) {
//			out << i << " ";
//			out << j << " ";
//			out << dm.Jx(i, j) << " ";
//			out << dm.Jy(i, j) << " ";
//			out << dm.Jz(i, j) << " ";
//			out << sqrt(pow(dm.Jx(i, j), 2) + pow(dm.Jy(i, j), 2)) << "\n";
//		}
//	}
//
//	out.close();
//}