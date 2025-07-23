#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cmath>
#include <limits>
#include "ConfocalParabolicBilliard.hpp"
#include "BoundaryWallMethod.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]){

    if(argc<4){
        std::cerr<<"Uso: "<<argv[0]<<" xi0 eta0 k2\n";
        return 1;
    }
    double xi0   = std::stod(argv[1]);
    double eta0  = std::stod(argv[2]);
    double k2    = std::stod(argv[3]);
    double angle_deg = std::stod(argv[4]);
    int    threads = std::atoi(argv[5]);
    double angle = angle_deg * (M_PI / 180.0);
    double gamma = std::numeric_limits<double>::infinity();
    int    num_b = 200;
    int    Ngrid = 200;

    // Crear el directorio para los resultados (ignora si ya existe)
    const std::string results_dir = "resultados";
    fs::create_directory(results_dir);

    // Construir la ruta de los archivos
    fs::path pathB = fs::path(results_dir) / ("boundary.dat");
    fs::path pathD = fs::path(results_dir) / ("density_" + std::to_string(k2) + "_" + std::to_string(angle_deg) + ".dat");
    fs::path pathP = fs::path(results_dir) / ("phase_" + std::to_string(k2) + "_" + std::to_string(angle_deg) + ".dat");


    // 1) frontera
    bwm::ConfocalParabolicBilliard billiard(xi0, eta0, num_b);
    auto boundary = billiard.getBoundary();
    std::ofstream foutB(pathB);
    for(auto &p: boundary) foutB<<p.x()<<" "<<p.y()<<"\n";
    foutB.close();

    // 2) solver BWM
    std::vector<double> gvec{gamma};
    bwm::BoundaryWallMethod solver(boundary, gvec, std::sqrt(k2), angle, threads);

    // 3) malla y volcado de density.dat y phase.dat
    std::ofstream foutD(pathD);
    std::ofstream foutP(pathP);
    double xmin=-8, xmax=8;
    for(int i=0;i<Ngrid;++i){
        double y = xmin + i*(xmax-xmin)/(Ngrid-1);
        for(int j=0;j<Ngrid;++j){
            double x   = xmin + j*(xmax-xmin)/(Ngrid-1);
            auto ψ     = solver.computeScatteredWave({{x,y}})[0];
            foutD<<std::norm(ψ)<<(j+1==Ngrid? "\n":" ");
            foutP<<std::arg(ψ)<<(j+1==Ngrid? "\n":" ");
        }
    }
    return 0;
}

