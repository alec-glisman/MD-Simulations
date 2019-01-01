#include "Simulation.h"


int main() {
    // Simulation variables
    int num_iter = 2000;
    int num_simulations = 1;
    unsigned long n_particle = 256;
    double temp = 1.0;
    double box = 100.0;
    double epsilon = 1.0;
    double sigma = 1.0;
    double dt = 1e-4;
    std::string filename = "test_file";
    std::string foldername = "test_folder";

    // delete files (if they previously existed)
    std::string deleteCommand = "exec rm -r ~/CLionProjects/MD-Simulations/" + foldername + "/*";
    system(deleteCommand.c_str());



    // Initialize class
    Simulation sim1{num_iter, num_simulations,
                    n_particle, temp,
                    box, epsilon,
                    sigma, dt,
                    filename, foldername};
    temp = 20;
    Simulation sim2{num_iter, num_simulations,
                    n_particle, temp,
                    box, epsilon,
                    sigma, dt,
                    filename, foldername};
    return 0;
}