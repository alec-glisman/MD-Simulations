#include "Simulation.h"

int main() {
    // Make class with uniform initialization
    Simulation sim_test_default{};

    // Simulation variables
    int num_iter = 2000;
    int num_simulations = 1;
    unsigned long n_particle = 256;
    double temp = 1.0;
    double box = 5.0;
    double epsilon = 1.0;
    double sigma = 1.0;
    double dt = 1e-4;
    std::string filename = "test_file";
    std::string foldername = "test_folder";

    // Initialize class
    Simulation sim_test_simVarConst{num_iter, num_simulations,
                                    n_particle, temp,
                                    box, epsilon,
                                    sigma, dt,
                                    filename, foldername};


    return 0;
}