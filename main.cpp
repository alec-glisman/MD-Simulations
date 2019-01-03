#include "Simulation.h"


int main() {
    // Simulation variables
    long num_iter = 1000;
    int num_simulations = 1;
    unsigned long n_particle = 80;
    double temp = 7.0;
    double box = 5.0;
    double epsilon = 2.0;
    double sigma = 1.0;
    double dt = 0.005;
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
    sim1.main();  // Run Simulation

}
