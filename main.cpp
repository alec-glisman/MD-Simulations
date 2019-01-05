#include "Simulation.h"

int main() {
    // Simulation variables
    long num_iter = 10000;
    int num_simulations = 1;
    unsigned long n_particle = 100;
    double temp = 4.0;
    double box = 20.0;
    double epsilon = 2.5;  // Depth of the potential well
    double sigma = 3.0;  // (Max) Length over which interactions can occur
    double dt = 0.0005;
    std::string filename = "01";
    std::string foldername = "Simulation-Output";
    int n_dump = 5; // Dump XYZ frame every 5 iterations

    // delete files (if they previously existed)
    std::string deleteCommand = "exec rm -r ~/CLionProjects/MD-Simulations/" + foldername;
    system(deleteCommand.c_str());

    // Initialize class
    Simulation sim1{num_iter, num_simulations,
                    n_particle, temp,
                    box, epsilon,
                    sigma, dt,
                    filename, foldername,
                    n_dump};
    sim1.main();  // Run Simulation

}
