#include "Simulation.h"

int main() {
    // Simulation variables
    long num_iter = 6000;
    int num_simulations = 1;
    unsigned long n_particle = 100;
    double temp = 1.0;
    double box = 15.0;
    double epsilon = 3.0;  // Depth of the potential well
    double sigma = 1.0;  // (Max) Length over which interactions can occur
    double dt = 0.005;
    std::string filename = "01";
    std::string foldername = "Simulation-Output";
    int n_dump = 10; // Dump XYZ frame every n iterations

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
