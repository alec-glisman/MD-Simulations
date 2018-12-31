//
// Created by Alec Glisman on 2018-12-28.
//

// Internal Project Dependencies
#include "Simulation.h"

// External Dependencies


// TYPE ALIASES
using doubleMatrix_t =  std::vector<std::vector<double> >;
using doubleVector_t =  std::vector<double>;
using intMatrix_t =  std::vector<std::vector<int> >;
using intVector_t =  std::vector<int>;
using string_t = std::string;


Simulation::Simulation() = default;


// Constructor to input simulation parameters
Simulation::Simulation(int num_iter, int num_simulations,
                       unsigned long n_particle, double temp,
                       double box, double epsilon,
                       double sigma, double dt,
                       string_t filename, string_t foldername) :
        m_num_iter(num_iter), m_num_simulations(num_simulations),
        m_temp(temp), m_n_particle(n_particle),
        m_box(box), m_epsilon(epsilon),
        m_sigma(sigma), m_dt(dt),
        m_t_total(dt * num_iter),
        m_vol(box * box * box), m_rho(n_particle / (box * box * box)),
        m_filename(std::move(filename)),
        m_foldername(std::move(foldername)) {
    // Initialize radii with FCC lattice positions
    fcc_lattice_init();
    Utilities::print(radii, "FCC Radii:");
}


void Simulation::fcc_lattice_init() {
    // Calculate the number of and size of unit cells needed
    double base = static_cast<double>(m_n_particle) / 4.0;
    double exponent = 1.0 / 3.0;
    auto cells = static_cast<int>(ceil(pow(base, exponent)));
    double cell_size = m_box / cells;

    // Init initial radius matrix
    doubleMatrix_t radius_{m_n_particle, doubleVector_t(m_n_dimensions)};
    // FCC atom centers in unit cell
    doubleMatrix_t r_fcc{{0.25, 0.25, 0.25},
                         {0.25, 0.75, 0.75},
                         {0.75, 0.75, 0.25},
                         {0.75, 0.25, 0.75}};

    // Get range for Cartesian product
    intVector_t range = rangeVect(0, cells);
    // Cartesian Product (Integers)
    intMatrix_t cartProd = cartesianProduct3(range);

    // Iterate through each row of Cartesian Product until n_particles reached
    int row = 0;
    for (int i = 0; i < (m_n_particle - 1); i++) {
        // 4 Atoms in a unit cell
        for (int a : rangeVect(0, 4)) {
            radius_[row][0] = (r_fcc[a][0] + cartProd[i][0]) * cell_size / m_box;
            radius_[row][1] = (r_fcc[a][1] + cartProd[i][1]) * cell_size / m_box;
            radius_[row][2] = (r_fcc[a][2] + cartProd[i][2]) * cell_size / m_box;
            row++;
            if (row == m_n_particle) {
                // Set Initial Radii positions
                radii = radius_;
                return;
            }
        }
    }
}


intMatrix_t Simulation::cartesianProduct3(intVector_t &vector) {
    auto numRows = static_cast<unsigned long>(pow(vector.size(), 3));
    unsigned long numCols = 3;
    intMatrix_t cartProd{numRows,
                         intVector_t(numCols)};
    int row = 0;
    // Triple for-loop to go through every permutation of range
    for (int i = 0; i < vector.size(); i++) {
        for (int j = 0; j < vector.size(); j++) {
            for (int k = 0; k < vector.size(); k++) {
                // Fill in values
                cartProd[row][0] = i;
                cartProd[row][1] = j;
                cartProd[row][2] = k;
                // Increment to next row in matrix to fill
                row++;
            }
        }
    }
    return cartProd;
}


intVector_t Simulation::rangeVect(int start, int stop) {
    // Get range variable to iterate through
    intVector_t range;
    boost::push_back(range, boost::irange(start, stop));
    //toPrint(range, "Cell Range:");
    return range;
}
