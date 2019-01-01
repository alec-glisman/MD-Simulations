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

    // Full filenames/paths
    m_filename_xyz = "../" + m_foldername + "/" + m_filename + ".xyz";
    m_filename_log = "../" + m_foldername + "/" + m_filename + ".txt";


    m_DoF = m_n_dimensions * (m_n_particle - 1);

    // Initialize radii with FCC lattice positions
    fcc_lattice_init();

    // Initialize Velocities
    velocity_init();
    //Utilities::print(velocities, "Velocity at t=0:");

    // Make initial FCC frame
    create_frame();
    // Output Log Data
    log_data();
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


void Simulation::velocity_init() {
    doubleVector_t vel_init_i(m_n_particle);
    doubleVector_t vel_init_j(m_n_particle);
    doubleVector_t vel_init_k(m_n_particle);
    // Initialize Velocities from Gaussian (Boltzman) Distribution
    std::seed_seq seed{43}; // Set seed
    std::minstd_rand0 generator(seed); // Make generator for random values
    double mean = 0;
    double std_dev = pow(m_temp, 0.5);
    std::normal_distribution<double> distribution(mean, std_dev); // Set distribution to sample from
    // Loop over all of velocity matrix to give values
    for (unsigned long p = 0; p < m_n_particle; p++) {
        vel_init_i.at(p) = distribution(generator);
        vel_init_j.at(p) = distribution(generator);
        vel_init_k.at(p) = distribution(generator);
    }

    // Calculate the mean of each column vector (i, j, k)
    double mean_i = std::accumulate(vel_init_i.begin(), vel_init_i.end(), 0.0) / m_n_particle;
    double mean_j = std::accumulate(vel_init_j.begin(), vel_init_j.end(), 0.0) / m_n_particle;
    double mean_k = std::accumulate(vel_init_k.begin(), vel_init_k.end(), 0.0) / m_n_particle;
    // Subtract the mean from all entries to ensure average linear momentum = 0
    for (unsigned long p = 0; p < m_n_particle; p++) {
        vel_init_i.at(p) -= mean_i;
        vel_init_j.at(p) -= mean_j;
        vel_init_k.at(p) -= mean_k;
    }
    // Ensure total sum = 0 once more by making row 0 the opposite sum of row 1 -> n
    double tot_i = std::accumulate(vel_init_i.begin(), vel_init_i.end(), -vel_init_i.at(0));
    double tot_j = std::accumulate(vel_init_j.begin(), vel_init_j.end(), -vel_init_j.at(0));
    double tot_k = std::accumulate(vel_init_k.begin(), vel_init_k.end(), -vel_init_k.at(0));
    vel_init_i.at(0) = -tot_i;
    vel_init_j.at(0) = -tot_j;
    vel_init_k.at(0) = -tot_k;

    // Make sure initialized velocities equal system Temp
    double v2_total{0.0};
    for (unsigned long p = 0; p < m_n_particle; p++) {
        v2_total += (vel_init_i.at(p) * vel_init_i.at(p))
                    + (vel_init_j.at(p) * vel_init_j.at(p))
                    + (vel_init_k.at(p) * vel_init_k.at(p));
    }
    v2_total = 0.5 * v2_total;  // Scale by 1/2 to get Kinetic Energy
    double kineticThermometer = m_temp * m_DoF * 0.5;
    double vel_KE_factor2 = kineticThermometer / v2_total;
    double vel_KE_factor = pow(vel_KE_factor2, 0.5);
    // Scale velocities by the thermometer KE
    for (unsigned long p = 0; p < m_n_particle; p++) {
        // Scale the Values
        vel_init_i.at(p) *= vel_KE_factor;
        vel_init_j.at(p) *= vel_KE_factor;
        vel_init_k.at(p) *= vel_KE_factor;
        // Enter values in the velocity array
        velocities.at(p).at(0) = vel_init_i.at(p);
        velocities.at(p).at(1) = vel_init_j.at(p);
        velocities.at(p).at(2) = vel_init_k.at(p);
    }
}


void Simulation::create_frame() {
    std::cout << std::setprecision(4) << std::fixed;
    // Convert n_particles to a string object
    std::stringstream ssPart; // now we can make a new ss
    ssPart << m_n_particle;
    string_t n_particleString = ssPart.str();
    // Convert box (double) to a string object
    std::stringstream ssBox; // now we can make a new ss
    ssBox << m_box;
    string_t boxString = ssBox.str();
    // Convert time (double) to a string object
    std::stringstream ssTime; // now we can make a new ss
    ssTime << (m_dt * m_iterationNumber);
    string_t timeString = ssTime.str();

    // File to write to
    std::ofstream file(m_filename_xyz.c_str(), std::fstream::app);
    if (file.is_open()) {
        string_t delimiter = "      ";
        // Header Information
        file << n_particleString + "\n"; // Number of particles
        file << "Lattice=\"" + boxString + " 0.0 0.0 0.0 "
                + boxString + " 0.0 0.0 0.0 "
                + boxString + "\" ";
        file << "Properties=\"species:S:1:pos:R:3:vel:R:3:diameter:R:1\" ";
        file << "Time=" + timeString + "\n";
        // Atom Information (1-Species, 3-Position, 3-Velocity, 1-Diameter)
        for (int i = 0; i < m_n_particle; i++) {
            file << m_atom + "\t\t";
            file << radii[i][0];
            file << delimiter;
            file << radii[i][1];
            file << delimiter;
            file << radii[i][2];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << velocities[i][0];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << velocities[i][1];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << velocities[i][2];
            file << delimiter;
            file << m_diameter;
            file << "\n";
        }
        file << "\n\n";
        file.close();
    } else {
        std::cout << "Unable to open XYZ File";
    }
}


void Simulation::log_data() {
    std::ofstream file(m_filename_log.c_str(), std::fstream::app);
    if (file.is_open()) {
        file << "\n\n";
        file << "Simulation: " << m_filename << "\n";
        file << "\n";
        file << "Input Parameters:\n";
        file << "dt " << std::setprecision(4) << m_dt << "\n";
        file << "Total time " << std::setprecision(4) << m_t_total << "\n";
        file << "Number of particles " << m_n_particle << "\n";
        file << "Number of steps " << m_num_iter << "\n";
        file << "Box size " << std::setprecision(4) << m_temp << "\n";
        file << "Initial temperature " << std::setprecision(4) << m_temp << "\n";
        file << "Epsilon " << std::setprecision(4) << m_temp << "\n";
        file << "Sigma " << std::setprecision(4) << m_sigma << "\n";
        file.close();
    } else {
        std::cout << "Unable to open log File";
    }
}
