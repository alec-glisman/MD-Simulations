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
Simulation::Simulation(long num_iter, int num_simulations,
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
    m_filename_xyz = "../" + m_foldername + "/" + m_filename + ".exyz";
    m_filename_log = "../" + m_foldername + "/" + m_filename + ".txt";
    m_filename_csv = "../" + m_foldername + "/" + m_filename + ".csv";


    m_DoF = m_n_dimensions * (n_particle - 1);

    // Initialize radii with FCC lattice positions
    fcc_lattice_init();
    // Initialize Velocities
    velocity_init();
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
    unsigned int row = 0;
    for (unsigned int i = 0; i < (m_n_particle - 1); i++) {
        // 4 Atoms in a unit cell
        for (auto a : rangeVect(0, 4)) {
            radius_.at(row).at(0) = (r_fcc.at(static_cast<unsigned int>(a)).at(0) + cartProd.at(i).at(0))
                    * cell_size / m_box;
            radius_.at(row).at(1) = (r_fcc.at(static_cast<unsigned int>(a)).at(1) + cartProd.at(i).at(1))
                    * cell_size / m_box;
            radius_.at(row).at(2) = (r_fcc.at(static_cast<unsigned int>(a)).at(2) + cartProd.at(i).at(2))
                    * cell_size / m_box;
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
    unsigned int row = 0;
    // Triple for-loop to go through every permutation of range
    for (int i = 0; i < vector.size(); i++) {
        for (int j = 0; j < vector.size(); j++) {
            for (int k = 0; k < vector.size(); k++) {
                // Fill in values
                cartProd.at(row).at(0) = i;
                cartProd.at(row).at(1) = j;
                cartProd.at(row).at(2) = k;
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
        string_t delimiter = "\t";
        // Header Information
        file << n_particleString + "\n"; // Number of particles
        file << "Lattice=\"" + boxString + " 0.0 0.0 0.0 "
                + boxString + " 0.0 0.0 0.0 "
                + boxString + "\" ";
        file << "Properties=pos:R:3:velo:R:3:species:S:1 "; // species:S:1:  :Radius:R:1
        file << "Time=" + timeString;
        file << "\n";
        // Atom Information (3-Position, 3-Velocity, 1-Radius)
        for (int i = 0; i < m_n_particle; i++) {
            file << std::fixed << std::setprecision(4) << radii[i][0];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << radii[i][1];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << radii[i][2];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << velocities[i][0];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << velocities[i][1];
            file << delimiter;
            file << std::fixed << std::setprecision(4) << velocities[i][2];
            file << delimiter;
            file << m_atom;
            /* file << delimiter;
            file << (m_diameter * 0.5); */
            file << "\n";
        }
        /* file << "\n\n"; */
        file.close();
    } else {
        std::cout << "Unable to open XYZ File";
    }
}


void Simulation::log_data() {
    std::ofstream file(m_filename_log.c_str(), std::fstream::app);
    if (file.is_open()) {
        file << "Simulation: " << m_filename << "\n";
        file << "\n";
        file << "Input Parameters:\n";
        file << "dt " << std::setprecision(4) << m_dt << "\n";
        file << "Total time " << std::setprecision(4) << m_t_total << "\n";
        file << "Number of particles " << m_n_particle << "\n";
        file << "Number of steps " << m_num_iter << "\n";
        file << "Box size " << std::setprecision(4) << m_temp << "\n";
        file << "Initial temperature " << std::setprecision(4) << m_temp << "\n";
        file << "Epsilon " << std::setprecision(4) << m_epsilon << "\n";
        file << "Sigma " << std::setprecision(4) << m_sigma << "\n";
        file.close();
    } else {
        std::cout << "Unable to open log File";
    }
}

void Simulation::forceAndEnergetics() {
    // Values to append to vectors at end of function
    double current_U{0};
    double current_W{0};
    doubleMatrix_t current_forces{m_n_particle, doubleVector_t(m_n_dimensions)};
    // Loop through all pairwise interactions where atom j > i
    for (unsigned long i = 0; i < m_n_particle; i++) {
        for (unsigned long j = (i + 1); j < m_n_particle; j++) {
            // Differential Distances
            double dx = radii.at(i).at(0) - radii.at(j).at(0);
            double dy = radii.at(i).at(1) - radii.at(j).at(1);
            double dz = radii.at(i).at(2) - radii.at(j).at(2);
            // Find smallest 'mirror' image using periodic boundary conditions
            dx = (dx - std::rint(dx)) * m_box;
            dy = (dy - std::rint(dy)) * m_box;
            dz = (dz - std::rint(dz)) * m_box;
            // Calculate |r_ij| = r2
            double r2 = (dx * dx) + (dy * dy) + (dz * dz);
            if (r2 < 0.00001) {
                r2 = 0.02;  // Prevent divide by zero
            }
            // Calculate (sigma/r_ij)^6
            double fr6 = pow(((m_sigma * m_sigma) / r2), 3);
            // Calculate overall LJ force
            double fpr = (48 * m_epsilon) * fr6 * (fr6 - 0.5) / r2;
            // Calculate component forces by multiplying with differential
            double fxi = fpr * dx;
            double fyi = fpr * dy;
            double fzi = fpr * dz;
            // Increment forces with current calculation, use anti-symmetry of force
            // f_ij = - f_ji
            current_forces.at(i).at(0) += fxi;
            current_forces.at(j).at(0) -= fxi;
            current_forces.at(i).at(1) += fyi;
            current_forces.at(j).at(1) -= fyi;
            current_forces.at(i).at(2) += fzi;
            current_forces.at(j).at(2) -= fzi;
            // Increment potential by interaction energy
            current_U += (4 * m_epsilon) * fr6 * (fr6 - 1);
            // Increment W by interaction energy as well
            current_W += -((dx * fxi) + (dy * fyi) + (dz * fzi));
        }
    }
    // Update class variables for current iteration
    E_potential.push_back(current_U);
    W.push_back(current_W);
    forces = current_forces;
}


void Simulation::kineticEnergy() {
    double current_KE{0};
    // Loop through all particles in system
    for (unsigned int i = 0; i < m_n_particle; i++) {
        current_KE += (velocities.at(i).at(0) * velocities.at(i).at(0))
                      + (velocities.at(i).at(1) * velocities.at(i).at(1))
                      + (velocities.at(i).at(2) * velocities.at(i).at(2));
    }
    // Divide by 2 to get 'proper' KE
    current_KE = current_KE / 2.0;
    // Append to KE vector
    E_kinetic.push_back(current_KE);
}


void Simulation::totalEnergy() {
    double current_E_total = E_kinetic.back() + E_potential.back();
    E_total.push_back(current_E_total);
}


void Simulation::simulationTemperature() {
    double current_T{};
    current_T = (2.0 / m_DoF) * E_kinetic.back(); // Use most recent KE calculation
    // Append to vector
    Temp_sim.push_back(current_T);
}


void Simulation::velocityVerlet() {
    // First half-timestep of algorithm
    for (unsigned int i = 0; i < m_n_particle; i++) {
        for (unsigned int j = 0; j < m_n_dimensions; j++) {
            velocities.at(i).at(j) += 0.5 * m_dt * forces.at(i).at(j);
            radii.at(i).at(j) += m_dt * velocities.at(i).at(j) / m_box;
            // Calculate 'unwrapped' radii
            unwrapped_radii.at(i).at(j) += m_dt * velocities.at(i).at(j) / m_box;
            // Applying Periodic Boundary Conditions to self.radii
            radii.at(i).at(j) -= std::rint(radii.at(i).at(j));
        }
    }
    // Second half-timestep
    forceAndEnergetics();  // Update forces
    for (unsigned int i = 0; i < m_n_particle; i++) {
        for (unsigned int j = 0; j < m_n_dimensions; j++) {
            velocities.at(i).at(j) += 0.5 * m_dt * forces.at(i).at(j);
        }
    }
    // Update system energetics
    kineticEnergy();
    simulationTemperature();
    totalEnergy();
    pressureCalc();
}


void Simulation::LJ_sim() {
    // Progress bar data
    ProgressBar progressBar(static_cast<unsigned int>(m_num_iter), 70);
    // NVE integration, Equilibration
    for ( ; m_iterationNumber < m_num_iter; m_iterationNumber++) {
        // Perform iteration step
        velocityVerlet();

        // Update Progress Bar
        ++progressBar;  // record the tick
        progressBar.display();  // display the bar

        // Dump frame (if at correct iteration number)
        if (m_iterationNumber % m_n_dump == 0) {
            create_frame();
        }
    }
    progressBar.done();  // tell the bar to finish
}


void Simulation::pressureCalc() {
    double P_id{m_rho * Temp_sim.back()};  // Ideal pressure using system Temp
    double P_ex{W.back() / (3.0 * m_vol)};  // Virial Pressure term
    double current_P{P_id + P_ex};
    Pressure.push_back(current_P); // Append current pressure
}


void Simulation::saveVars() {
    // Variables to save (csv)
    std::ofstream file(m_filename_csv.c_str()); // , std::fstream::app);
    if (file.is_open()) {
        // Header information
        string_t delimiter{","};
        file << "Iteration" << delimiter;
        file << "E_total" << delimiter;
        file << "E_kinetic" << delimiter;
        file << "E_potential" << delimiter;
        file << "Temperature" << delimiter;
        file << "Pressure" << delimiter;
        file << "\n";
        // Loop through vectors to output information
        for (unsigned int i = 0; i < E_total.size(); i++) {
            file << i << delimiter;
            file << std::setprecision(6) << std::fixed << E_total.at(i) << delimiter;
            file << std::setprecision(6) << std::fixed << E_kinetic.at(i) << delimiter;
            file << std::setprecision(6) << std::fixed << E_potential.at(i) << delimiter;
            file << std::setprecision(6) << std::fixed << Temp_sim.at(i) << delimiter;
            file << std::setprecision(6) << std::fixed << Pressure.at(i) << delimiter;
            file << "\n";
        }
        file.close();
    } else {
        std::cout << "Unable to open save file";
    }
}


void Simulation::plots() {
    // Energy vs time
    std::string energyPlot = "cd ../gnu-plots; gnuplot  -p energy_time.gnuplot";
    system(energyPlot.c_str());
}


void Simulation::main() {
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "\n"; // Just for aesthetics on output

    // Main simulation with Lennard-Jones potential
    LJ_sim();
    // Save variables of interest
    saveVars();
    // Make Plots
    plots();

    // End timer
    auto stop = std::chrono::high_resolution_clock::now();
    auto durationMin = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
    auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by main: " << durationMin.count() << " minutes, ";
    std::cout << durationSec.count() % 60 << " seconds \n";
    std::cout << "\t" << m_num_iter / durationSec.count() << " iterations/second \n";
    std::cout << "\n"; // Just for aesthetics on output

    // Output timing data to log file
    std::ofstream file(m_filename_log.c_str(), std::fstream::app);
    if (file.is_open()) {
        file << "Time taken by main: " << durationMin.count() << " minutes, ";
        file << durationSec.count() % 60 << " seconds \n";
        file << "\t" << m_num_iter / durationSec.count() << " iterations/second \n";
        file << "\n\n";
        file.close();
    } else {
        std::cout << "Unable to open log File";
    }
}
