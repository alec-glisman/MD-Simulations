//
// Created by Alec Glisman on 2018-12-28.
//

// Header file with definitions
#include "Simulation.h"


Simulation::Simulation() = default;


// Constructor to input simulation parameters
Simulation::Simulation(long num_iter, int num_simulations,
                       unsigned long n_particle, double temp,
                       double box, double epsilon,
                       double sigma, double dt,
                       string_t filename, string_t foldername,
                       int n_dump) :
        m_num_iter(num_iter), m_num_simulations(num_simulations),
        m_temp(temp), m_n_particle(n_particle),
        m_box(box), m_epsilon(epsilon),
        m_sigma(sigma), m_dt(dt),
        m_t_total(dt * num_iter),
        m_vol(box * box * box), m_rho(n_particle / (box * box * box)),
        m_filename(std::move(filename)),
        m_foldername(std::move(foldername)),
        m_n_dump(n_dump){

    // Full filenames/paths
    m_filename_xyz = "../" + m_foldername + "/" + m_filename + ".exyz";
    m_filename_log = "../" + m_foldername + "/" + m_filename + ".txt";
    m_filename_csv = "../" + m_foldername + "/" + m_filename + ".csv";
    // Make directory if it does not already exist
    std::string mkdir = "mkdir ../" + m_foldername;
    system(mkdir.c_str());

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
        // Divide by 2 to get atoms to all start in a subset of the box in X-axis
        for (auto a : rangeVect(0, 4)) {
            radius_[row][0] = (r_fcc[static_cast<unsigned int>(a)][0] + cartProd[i][0])
                              * cell_size / m_box;
            radius_[row][1] = (r_fcc[static_cast<unsigned int>(a)][1] + cartProd[i][1])
                              * cell_size / m_box;
            radius_[row][2] = (r_fcc[static_cast<unsigned int>(a)][2] + cartProd[i][2])
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
    return range;
}


void Simulation::velocity_init() {
    doubleVector_t vel_init_i(m_n_particle);
    doubleVector_t vel_init_j(m_n_particle);
    doubleVector_t vel_init_k(m_n_particle);
    // Initialize Velocities from Gaussian (Boltzmann) Distribution
    std::seed_seq seed{43}; // Set seed
    std::minstd_rand0 generator(seed); // Make generator for random values
    double mean = 0;
    double std_dev = pow(m_temp, 0.5);
    std::normal_distribution<double> distribution(mean, std_dev); // Set distribution to sample from
    // Loop over all of velocity matrix to give values
    for (unsigned long p = 0; p < m_n_particle; p++) {
        vel_init_i[p] = distribution(generator);
        vel_init_j[p] = distribution(generator);
        vel_init_k[p] = distribution(generator);
    }

    // Calculate the mean of each column vector (i, j, k)
    double mean_i = std::accumulate(vel_init_i.begin(), vel_init_i.end(), 0.0) / m_n_particle;
    double mean_j = std::accumulate(vel_init_j.begin(), vel_init_j.end(), 0.0) / m_n_particle;
    double mean_k = std::accumulate(vel_init_k.begin(), vel_init_k.end(), 0.0) / m_n_particle;
    // Subtract the mean from all entries to ensure average linear momentum = 0
    for (unsigned long p = 0; p < m_n_particle; p++) {
        vel_init_i[p] -= mean_i;
        vel_init_j[p] -= mean_j;
        vel_init_k[p] -= mean_k;
    }
    // Ensure total sum = 0 once more by making row 0 the opposite sum of row 1 -> n
    double tot_i = std::accumulate(vel_init_i.begin(), vel_init_i.end(), -vel_init_i[0]);
    double tot_j = std::accumulate(vel_init_j.begin(), vel_init_j.end(), -vel_init_j[0]);
    double tot_k = std::accumulate(vel_init_k.begin(), vel_init_k.end(), -vel_init_k[0]);
    vel_init_i[0] = -tot_i;
    vel_init_j[0] = -tot_j;
    vel_init_k[0] = -tot_k;

    // Make sure initialized velocities equal system Temp
    double v2_total{0.0};
    for (unsigned long p = 0; p < m_n_particle; p++) {
        v2_total += (vel_init_i[p] * vel_init_i[p])
                    + (vel_init_j[p] * vel_init_j[p])
                    + (vel_init_k[p] * vel_init_k[p]);
    }
    v2_total = 0.5 * v2_total;  // Scale by 1/2 to get Kinetic Energy
    double kineticThermometer = m_temp * m_DoF * 0.5;
    double vel_KE_factor2 = kineticThermometer / v2_total;
    double vel_KE_factor = pow(vel_KE_factor2, 0.5);
    // Scale velocities by the thermometer KE
    for (unsigned long p = 0; p < m_n_particle; p++) {
        // Scale the Values
        vel_init_i[p] *= vel_KE_factor;
        vel_init_j[p] *= vel_KE_factor;
        vel_init_k[p] *= vel_KE_factor;
        // Enter values in the velocity array
        velocities[p][0] = vel_init_i[p];
        velocities[p][1] = vel_init_j[p];
        velocities[p][2] = vel_init_k[p];
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
        for (unsigned int i = 0; i < m_n_particle; i++) {
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
            file << "\n";
        }
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
        file << "Volume " << std::setprecision(4) << m_vol << "\n";
        file << "Density " << std::setprecision(4) << m_rho << "\n";
        file << "Initial temperature " << std::setprecision(4) << m_temp << "\n";
        file << "Epsilon " << std::setprecision(4) << m_epsilon << "\n";
        file << "Sigma " << std::setprecision(4) << m_sigma << "\n";
        file.close();
    } else {
        std::cout << "Unable to open log File";
    }
}

void Simulation::forceAndEnergetics() {
    // Variables for for-loop
    unsigned long n = m_n_particle; // Number of bodies
    unsigned long r = n * (n - 1) / 2; // Number of interactions to count
    // New time-step variables
    double current_U{0};
    double current_W{0};
    doubleMatrix_t current_forces{n, doubleVector_t(n)};

#pragma omp parallel
    {
        // Values to append to vectors at end of function
        double local_U{0};
        double local_W{0};
        doubleMatrix_t local_forces{n, doubleVector_t(n)};
        doubleMatrix_t local_radii{radii};

        // Collapsed lower-triangle matrix double for-loop: Loop through all pairwise interactions where atom j > i
        #pragma omp for schedule(static)
        for (unsigned k = 0; k < r; k++) {
            // i and j convention.
            // Source: https://stackoverflow.com/questions/33810187/openmp-and-c-private-variables/33836073#33836073
            unsigned long i = k / n;
            unsigned long j = k % n;
            if (j <= i) {
                i = n - i - 2;
                j = n - j - 1;
            }
            // Differential Distances
            double dx = local_radii[i][0] - local_radii[j][0];
            double dy = local_radii[i][1] - local_radii[j][1];
            double dz = local_radii[i][2] - local_radii[j][2];
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
            local_forces[i][0] += fxi;
            local_forces[j][0] -= fxi;
            local_forces[i][1] += fyi;
            local_forces[j][1] -= fyi;
            local_forces[i][2] += fzi;
            local_forces[j][2] -= fzi;
            // Increment potential by interaction energy
            local_U += (4.0 * m_epsilon) * fr6 * (fr6 - 1);
            // Increment W by interaction energy as well
            local_W += -((dx * fxi) + (dy * fyi) + (dz * fzi));
        } // End of lower-triangular for-loop

        // Sum local variables after iterations complete
        #pragma omp critical
        {
            for (unsigned long i = 0; i < m_n_particle; i++) {
                current_forces[i][0] += local_forces[i][0];
                current_forces[i][1] += local_forces[i][1];
                current_forces[i][2] += local_forces[i][2];
            }
            current_U += local_U;
            current_W += local_W;
        }
    } // End of pragma omp parallel

    // Update class variables for current iteration
    E_potential.push_back(current_U);
    W.push_back(current_W);
    forces = current_forces;
}


void Simulation::kineticEnergy() {
    double current_KE{0};
    // Loop through all particles in system
    for (unsigned int i = 0; i < m_n_particle; i++) {
        current_KE += (velocities[i][0] * velocities[i][0])
                      + (velocities[i][1] * velocities[i][1])
                      + (velocities[i][2] * velocities[i][2]);
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
            velocities[i][j] += 0.5 * m_dt * forces[i][j];
            radii[i][j] += m_dt * velocities[i][j] / m_box;
            // Calculate 'unwrapped' radii
            unwrapped_radii[i][j] += m_dt * velocities[i][j] / m_box;
            // Applying Periodic Boundary Conditions to self.radii
            radii[i][j] -= std::rint(radii[i][j]);
        }
    }
    // Second half-timestep
    forceAndEnergetics();  // Update forces
    for (unsigned int i = 0; i < m_n_particle; i++) {
        for (unsigned int j = 0; j < m_n_dimensions; j++) {
            velocities[i][j] += 0.5 * m_dt * forces[i][j];
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
    ProgressBar progressBar(static_cast<unsigned int>(m_num_iter), 60);
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
            file << std::setprecision(6) << std::fixed << E_total[i] << delimiter;
            file << std::setprecision(6) << std::fixed << E_kinetic[i] << delimiter;
            file << std::setprecision(6) << std::fixed << E_potential[i] << delimiter;
            file << std::setprecision(6) << std::fixed << Temp_sim[i] << delimiter;
            file << std::setprecision(6) << std::fixed << Pressure[i] << delimiter;
            file << "\n";
        }
        file.close();
    } else {
        std::cout << "Unable to open save file";
    }
}


void Simulation::plots() {
    // Energy vs time
    std::string energyPlot = "cd ../gnu-plots; gnuplot -p -c energy_time.gnuplot \'../" + m_foldername + "/" + m_filename + "_energy_time.png\' "
                                + "\'../" + m_foldername + "/" + m_filename + ".csv\' ";
    system(energyPlot.c_str());
    // Pressure vs time
    std::string pressurePlot = "cd ../gnu-plots; gnuplot -p -c pressure_time.gnuplot \'../" + m_foldername + "/" + m_filename + "_pressure_time.png\' "
                                + "\'../" + m_foldername + "/" + m_filename + ".csv\' ";
    system(pressurePlot.c_str());
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
    // Compress .exyz file
    std::string compress = "cd ../" + m_foldername + "; gzip " + m_filename + ".exyz";
    system(compress.c_str());


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
