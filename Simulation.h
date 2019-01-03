//
// Created by Alec Glisman on 2018-12-28.
//

#ifndef MD_SIMULATIONS_SIMULATION_CLASS_H
#define MD_SIMULATIONS_SIMULATION_CLASS_H


// Internal Project Dependencies
#include "Utilities.h"
#include "ProgressBar.hpp"

// External Dependencies
#include <cmath> // pow, ceil
#include <iostream> // Printing to terminal
#include <iomanip> // Set precision of printing
#include <fstream> // for file writing
#include <sstream> // to convert int to string_t
#include <vector>
#include <boost/range/irange.hpp> // i range
#include <boost/range/algorithm_ext/push_back.hpp> // push_back
#include <random>  // Normal distribution for velocity init
#include <chrono>  // Timing program execution
// #include "python/Python.h"  // To use GSD Output function


// TYPE ALIASES
using doubleMatrix_t =  std::vector<std::vector<double> >;
using doubleVector_t =  std::vector<double>;
using intMatrix_t =  std::vector<std::vector<int> >;
using intVector_t =  std::vector<int>;
using string_t = std::string;


class Simulation {

public:
    string_t m_filename_xyz;
    string_t m_filename_log;
    string_t m_filename_csv;

private:
    // Log and other filesystem information
    const string_t m_foldername = "Q3.1-Energy";
    const string_t m_filename = "sim00";

    string_t m_atom = "C"; // Atomic species used in simulation
    const double m_diameter = 0.5; // Particle diameter

    // Simulation parameters that must be input by the constructor or will take these default values
    const long m_num_iter = 4000;
    const int m_num_simulations = 6;
    const unsigned long m_n_particle = 50;
    const double m_temp = 1.8;
    const double m_box = 1.0;
    const double m_epsilon = 1.0;
    const double m_sigma = 1.0;
    const double m_dt = 0.005;
    const double m_t_total = 0.8;  // Overwritten later

    const static unsigned long m_n_dimensions = 3; // Spatial dimensions in system
    long m_DoF = m_n_dimensions * (m_n_particle - 1);  // Degrees of freedom
    long m_iterationNumber = 0;  // Current iteration number that is incremented
    const int m_n_dump = 1; // Dump an XYZ frame every m_n_dump iterations

    // Constant box properties
    const double m_vol = m_box * m_box * m_box;
    const double m_rho = m_n_particle / m_vol;

    // Kinematic Variables
    doubleMatrix_t radii{m_n_particle, doubleVector_t(m_n_dimensions)};
    doubleMatrix_t unwrapped_radii{m_n_particle, doubleVector_t(m_n_dimensions)};
    doubleMatrix_t velocities{m_n_particle, doubleVector_t(m_n_dimensions)};

    // Dynamic Variables
    doubleMatrix_t forces{m_n_particle, doubleVector_t(m_n_dimensions)};

    // Thermodynamic Variables
    doubleVector_t W{};  // 3X Internal Virial Energy
    doubleVector_t Temp_sim{};
    doubleVector_t Pressure{};

    // Energetic Variables
    doubleVector_t E_potential{};
    doubleVector_t E_kinetic{};
    doubleVector_t E_total{};


public:  // Default trivial constructor
    Simulation();

public: // Constructor with ``simulation variables''
    Simulation(long num_iter, int num_simulations,
               unsigned long n_particle, double temp,
               double box, double epsilon,
               double sigma, double dt,
               string_t filename, string_t foldername);

private: // Initializes lattice (radii) with FCC configuration
    void fcc_lattice_init();

private: // Helper function that calculates triple Cartesian product
    intMatrix_t cartesianProduct3(intVector_t &vector);

private: // Helper function to calculate an integer range (vector) between start (inclusive) and stop (exclusive)
    intVector_t rangeVect(int start, int stop);

private: // Initializes the velocities to have 0 lin momentum and KE that agrees with Temp
    void velocity_init();

private: // Saves and extended XYZ file format to be read into OVITO
    void create_frame();

private: // Writes simulation variables to log file (.txt)
    void log_data();

private: // Calculate forces, U, W
    void forceAndEnergetics();

private: // Calculate KE
    void kineticEnergy();

private: // Calculate Total E
    void totalEnergy();

private: // Calculate T_sim
    void simulationTemperature();

private: // Verlet Integration Step Method (Velocity)
    void velocityVerlet();

private: // Perform complete NVE integration
    void LJ_sim();

private: // Calculate system pressure
    void pressureCalc();

private: // Function to save variables
    void saveVars();

private: // Plot values
    void plots();

public: // Run entire program
    void main();


};


#endif //MD_SIMULATIONS_SIMULATION_CLASS_H
