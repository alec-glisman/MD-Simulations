//
// Created by Alec Glisman on 2018-12-28.
//

#ifndef MD_SIMULATIONS_SIMULATION_CLASS_H
#define MD_SIMULATIONS_SIMULATION_CLASS_H

// Internal Project Dependencies
#include "Utilities.h"

// External Dependencies
#include <cmath> // pow, ceil
#include <iostream>
#include <vector>
#include <boost/range/irange.hpp> // i range
#include <boost/range/algorithm_ext/push_back.hpp> // push_back
#include <random>  // Normal distribution for velocity init


// TYPE ALIASES
using doubleMatrix_t =  std::vector<std::vector<double> >;
using doubleVector_t =  std::vector<double>;
using intMatrix_t =  std::vector<std::vector<int> >;
using intVector_t =  std::vector<int>;
using string_t = std::string;


class Simulation {

private:
    // Log and other filesystem information
    string_t m_foldername = "Q3.1-Energy";
    string_t m_filename = "sim00";

    // Simulation parameters that must be input by the constructor or will take these default values
    int m_num_iter = 4000;
    int m_num_simulations = 6;
    unsigned long m_n_particle = 50;
    double m_temp = 1.8;
    double m_box = 20.0;
    double m_epsilon = 1.0;
    double m_sigma = 1.0;
    double m_dt = 0.005;
    double m_t_total = 0.8;  // Overwritten later

    unsigned long m_n_dimensions = 3;
    long m_DoF = m_n_dimensions * (m_n_particle - 1);

    // Constant box properties
    double m_vol = m_box * m_box * m_box;
    double m_rho = m_n_particle / m_vol;

    // Kinematic Variables
    doubleMatrix_t radii{m_n_particle, doubleVector_t(m_n_dimensions)};
    doubleMatrix_t velocities{m_n_particle, doubleVector_t(m_n_dimensions)};

    // Dynamic Variables
    doubleMatrix_t forces{m_n_particle, doubleVector_t(m_n_dimensions)};

    // TODO: Energetic Variables


public:  // Default trivial constructor
    Simulation();

public: // Constructor with ``simulation variables''
    Simulation(int num_iter, int num_simulations,
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

};


#endif //MD_SIMULATIONS_SIMULATION_CLASS_H
