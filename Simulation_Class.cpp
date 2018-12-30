//
// Created by Alec Glisman on 2018-12-28.
//

#include <iostream>
#include <vector>

// TYPE ALIASES
using doubleMatrix_t =  std::vector<std::vector<double> >;
using doubleVector_t =  std::vector<double>;
using intMatrix_t =  std::vector<std::vector<int> >;
using intVector_t =  std::vector<int>;



class Simulation
{

private:

    // Log and other filesystem information
    std::string m_folder = "Q3.1-Energy";
    std::string m_filename = "sim00";

    unsigned long m_n_dimensions = 3;

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

    // Constant box properties
    double m_vol = m_box * m_box * m_box;
    double m_rho = m_n_particle / m_vol;

    // Kinematic Variables
    doubleMatrix_t radii{m_n_particle, doubleVector_t(m_n_dimensions)};
    doubleMatrix_t velocities{m_n_particle, doubleVector_t(m_n_dimensions)};

    // Dynamic Variables
    doubleMatrix_t forces{m_n_particle, doubleVector_t(m_n_dimensions)};

    // Energetic Variables



    // Trivial Default Constructor
public:
    Simulation() = default;


    // Constructor to input simulation parameters
public:
    Simulation(int num_iter, int num_simulations,
               unsigned long n_particle, double temp,
               double box, double epsilon,
               double sigma, double dt) :
            m_num_iter(num_iter), m_num_simulations(num_simulations),
            m_temp(temp), m_box(box),
            m_epsilon(epsilon), m_sigma(sigma),
            m_dt(dt), m_t_total(dt * m_num_iter) {
        // Nothing to do here
    }


    // Print values to test
private:
    void toPrint(doubleMatrix_t &vector)
    {
        // Loop through all rows
        for (intMatrix_t::size_type i = 0; i < vector.size(); i++)
        {
            // Loop through all columns (in each row)
            for (intVector_t::size_type j = 0; j < vector[i].size(); j++)
            {
                std::cout << vector[i][j] << "\t"; // Print the value at that location (i,j)
            }
            std::cout << "\n"; // End-line character for next line
        }
    }


};
