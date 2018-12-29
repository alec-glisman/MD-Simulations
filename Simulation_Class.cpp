//
// Created by Alec Glisman on 2018-12-28.
//

#include <iostream>
#include <vector>


class Simulation
{

public:
    // Public variables and methods


private:
    // Log other filesystem information
    std::string m_folder = "Q3.1-Energy";
    std::string m_filename = "sim00";
    // Simulation parameters that must be input by the constructor or will take these default values
    unsigned long m_n_dimensions = 3;
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
    std::vector<std::vector<double> > radii {m_n_particle, std::vector<double>(m_n_dimensions)};
    std::vector<std::vector<double> > velocities {m_n_particle, std::vector<double>(m_n_dimensions)};
    // Dynamic Variables

    // Trivial Default Constructor
public:
    Simulation() = default;


    // Print values to test
    void toPrint( std::vector<std::vector<double> > &vector)
    {
        // Loop through all rows
        for ( std::vector<std::vector<int>>::size_type i = 0; i < vector.size(); i++ )
        {
            // Loop through all columns (in each row)
            for ( std::vector<int>::size_type j = 0; j < vector[i].size(); j++ )
            {
                std::cout << vector[i][j] << "\t"; // Print the value at that location (i,j)
            }
            std::cout << "\n"; // End-line character for next line
        }
    }


};
