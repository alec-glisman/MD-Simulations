//
// Created by Alec Glisman on 2018-12-28.
//

// Dependencies
#include <cmath> // pow
#include<numeric> // iota
#include <iostream>
#include <vector>

#include <boost/range/irange.hpp> // i range
#include <boost/range/algorithm_ext/push_back.hpp> // push_back


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
               double sigma, double dt,
               string_t filename, string_t foldername) :
            m_num_iter(num_iter), m_num_simulations(num_simulations),
            m_temp(temp), m_n_particle(n_particle),
            m_box(box), m_epsilon(epsilon),
            m_sigma(sigma), m_dt(dt),
            m_t_total(dt * num_iter),
            m_vol(box * box * box), m_rho(n_particle / (box * box * box)),
            m_filename(std::move(filename)), m_foldername(std::move(foldername)) {
        fcc_lattice_init(); // Initialize radii

        toPrint(velocities, "Init Velocities:");
    }


private:
    void fcc_lattice_init() {
        double base = static_cast<double>(m_n_particle) / 4.0;
        double exponent = 1.0 / 3.0;
        auto cells = static_cast<int>(ceil(pow(base, exponent)));
        double cell_size = m_box / cells;

        doubleMatrix_t radius_{m_n_particle, doubleVector_t(m_n_dimensions)};


        cartesianProductFCC(radius_, cells);
        // TODO: Cartesian product and initial FCC positions

    }

private:
    doubleMatrix_t cartesianProductFCC(doubleMatrix_t radius_, int cells) {
        doubleMatrix_t r_fcc{{0.25, 0.25, 0.25},
                             {0.25, 0.75, 0.75},
                             {0.75, 0.75, 0.25},
                             {0.75, 0.25, 0.75}};

        std::vector<std::pair<double, double> > cartProd;

        // Get range variable to iterate through
        intVector_t cell_range;
        boost::push_back(cell_range, boost::irange(0, cells));
        toPrint(cell_range, "Cell Range:");

        // Create cartesian product with 3 repeats

        // TODO: Fill out function body


        return radius_;

    }


    // Print values to test
private:
    void toPrint(doubleMatrix_t &vector, const string_t &string)
    /// Method to print out a 2D Matrix
    /// \param vector a 2D vector with type double entries
    {
        std::cout << string + "\n";
        // Loop through all rows
        std::cout << "[";

        for (intMatrix_t::size_type i = 0; i < vector.size(); i++) {
            std::cout << "[";

            // Loop through all columns (in each row)
            for (intVector_t::size_type j = 0; j < vector[i].size(); j++) {
                if (j == 0 || j == 1) {
                    std::cout << vector[i][j] << ",\t"; // Print the value at that location (i,j)
                } else {
                    std::cout << vector[i][j] << "]";
                }
            }

            if (i < vector.size() - 1) {
                std::cout << "\n";  // End-line character for next line
            }
        }
        std::cout << "]\n\n";
    }

private:
    void toPrint(intVector_t &vector, const string_t &string) {
        std::cout << string + "\n";
        std::cout << "[";
        for (auto i = vector.begin(); i != vector.end(); ++i) {
            std::cout << *i;
            if (i < vector.end() - 1) {
                std::cout << ",\t";
            }
        }
        std::cout << "]\n\n";
    }

};
