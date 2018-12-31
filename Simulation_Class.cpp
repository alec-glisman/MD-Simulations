//
// Created by Alec Glisman on 2018-12-28.
//

// Dependencies
#include <cmath> // pow, ceil
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

    // TODO: Energetic Variables



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
            m_filename(std::move(filename)),
            m_foldername(std::move(foldername)) {

        // Initialize radii with FCC lattice positions
        fcc_lattice_init();


    }


private:
    void fcc_lattice_init() {
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



        // initial FCC positions

    }

private:
    intMatrix_t cartesianProduct3(intVector_t &vector) {

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

private:
    intVector_t rangeVect(int start, int stop) {
        // Get range variable to iterate through
        intVector_t range;
        boost::push_back(range, boost::irange(start, stop));
        //toPrint(range, "Cell Range:");
        return range;
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
    void toPrint(intMatrix_t &vector, const string_t &string)
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




/*
      // TODO: TESTS (Delete later)
      toPrint(velocities, "Init Velocities:");
       */


/*
// TODO: rangeVect() TESTS (Delete later)
intVector_t rn1 = rangeVect(0, cells);
toPrint(rn1, "0 -> Cells:");
intVector_t rn2 = rangeVect(1, cells);
toPrint(rn2, "1 -> Cells:");
intVector_t rn3 = rangeVect(7, 20);
toPrint(rn3, "7 -> 20");
intVector_t rn4 = rangeVect(1, 2);
toPrint(rn4, "1 -> 2");
*/





/*
   // TODO: cartesianProduct3() TESTS (Delete later)
   intMatrix_t cartProd1 = cartesianProduct3(range);
   toPrint(cartProd1, "range - Cart Prod:");

   intVector_t v2 {0, 1};
   intMatrix_t cartProd2 = cartesianProduct3(v2);
   toPrint(cartProd2, "1 - Cart Prod:");

   intVector_t v3 {0};
   intMatrix_t cartProd3 = cartesianProduct3(v3);
   toPrint(cartProd3, "0 - Cart Prod:");

   intVector_t v4 {0, 1, 2, 3, 4, 5};
   intMatrix_t cartProd4 = cartesianProduct3(v4);
   toPrint(cartProd4, "5 - Cart Prod:");
   */


