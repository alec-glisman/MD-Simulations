//
// Created by Alec Glisman on 2018-12-31.
//

// Internal Project Dependencies
#include "Utilities.h"

// External Dependencies
#include <cmath> // pow, ceil

// TYPE ALIASES
using doubleMatrix_t =  std::vector<std::vector<double> >;
using doubleVector_t =  std::vector<double>;
using intMatrix_t =  std::vector<std::vector<int> >;
using intVector_t =  std::vector<int>;
using string_t = std::string;

// TODO: Make printing method a function template


void Utilities::print(doubleMatrix_t &vector, const string_t &string)
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


void Utilities::print(intMatrix_t &vector, const string_t &string)
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


void Utilities::print(intVector_t &vector, const string_t &string) {
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

