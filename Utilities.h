//
// Created by Alec Glisman on 2018-12-31.
//

#ifndef MD_SIMULATIONS_UTILITIES_H
#define MD_SIMULATIONS_UTILITIES_H

#include <iostream>
#include <vector>
#include <iomanip> // Set precision of printing
#include <fstream> // for file writing

// TYPE ALIASES
using doubleMatrix_t =  std::vector<std::vector<double> >;
using doubleVector_t =  std::vector<double>;
using intMatrix_t =  std::vector<std::vector<int> >;
using intVector_t =  std::vector<int>;
using string_t = std::string;


class Utilities {

/*
public:
    template <typename Type>
    static const void print(const Type &vector, const string_t &string);
*/

public:
    static void print(doubleMatrix_t &vector, const string_t &string);

public:
    static void print(intMatrix_t &vector, const string_t &string);

public:
    static void print(intVector_t &vector, const string_t &string);

public:
    static void saveCSV(doubleVector_t &vector, const string_t &filename);

};


#endif //MD_SIMULATIONS_UTILITIES_H
