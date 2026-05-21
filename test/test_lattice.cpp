#include "Lattice.hpp"
#include <iostream>
#include <cassert>

/*
File: lattice_test.cpp
Description: Unit testing of implemnted methods.
Author: Marcel Wilanowicz
Date: 2026-04-01
*/

int main() {
    Lattice lattice; // Create an instance of the Lattice class

    // Index mapping tests
    std::cout << "Testing: index maping..." << std::endl;
    assert(lattice.map_idx(0, 0) == 0); // First element in first row
    assert(lattice.map_idx(1, 0) == 1); // second element in first row
    assert(lattice.map_idx(0, 1) == LBM::Config::width); // Start of the second row (after one full width step)
    
    // mappped(last element in x, last element in y) should be equal to last element in last row of mapped/linearized grid
    assert(lattice.map_idx(LBM::Config::width - 1, LBM::Config::height - 1) == (LBM::Config::width * LBM::Config::height) - 1); // Last element in last row
    std::cout << "Passed!" << std::endl;

    // Memory allocation tests
    std::cout << "Testing: memory allocation..." << std::endl;
    for (int i = 0; i < LBM::Q; ++i) { // All 9 directions
        assert(lattice.get_component_size(i) == static_cast<size_t>(LBM::Config::width * LBM::Config::height));
    }
    std::cout << "Passed!" << std::endl;

    // Double buffering test
    std::cout << "Testing: double buffering..." << std::endl;
    lattice.set_f_new(100, 50, 4, 1.2345); // Random value for this vertex
    lattice.swap();
    assert(lattice.get_f_old(100, 50, 4) == 1.2345); // Is the same after swap?
    std::cout << "Passed!" << std::endl;

    // Checking rho value
    std::cout << "Testing: rho value..." << std::endl;
    lattice.initialize();
    assert(lattice.get_rho(10, 10) - 1.0 < 1e-9);
    std::cout << "Passed!" << std::endl;

    // Checking ux, uy values 
    std::cout << "Testing: u value..." << std::endl;
    auto u = lattice.get_u(10, 10);
    assert(std::abs(u.first) < 1e-12);
    assert(std::abs(u.second) < 1e-12);
    std::cout << "Passed!" << std::endl;

    return 0;
}
