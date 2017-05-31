#include "basis.h"
#include "hartree_fock_equations.h"
#include "coulomb_functions.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
#include <typeinfo>
#include <chrono>

using namespace std;
using namespace arma;

typedef std::chrono::high_resolution_clock Clock;

int main()
{
    int particles = 2;
    int shells = 3;
    double hw = 1;

    double energy;

    auto t1 = Clock::now();

    Basis first(shells);
    mat mapping = first.map_quantum_numbers(first.number_of_particles);
    Hartree_fock_equations trial_1(shells, particles);
    trial_1.GetCoulombIntegrals(mapping, hw);
    vec energies = trial_1.hartree_fock_method(mapping);


    energy = trial_1.total_energy(mapping,energies);

    auto t2 = Clock::now();

    ofstream myfile;
    myfile.open ("/Users/benedicte/Programs/fys4410/Hartree_fock/build-Hartree_Fock-Desktop-Release/oooo1o.txt", ios::out);


    if (myfile.is_open()) {
        cout << "Wrote to file" << endl;
        myfile << "Number of particles: " << particles << "   Number of shells: " << shells <<"\n";
        myfile << "Ground state energy: " << std::fixed << std::setprecision(8) << energy <<"\n \n";

    }else {
            cout << "nei" << endl;
    }

    myfile.close();

    cout << "Program took "
        << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
        << " seconds" << std::endl;

    return 0;
}
