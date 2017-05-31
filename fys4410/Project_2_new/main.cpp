#include "basis.h"
#include "hartree_fock.h"
#include "coulomb_functions.h"
#include "ccd.h"
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

    Basis first(shells);
    mat mapping = first.map_quantum_numbers(first.number_of_states);
    Hartree_fock_equations trial_1(shells, particles);
    trial_1.GetCoulombIntegrals(mapping, hw);
    mat fock_matrix = trial_1.hartree_fock_method(mapping);
    vec TBME = trial_1.new_basis(mapping, fock_matrix);
    double total_energy = trial_1.total_energy(mapping, trial_1.energies, TBME);
    double total_energy_1 = trial_1.total_energy_test(mapping, trial_1.energies, TBME);


    vec SP_energies = trial_1.energies;

    cout << "Single particle energies" << size(SP_energies) << endl;

    CCD trial(shells, particles, TBME, SP_energies);
    trial.intial_amplitudes();
    trial.CCD_solver();


    return 0;
}
