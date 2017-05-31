#ifndef BASIS_H
#define BASIS_H
#include <armadillo>
using namespace arma;

class Basis{

public:

    Basis (int sh_number);

    int shell_number;
    bool spin;
    int energy;
    int projection;
    int number_of_particles;

    void set_energy();
    mat get_quantum_numbers(int shell);
    mat map_quantum_numbers(int number_of_particles);
};

#endif // BASIS_H
