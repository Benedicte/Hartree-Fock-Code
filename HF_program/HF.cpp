/**
    HF.cpp
    Purpose: Programming a very basic HF program

    @author Benedicte Ofstad
    @version 1.1 01/11/16
*/

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;
 
/**
    @param 	bond_length the bond length og the diatomic.
			zeta_atom_1 The slater orbital exsponent for function 1.
			zeta_atom_2 The slater orbital exponent for function 2.
			atom_1_atomic_number The atomic number of atom 1
			atom_2_atomic_number The atomic number of atom 2
    @return The one and two electron integrals
*/     

double one_and_two_integrals_calculation(int STO_n, double bond_length, double zeta_atom_1,
		double zeta_atom_2, double atom_1_atomic_number, 
		double atom_2_atomic_number)
{
		
		return 0;
}

double initial_fock_matrix() //This calculates the Fock Matrix
{
	return 0;
}

double error_function() //calulcates the error function
{
	return 0;
}

double overlap_integral() //calculated the overlap integral
{
	return 0;
}

double kinetic_energy_integral() //calculates the kinetic energy integral
{
	return 0;
}

double nuclear_attraction_integral() //calculates the nuclear attraction
{
	return 0;
}

double two_electron_integrals() //calcilates the two electron integrals
{
	return 0;
}

double assemble_integrals() // Assembles all the integrals
{
	return 0;
}

double SCF() //performs the SCF iterations
{
	return 0;
}

double create_G_matrix() //creates the G matrix
{
	return 0;
}

/**
    Function that returns the HF calcuation for a two electron 
    diatomic. Calls all the subroutines needed to calculate the
     Hartree Fock calculations

    @param 	bond_length the bond length og the diatomic.
			zeta_atom_1 The slater orbital exponent for function 1.
			zeta_atom_2 The slater orbital exponent for function 2.
			atom_1_atomic_number The atomic number of atom 1
			atom_2_atomic_number The atomic number of atom 2
    @return A hartree fock matrix?
*/

double hartree_fock_calculation(int STOn, double bond_length, double zeta_atom_1,
		double zeta_atom_2, double atom_1_atomic_number, 
		double atom_2_atomic_number)
{
	int not_sure_what_is_returned = one_and_two_integrals_calculation(STOn,
		bond_length, zeta_atom_1, zeta_atom_2, atom_1_atomic_number, 
		atom_2_atomic_number);
		
	//int not_sure_what_is_returned1 = assemble_integrals( 
	//	bond_length, zeta_atom_1, zeta_atom_2, atom_1_atomic_number, 
	//	atom_2_atomic_number);
		
	//int not_sure_what_is_returned2 = SCF( 
	//	bond_length, zeta_atom_1, zeta_atom_2, atom_1_atomic_number, 
	//	atom_2_atomic_number);
		
	double test_number = bond_length + zeta_atom_1;
	return test_number;
	
}

int main()
{
	int STOn = 1;
	double bond_length = 1.4632;
	double zeta_atom_1 = 2.0925;
	double zeta_atom_2 = 1.24; 
	double atom_1_atomic_number = 2; 
	double atom_2_atomic_number = 1;
	mat A = eye<mat>(5,5);
	
		
	double x = hartree_fock_calculation(STOn, bond_length, zeta_atom_1,
				zeta_atom_2, atom_1_atomic_number, atom_2_atomic_number);
	
	cout << x << endl;
	cout << "Program finished!" << endl;
	
	return 0;
}

