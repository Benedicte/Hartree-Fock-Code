#include "hartree_fock_equations.h"

Hartree_fock_equations::Hartree_fock_equations(int sh_number, int number_of_electrons){

    int shell_number = sh_number;
    this->number_of_electrons = number_of_electrons;

    for(int i = shell_number; i > 0; i--){
        number_of_particles = number_of_particles + i*2;
    }
}

void Hartree_fock_equations::GetCoulombIntegrals(mat mapping, double hw){

    coloumb_matrix_direct = new double[number_of_particles*number_of_particles*number_of_particles*number_of_particles];
    coloumb_matrix_exchange = new double[number_of_particles*number_of_particles*number_of_particles*number_of_particles];

    int n_alpha, m_alpha, s_alpha;
    int n_beta, m_beta, s_beta;
    int n_gamma, m_gamma, s_gamma;
    int n_delta, m_delta, s_delta;

    int np = number_of_particles;
    int np2 = np*np;
    int np3 = np2*np;

    for( int alpha = 0; alpha < number_of_particles; alpha++){

        //cout << "alpha "  << alpha << endl;

        n_alpha = mapping(alpha,0);
        m_alpha = mapping(alpha,1);
        s_alpha = mapping(alpha,2);

        for( int beta = 0; beta < number_of_particles; beta++){

            n_beta = mapping(beta,0);
            m_beta = mapping(beta,1);
            s_beta = mapping(beta,2);

            for (int gamma = 0; gamma < number_of_particles; gamma++){

                n_gamma = mapping(gamma,0);
                m_gamma = mapping(gamma,1);
                s_gamma = mapping(gamma,2);

                for( int delta = 0; delta < number_of_particles; delta++){

                   // cout << "inner "  << alpha << endl;

                    n_delta = mapping(delta,0);
                    m_delta = mapping(delta,1);
                    s_delta = mapping(delta,2);

                    if(s_gamma == s_delta && s_alpha == s_beta){

                        coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta] =
                                Coulomb_HO(hw, n_alpha, m_alpha, n_gamma, m_gamma, n_beta, m_beta, n_delta, m_delta);

                    } else{
                         coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta] = 0.0;
                        }


                    if(s_alpha == s_delta && s_gamma == s_beta){
                        coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta] =
                                Coulomb_HO(hw, n_alpha, m_alpha, n_gamma, m_gamma, n_delta, m_delta, n_beta, m_beta);
                    } else{
                        coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta] = 0.0;
                        }
                }
            }
        }
    }
}

vec Hartree_fock_equations::hartree_fock_method(mat mapping){

    cout << "number of particles " << number_of_particles << endl;
    slater_determinant = eye<mat>(number_of_particles, number_of_particles);
    density_matrix = zeros<mat>(number_of_particles, number_of_particles);
    mat fock_matrix = zeros<mat>(number_of_particles,number_of_particles);

    int np = number_of_particles;
    int np2 = np*np;
    int np3 = np2*np;

    double direct_elements = 0;
    double exchange_elements = 0;
    double hw = 1;

    int maxIterator = 200;
    double epsilon = 1.0e-9;
    double difference = 1.0;
    int hf_count = 1;
    double sumFockTerm = 0.0;

    int s_alpha, s_beta, s_gamma, s_delta;

    vec oldenergies(number_of_particles);
    vec newenergies(number_of_particles);

    for( int gamma = 0; gamma < number_of_particles; gamma++){
         for( int delta = 0; delta < number_of_particles; delta++){
            double sum = 0.0;
            for(int i = 0; i < number_of_electrons; i++){
                sum += slater_determinant(gamma,i)*slater_determinant(delta,i);
            }
            density_matrix(gamma, delta) = sum;
        }
    }

    while( hf_count < maxIterator && difference >= epsilon){
        cout << "############### Iteration " << hf_count <<" ###############" << endl;
        fock_matrix.zeros();

        for( int alpha = 0; alpha < number_of_particles; alpha++){

            double energy_alpha = Hartree_fock_equations::get_energy(mapping(alpha,0), mapping(alpha,1), hw);

            s_alpha = mapping(alpha,2);

            fock_matrix(alpha,alpha) = energy_alpha;

            for( int beta = 0; beta < number_of_particles; beta++){

                sumFockTerm = 0;
                s_beta = mapping(beta,2);

                for (int gamma = 0; gamma < number_of_particles; gamma++){

                    s_gamma = mapping(gamma,2);

                    for( int delta = 0; delta < number_of_particles; delta++){

                        s_delta = mapping(delta,2);
                        direct_elements = 0;
                        exchange_elements = 0;

                        if(s_gamma == s_delta && s_alpha == s_beta){
                            direct_elements = density_matrix(gamma,delta)
                                *coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                            sumFockTerm += direct_elements;
                        }

                        if(s_alpha == s_delta && s_gamma == s_beta){
                            exchange_elements = density_matrix(gamma,delta)
                                  *coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta];
                            sumFockTerm -= exchange_elements;
                        }
                    }               
                }
                //fock_matrix(alpha, beta) += sumFockTerm;
                 fock_matrix(beta, alpha) += sumFockTerm;

            }
        }

        vec energies;

        eig_sym(energies, slater_determinant, fock_matrix);

        /// Setting up new density matrix

        for( int gamma = 0; gamma < number_of_particles; gamma++){
            for( int delta = 0; delta < number_of_particles; delta++){
                double sum = 0.0;

                for(int i = 0; i < number_of_electrons; i++){
                    sum += slater_determinant(gamma,i)*slater_determinant(delta,i);
                }
                density_matrix(gamma, delta) = sum;
            }
        }

        newenergies = energies;

        //Brute force computation of difference between previous and new sp HF energies """
        double sum = 0.0;

        for( int i = 0; i < number_of_particles; i++){
            sum += (std::abs(newenergies(i)-oldenergies(i)))/number_of_particles;
            difference = sum;
            oldenergies = newenergies;
        }

        hf_count ++;
        cout << "The difference " << difference << endl;
    }

    std::setprecision(7);
    cout.setf(ios::fixed);

    fock_matrix.print();

    for(int i = 0; i < number_of_particles; i++){
        for(int j = 0; j < number_of_particles; j++){

            if( fabs(fock_matrix(i,j)) < 1E-8)
                cout << 0 << " ";
            else
                cout << fock_matrix(i,j) << " ";

        }
        cout << endl;
    }

    newenergies.raw_print(std::cout);

    return newenergies;
}

double Hartree_fock_equations::get_energy(int n, int m, double omega){

    double energy;

    energy = omega*(2.0*n + abs(m) + 1);

    return energy;
}

double Hartree_fock_equations::kroneckerDelta(int x, int y) {

  if(x==y) {
      return 1;
  } else {
      return 0;
  }
}

double Hartree_fock_equations::total_energy(mat mapping, vec energies){

    int fermi_level = number_of_electrons;
    double ground_state_energy;
    double single_particle_contributions = 0.0;
    double double_particle_contributions = 0.0;
    double direct_elements;
    double exchange_elements;
    double anti_sym;



    int np = number_of_particles;
    int np2 = np*np;
    int np3 = np2*np;


    for(int i = 0; i < fermi_level; i++){
        single_particle_contributions += energies(i);
    }

    for( int alpha = 0; alpha < number_of_particles; alpha++){
        int s_alpha = mapping(alpha,2);

        for( int beta = 0; beta < number_of_particles; beta++){
            int s_beta = mapping(beta,2);

            for (int gamma = 0; gamma < number_of_particles; gamma++){

                int s_gamma = mapping(gamma,2);

                for( int delta = 0; delta < number_of_particles; delta++){

                    int s_delta = mapping(delta,2);

                    direct_elements = 0.0;
                    exchange_elements = 0.0;

                    if(s_alpha == s_beta && s_gamma == s_delta){
                        direct_elements =
                                coloumb_matrix_direct[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                    }

                    if(s_alpha == s_delta && s_gamma == s_beta){

                    exchange_elements =
                            coloumb_matrix_exchange[np3*(alpha) + np2*(beta) + np*(gamma) + delta];

                    }

                    anti_sym = direct_elements - exchange_elements;

                    if(anti_sym != 0){
                    double_particle_contributions += density_matrix(alpha,beta)*density_matrix(gamma,delta)*anti_sym;
                    }
                }
            }
        }
    }


    //double_particle_contributions = 19.527
    //Ground state energy 20.76692

    ground_state_energy = single_particle_contributions - 0.5 * double_particle_contributions;

    cout << "number of particles: " << number_of_electrons << endl;
    cout << "Ground state energy" << endl;
    cout << ground_state_energy << endl;


    return ground_state_energy;
}

/*mat C_igamma = trans(slater_determinant);
mat C_gammai = C_igamma.submat( 0,number_of_particles, 0,number_of_electrons);
mat C_idelta = slater_determinant.submat(0, number_of_electrons, 0, number_of_particles);

density_matrix_new = C_gammai*C_idelta;

cout << density_matrix << "and the new one: " << density_matrix_new << endl;
*/


