//
//  Mainlib.h
//  Exact_diagonalization
//
//  Created by Edoardo centamori on 28/10/2019.
//  Copyright © 2019 Edoardo centamori. All rights reserved.
//

#ifndef Mainlib_h
#define Mainlib_h

#include <iostream>
#include <cmath>
#define ARMA_DONT_USE_WRAPPER
#define  ARMA_USE_LAPACK
#include <armadillo>


using namespace std;
using namespace arma;

namespace ML {
    
    /**
     * Exact ground state energy for (or h=0 or g=1).
     * Used to check validity of the code.
     * @param n : dimension of the chain.
     * @return: analitical ground state energy.
     */
    double expected_en(int n){
        
        return 1.0-1.0/sin(M_PI/(2*(2*n+1)));
    }
    
    /**
     * Coding from integer (index) to states of the computational base
     * @param n : dimension of the chain.
     * @return Mat<int> (n,2^n): matrix of states encoded in binary notation
     */
    Mat<int> codestate(int n){
        
        int N = pow(2,n);
        Mat<int> A(N, n);
        for(int i=0; i<N; i++){
            for( int j=0, i_temp=i; j<n; j++){
                A(i, n-1-j) = i_temp % 2; // Professor uses n-1-j -> j, not sure why!
                i_temp /= 2;
            }
        }
        return A;
    }
    
    
    /**
     * Create system matrix
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return Mat<double> (2^n,2^n): Hamiltonian of the Ising model in the computational base
     */
    Mat<double> ham(int n, double h, double g, Mat<int> states, bool PBC){
        
        int N = pow(2,n);
        Mat<double> A(N,N, fill::zeros);
        //sigma_z * sigma_z coupling
        // notice that it's diagonal thus can affect only terms of A(i,i)!
        for(int i_site=0; i_site<n-1; i_site++){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,i_site) == states(i_state,i_site+1))
                    A(i_state,i_state) -= 1.0;
                else
                    A(i_state,i_state) += 1.0;
            }
        }
        if(PBC){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,n-1) == states(i_state,0))
                    A(i_state, i_state) -= 1.0;
                else
                    A(i_state, i_state) += 1.0;
            }
        }
        //sigma_z
        for(int i_site=0; i_site<n; i_site++){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,i_site) == 1)
                    A(i_state,i_state) -= h;
                else
                    A(i_state,i_state) += h;
            }
        }
        //sigma_x
        int i_coupled; //determine wich state is coupled to i_state by sigma_x
        for(int i_site=0; i_site<n; i_site++){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,i_site) == 1)
                    i_coupled = i_state - pow(2, n-1-i_site);
                    // e.g. (i_site = 0) |100> (4) is coupled with |000> (0 = 4 - 2^(3-1-0))
                else
                    i_coupled = i_state + pow(2, n-1-i_site);
                    // e.g. (i_site = 0) |010> (2) is coupled with |110> (6 = 2 + 2^(3-1-0))
                A(i_coupled, i_state) += g;
            }
        }
        return A;
    }
    
    
    /**
     * Create system matrix in sparse form
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return Mat<double> (2^n,2^n): Hamiltonian of the Ising model in the computational base (sparse matrix form)
     */
    sp_mat sp_ham(int n, double h, double g, Mat<int> states, bool PBC){
        // PBC need to be implemented
        int N = pow(2,n);
        sp_mat A(N,N);
        //sigma_z * sigma_z coupling
        // notice that it's diagonal thus can affect only terms of A(i,i)!
        for(int i_site=0; i_site<n-1; i_site++){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state, i_site) == states(i_state, i_site+1))
                    A(i_state,i_state) -= 1.0;
                else
                    A(i_state,i_state) += 1.0;
            }
        }
        if(PBC){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,n-1)==states(i_state,0))
                    A(i_state,i_state) -= 1.0;
                else
                    A(i_state,i_state) += 1.0;
            }
        }
        //sigma_z
        for(int i_site=0; i_site<n; i_site++){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,i_site) == 1)
                    A(i_state,i_state) -= h;
                else
                    A(i_state,i_state) += h;
            }
        }
        //sigma_x
        int i_coupled; //determine wich state is coupled to i_state by sigma_x
        for(int i_site=0; i_site<n; i_site++){
            for(int i_state=0; i_state<N; i_state++){
                if(states(i_state,i_site) == 1)
                    i_coupled = i_state - pow(2, n-1-i_site);
                // e.g. (i_site = 0) |100> (4) is coupled with |000> (0 = 4 - 2^(3-1-0))
                else
                    i_coupled = i_state + pow(2, n-1-i_site);
                // e.g. (i_site = 0) |010> (2) is coupled with |110> (6 = 2 + 2^(3-1-0))
                A(i_coupled, i_state) += g;
            }
        }
        return A;
    }
    
    
    /**
     * Evaluate transverse magnetization
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return avarage magnetization of the chain in the z axis
     */
    double magnetization(int n, double h, double g, Mat<int> states, bool PBC){
        int N = pow(2,n);
        mat H=ham(n, h, g, states, PBC);
        vec eigenvalues;
        mat eigenvectors;
        eig_sym(eigenvalues,eigenvectors, H);
        vec ground = eigenvectors.col(0);
        vec mag(N);
        double magz = 0.0;
        for(int i_state = 0; i_state<N; i_state++){
            mag(i_state) = 0.0;
            for(int i_site = 0; i_site<n; i_site++){
                if(states(i_state,i_site)==1) mag(i_state)+= 1.0;
                else                          mag(i_state)-= 1.0;
            }
            /*
             In the next line:
             -use abs when working with h == 0 to deal with degeneration of ground state
             -use version without abs to deal with h != 0
             */
            //magz += abs(mag(i_state))*pow(abs(ground(i_state)),2);
            magz += mag(i_state)*pow(ground(i_state),2);
        }
        return magz/n; // '/n' provide the mean magnetization
    }
    
    /**
     * Evaluate longitudinal magnetization
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return avarage magnetization of the chain in the x axis
     */
    double x_magnetization(int n, double h, double g, Mat<int> states, bool PBC){
        int N = pow(2,n);
        mat H=ham(n, h, g, states, PBC);
        vec eigenvalues;
        mat eigenvectors;
        eig_sym(eigenvalues,eigenvectors, H);
        vec ground = eigenvectors.col(0);
        vec mag(n);
        int i_coupled;
        double magx = 0.0;
        for(int i_site=0; i_site<n; i_site++){
            mag(i_site)=0.0;
            for(int i_state=0; i_state<N;i_state++){
                if(states(i_state,i_site)==1)
                    i_coupled = i_state - pow(2,n-1-i_site);
                else
                    i_coupled = i_state + pow(2,n-1-i_site);
                mag(i_site)+=ground(i_state)*ground(i_coupled);
            }
            magx += mag(i_site);
        }
        return magx/n; // '/n' provide the mean magnetization
    }

    
    /**
     * Evaluate suscettibility.
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return avarage magnetic suscettibility (z axis)
     */
    double chi_z_h(int n, double h, double g, Mat<int> states, bool PBC){
        double delta = 0.001;
        double magz_m = magnetization(n, h-delta/2., g, states, PBC);
        double magz_p = magnetization(n, h+delta/2., g, states, PBC);
        return (magz_p-magz_m)/delta;
    }

    
    /**
     * Evaluate energy gaps
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return first and second energy gap of the system
     */
    vec gap(int n, double h, double g, Mat<int> states, bool PBC){
        mat H=ML::ham(n, h, g, states, PBC);
        vec eigenvalues;
        mat eigenvectors;
        eig_sym(eigenvalues,eigenvectors, H);
        vec gaps(2);
        gaps(0)=eigenvalues[1]-eigenvalues[0];
        gaps(1)=eigenvalues[2]-eigenvalues[0];
        return gaps;
    }
    
    
    /**
     * Evaluate low energies.
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return several energy levels
     */
    vec energ(int n, double h, double g, Mat<int> states, bool PBC){
        mat H=ML::ham(n, h, g, states, PBC);
        vec eigenvalues;
        mat eigenvectors;
        eig_sym(eigenvalues,eigenvectors, H);
        vec energ(4);
        energ(0)=eigenvalues[0];
        energ(1)=eigenvalues[1];
        energ(2)=eigenvalues[2];
        energ(3)=eigenvalues[3];
        return energ;
    }
    
    
    /**
     * Put many relvant information in the output file.
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @param k : number of relevant eigenstates of the system
     * @param outputstream : file stream to output file : it assume that the file is already open so do it outside
     * Side effect: print several usefull information in a file for external analysis in python
     */
    void all(int n, double h, double g, Mat<int> states, bool PBC, int k, fstream & outputstream){
        sp_mat sp_H = sp_ham(n,h,g,states,PBC);
        vec eigenvalues;
        mat eigenvectors;
        eigs_sym(eigenvalues, eigenvectors, sp_H, k);
        outputstream << n << "    " << h << "    " << g << "    " << eigenvalues.t();
        outputstream.seekp(-1,ios_base::cur); //remove the undesired '\n' coming from << eigenvalues.t()
        outputstream << "    " << eigenvectors.as_row() << endl;
    }
    
    
    /**
     * Evaluate transverse magnetization (Lanczos routine)
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return : Magnetization in the z direction evaluated using the Lanczos routine
     * Can handle errors from Lanczos
     */
    double s_magnetization(int n, double h, double g, Mat<int> states, bool PBC){
        int N = pow(2,n);
        sp_mat sp_H = sp_ham(n,h,g,states,PBC);
        vec eigenvalues;
        mat eigenvectors;
        if(eigs_sym(eigenvalues,eigenvectors, sp_H, 1,"sa")){
            vec ground = eigenvectors.col(0);
            vec mag(N);
            double magz = 0.0;
            for(int i_state = 0; i_state<N; i_state++){
                mag(i_state) = 0.0;
                for(int i_site = 0; i_site<n; i_site++){
                    if(states(i_state,i_site)==1) mag(i_state)+= 1.0;
                    else                          mag(i_state)-= 1.0;
                }
                /*
                In the next line:
                -use abs when working with h == 0 to deal with degeneration of ground state
                -use version without abs to deal with h != 0
                 */
                //magz += abs(mag(i_state))*pow(abs(ground(i_state)),2);}
                magz += (mag(i_state))*pow(ground(i_state),2);}
            return magz/n; // '/n' provide the mean magnetization
        }
        else
            return -10.0; //error, used to detect when eigs_sym fail for problems of Armadillo
    }
    
    
    /**
     * Evaluate transverse magnetization (Lanczos routine)
     * @param n : dimension of the chain.
     * @param h : external magnetic field
     * @param g : 'temperature'
     * @param states : output of codestate
     * @param PBC : Flag, if false the open boundary conditions, if true then periodic boundary conditions.
     * @return : Magnetic susceptibility (Lanczos)
     */
    double s_chi_z_h(int n, double h, double g, Mat<int> states, bool PBC){
        double delta = 0.001;
        double magz_m = s_magnetization(n, h-delta/2., g, states, PBC);
        double magz_p = s_magnetization(n, h+delta/2., g, states, PBC);
        return (magz_p-magz_m)/delta;
    }
    
    
    /**
     * Print matrix to debug
     * @param matrix : matrix.
     * @param dim : dimension of the matrix
     */
    void print_matrix(sp_mat matrix, int dim){
        for(int i=0; i<dim; i++){
            for(int j=0; j<dim; j++){
                cout << matrix(i,j) << " ";
            }
            cout << endl;
        }
        return;
    }
    
    
}


#endif /* Mainlib_h */
