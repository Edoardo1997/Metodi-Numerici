#include <iostream>
#include <fstream>
#include <algorithm>
#define ARMA_DONT_USE_WRAPPER
#define  ARMA_USE_LAPACK
#include <armadillo>
#include "Mainlib.h"
#include <cstdlib>

#define get_data(instream, strstream, stri){ \
            strstream.clear();\
            getline(instream, stri); \
            getline(instream, stri); \
            strstream.str(stri); \
        }

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    string run;
    if(argc > 1)
        run.assign(argv[1]);
    else
        run.assign("1");
    fstream inputstream("/Users/edoardo/Desktop/run" + run + ".txt", ios_base::in); //replace with your path
    vec H,G;
    Col<int> N;
    bool PBC_e;
    bool lanczos;
    int k;
    
    string line;
    
    if(inputstream.is_open()){
        cout << "Please work" << endl;
        istringstream iss;
        get_data(inputstream, iss, line)
        int a, index = 0;
        size_t s = count(line.begin(), line.end(),' ') + 1; //conta quanti dati deve estrarre (contando gli spazi)
        int T[s];
        while(iss >> a){
            T[index] = a;
            index++;
        }
        N=zeros<Col<int>>(s); //Col è una template da Armadillo (dense column vector)
        for(int i=0; i<s ; i++){
            N(i) = T[i]; // ha creato un array con tutti gli N (lunghezza della catena) che voglio testare
        }
        get_data(inputstream, iss, line)
        bool fix_g = line == "true" ? true : false; // decide se scorro su g (o h)
        get_data(inputstream, iss, line)
        double low_l, high_l; // credo, ma vedere dopo (che siano i limiti per cui scorre g (h))
        iss >> low_l;
        iss >> high_l;
        get_data(inputstream, iss, line)
        int num_datas; // quanti dati voglio (?)(frammentazione intervallo g)
        iss >> num_datas;
        cout << num_datas << endl;
        get_data(inputstream, iss, line)
        PBC_e = line == "true" ? true : false; // periodic boundary conditions or not
        get_data(inputstream, iss, line)
        lanczos = line == "true" ? true : false;
        if(lanczos){
            get_data(inputstream, iss, line)
            iss >> k;
        }
        if(fix_g){
            H=linspace<vec>(low_l, high_l, num_datas);
            G=ones<vec>(1);
        }
        else{
            G=linspace<vec>(low_l, high_l, num_datas);
            H=zeros<vec>(1);
        }
    }
    else{
        cout << "Input file has not ben opened properly" << endl;
        exit(1);
    }
    inputstream.close();
    
    cout << "G : " << G << endl;
    cout << "H : " << H << endl;
    cout << "N : " << N << endl;
    cout << "PBC : " << PBC_e << endl;
    cout << "Lanczos : " << lanczos << endl;
    if(!lanczos)
        k=1; // just to make it work
    cout << "k : " << k << endl;
    
    fstream outputstream("/Users/edoardo/Desktop/University/Metodi_numerici/Modulo2/Progetto_modulo2/Exact_diagonalization/Exact_diagonalization/results"+run+".txt", ofstream::out | ofstream::trunc);
    outputstream << "# h    #|M_z|"<< endl;
   
    vec G1;
    for(int n=2; n<13; n++){
        G1 = linspace<vec>(-1/pow(n,15./8.), 1/pow(n, 15./8.), 50);
        Mat<int> A = ML::codestate(n);
        for(int i=0; i<50; i++){
            cout << i << "    " << G1(i) << endl;
            //outputstream << n << "    " << H1(i);
            //eigs_sym(eigenvalues, eigenvectors, sp_H, 2, "sa"); //when failing give bool false as a result
            
            double mz = ML::s_magnetization(n, 0.0, G1(i), A, true);
            //double mz = ML::chi_z_h(n, H1(i), 1.0, A, true);
            //double chi = ML::s_chi_z_h(n, 0.0, G1(i), A, true);
            //outputstream << "    " << chi;
            //outputstream << "    " << mz;
            outputstream << "    " << mz;
            //vec gaps(2);
            //gaps = ML::gap(n, 0.0, G1(i), A, true);
            //outputstream << "    " << gaps(0) << "    " << gaps(1);
            outputstream << endl;
        }
    }
    
    return 0;
}
