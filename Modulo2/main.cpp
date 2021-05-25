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
//UNDERSTAND : CREDO PRENDA DATI DA UN FILE IN CUI SI ALTERNANO RIGHE DI COMMENTI E RIGHE DI INPUT (BUTTA VIA I COMMENTI)


/*
 #define get_data(instream, strstream, str){ \
 getline(instream, str); \
 getline(instream, str); \
 strstream = (istringstream) str; \
 }
*/


using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    string run;
    if(argc > 1)
        run.assign(argv[1]);
    else
        run.assign("1");
    fstream inputstream("/Users/edoardo/Desktop/run"+run+".txt", ios_base::in); //replace with your path
    //more in general it would be helpful to understand how to make this part portable!
    //UNDERSTAND : TAKE INPUT DATA FROM DESKTOP
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
        int a, index=0;
        size_t s = count(line.begin(), line.end(),' ')+1; //conta quanti dati deve estrarre (contando gli spazi)
        int T[s];
        while(iss >> a){
            T[index]=a;
            index++;
        }
        N=zeros<Col<int> >(s); //Col è una template da Armadillo (dense column vector)
        for(int i=0; i<s ; i++){
            N(i)=T[i]; // ha creato un array con tutti gli N (lunghezza della catena) che voglio testare
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
    
    //cout << num_datas << endl; // debug
    
    //int n = 6;  // 12 take order of few seconds
    //bool PBC = true;
    //Mat<int> A = ML::codestate(n);
    
    //double m = ML::magnetization(n, 0.0, 0.3, A, PBC=true);
    //vec gaps = ML::gap(n, 0.0, 1.1, A, PBC=true);
    //cout << "magnetizazion : " << m << endl;
    
    //cout << "gaps : \n" << gaps << endl;
    
    
    fstream outputstream("/Users/edoardo/Desktop/University/Metodi_numerici/Modulo2/Progetto_modulo2/Exact_diagonalization/Exact_diagonalization/results"+run+".txt", ofstream::out | ofstream::trunc); //replace with your path
    //more in general it would be helpful to understand how to make this part portable!
    outputstream << "# h    #|M_z|"<< endl;
    /*
    for(int i=0; i< H.n_elem; i++){
        outputstream << H(i);
        for(int n=2; n<13; n++){
            Mat<int> A = ML::codestate(n);
            double mz = ML::magnetization(n, H(i), 1.0, A, true);
            outputstream << "    " << mz;
            //sp_mat sp_H = ML::sp_ham(n,0.0,G(i),A,true);
            //ML::print_matrix(sp_H,pow(2,n));
            //cout << "%%%%%%%%%%%%%%%%%%%%%%" << endl;
        }
        outputstream << endl;
    }*/
    //Mat<int> A = ML::codestate(4);
    //cout << A << endl;
    //cout <<ML::magnetization(4, 0.0, 0.5, A, true) << endl;
    /*
    vec H1;
    for(int n=2; n<13; n++){
        H1=linspace<vec>(-1,1, 80);
        Mat<int> A = ML::codestate(n);
        for(int i=0; i<80; i++){
            //cout << i << "    " << G1(i) << endl;
            outputstream << n << "    " << H1(i);
            double mz = ML::magnetization(n, H1(i), 0.9, A, true);
            //double chi = ML::chi_z_h(n, 0.0, G1(i), A, true);
            //outputstream << "    " << chi;
            outputstream << "    " << mz;
            //vec gaps(2);
            //gaps = ML::gap(n, 0.0, G1(i), A, true);
            //outputstream << "    " << gaps(0) << "    " << gaps(1);
            outputstream << endl;
        }
    }
    */
    
   /*
    vec G1;
    for(int n=2; n<18; n++){
        G1=linspace<vec>(0,2, 40);
        Mat<int> A = ML::codestate(n);
        for(int i=0; i<40; i++){
            //cout << i << "    " << G1(i) << endl;
            outputstream << n << "    " << G1(i);
            //eigs_sym(eigenvalues, eigenvectors, sp_H, 2, "sa"); //when failing give bool false as a result
            
            double mz = ML::s_magnetization(n, 0.0, G1(i), A, true);
            //double chi = ML::chi_z_h(n, 0.0, G1(i), A, true);
            //outputstream << "    " << chi;
            outputstream << "    " << mz;
            //vec gaps(2);
            //gaps = ML::gap(n, 0.0, G1(i), A, true);
            //outputstream << "    " << gaps(0) << "    " << gaps(1);
            outputstream << endl;
        }
    }
    
    */
    vec H1;
    for(int n=2; n<13; n++){
        H1=linspace<vec>(-1/pow(n,15./8.),1/pow(n,15./8.), 50);
        Mat<int> A = ML::codestate(n);
        for(int i=0; i<50; i++){
            //cout << i << "    " << G1(i) << endl;
            outputstream << n << "    " << H1(i);
            //eigs_sym(eigenvalues, eigenvectors, sp_H, 2, "sa"); //when failing give bool false as a result
            
            //double mz = ML::s_magnetization(n, 0.0, G1(i), A, true);
            double mz = ML::chi_z_h(n, H1(i), 1.0, A, true);
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
    /*
     vec H1,G1;
     H1=linspace<vec>(-1, 1, 100);
     G1=linspace<vec>(0, 2, 100);
     Mat<int> A = ML::codestate(12);
     for(int i=0; i<100; i++){
     for(int j=0; j<100; j++){
     outputstream << H1(i) << "   " << G1(j);
     double mz = ML::s_magnetization(12, H1(i), G1(j), A, true);
     outputstream << "    " << mz;
     outputstream << endl;
     }
     }
     
     */
    
    
    /*
     %%%%%%%%%%%%%%%%%%%%%%%%%%%    moved for debug      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for(int i=0; i<N.n_elem; i++){
        Mat<int> States = ML::codestate(N(i));
        if(G.n_elem == 1){
            for(int j=0; j<H.n_elem; j++){
                ML::all(N(i), H(j), G(0), States, PBC_e, k, outputstream);
            }
        }
        else{
            for(int j=0; j<G.n_elem; j++){
                ML::all(N(i), H(0), G(j), States, PBC_e, k, outputstream);
            }
        }
    }
     %%%%%%%%%%%%%%%%%%%%%%%%%%%    moved for debug      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     */
    /* //single line evaluator
    outputstream << "# h    #|M_z|"<< endl;
    for(int i=0; i<repetition; i++){
        outputstream << G(i);
        for(int n=2; n<10; n++){
            Mat<int> A = ML::codestate(n);
            double mz = ML::magnetization(n, H(i), 1.0, A, PBC=true);
            outputstream << "    " << mz;
        }
        outputstream << endl;
    }
     */
    
    /*
    outputstream << "# h    # g    #|M_z|"<< endl;
    for(int i=0; i<repetition; i++){
        for(int j=0; j<repetition; j++){
            outputstream << H(i) << "   " << G(j);
            double mz = ML::magnetization(n, H(i), G(j), A, PBC=true);
            outputstream << "    " << mz;
            outputstream << endl;
        }
    }
     */
    
    
    // outputstream << "# g    #|dE|"<< endl;
    
    /*
    for(int i=0; i<repetition; i++){
        outputstream << G(i);
        for(int n=2; n<10; n++){
            Mat<int> A = ML::codestate(n);
            vec gaps = ML::gap(n, 0.0, G(i), A, PBC=true);
            outputstream << "    " << gaps(0) << "    " << gaps(1);
        }
        outputstream << endl;
    }
     */
    
    /*
    for(int i=0; i<repetition; i++){
        outputstream << G(i);
        for(int n=2; n<10; n++){
            Mat<int> A = ML::codestate(n);
            vec ener = ML::energ(n, 0.0, G(i), A, PBC=true);
            outputstream << "    " << ener(0) << "    " << ener(1) << "    " << ener(2) << "    " << ener(3);
        }
        outputstream << endl;
    }
    */
    
    //ML::all(n, 0.2, 0.3, A, true, 5, outputstream);
    //ML::all(n, 0.4, 1.7, A, true, 5, outputstream);
    //outputstream.close(); //DEBUG
    return 0;
}

//Might be usefull in future, keep safe
/*
 if(inputstream.is_open()){
 while(getline(inputstream, line)){
 line.erase(find(line.begin(),line.end(),'#'),line.end()); //ignore comment starting with #
 //It might have problem with a line starting with #, check it out!
 istringstream iss(line); // create a stream from a string
 if(find(line.begin(),line.end(),'.') != line.end()){ //handling double
 double a;
 while(iss >> a)
 cout << a << endl; //It will become an initialization instead of being printed
 }
 else{ //handling int
 int a;
 while(iss >> a)
 cout << a << endl; //It will become an initialization instead of being printed
 }
 }
 }
 else
 cout << "Input file has not ben opened properly" << endl;
 inputstream.close();
 */
