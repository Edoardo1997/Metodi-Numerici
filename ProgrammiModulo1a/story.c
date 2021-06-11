//
//  ising_project.c
//  ising
//
//  
//

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
int field[1000][1000], npp[1000], nmm[1000];
int iflag;          //partenza caldo(1)/freddo(0)/precedente(altro)
int measures ;      //numero di step di questo run
int i_decorrel;     //updating fra una misura e l'altra
double extfield;      //valore del campo magnetico esterno
double beta;           //valore di 1/(kT) = beta
int nlatt;     //grandezza reticolo
long int seed;// seed da dare al rand di C. Occhio che non va bene farlo con il timeNULL) altrimenti per simulazioni troppo brevi non si campiona bene (stesso seed)

FILE *finput;
FILE *flattice;
FILE *fmisure;


/*
 Prototipi delle funzioni che saranno utilizzate nel codice
 */
void inizializza_field();
void crea_lattice();
double magnetizzazion();
double energy();
void geometry();
void update_metropolis();

int main()
{
    long int time_inizio;
    time_inizio = time(NULL);
   //leggo file
    
    finput = fopen("input.txt", "r");
    if (finput==NULL)
    {
        perror("Errore in apertura del file");
        exit(1);
    }
    
    flattice = fopen("input_lattice.txt", "r");
    if (flattice==NULL)
    {
        perror("Errore in apertura del file");
        exit(1);
    }

    fmisure = fopen("misure.txt", "w");
    if (fmisure==NULL)
    {
        perror("Errore in apertura del file");
        exit(1);
    }
    
    
 
   //leggo gli input
    fscanf(finput,"%ld\n" ,&seed);
    fscanf(finput,"%d\n" ,&nlatt);
    fscanf(finput,"%lf\n" ,&beta);
    fscanf(finput,"%lf\n" ,&extfield);
    fscanf(finput,"%d\n" ,&measures);
    fscanf(finput,"%d\n" ,&i_decorrel);
    fscanf(finput,"%d\n" ,&iflag);
    
    
    srand(seed);
    geometry();
    
    inizializza_field();
    for(int i=0;i<nlatt*nlatt;i++)//serve per la termalizzazione
        update_metropolis();
    int step;
    long double magn, ene;
    for(step=0;step<measures;step++)
    {
        for(int i=0;i<i_decorrel;i++)
            update_metropolis();
        magn=magnetizzazion();
        ene=energy();
        fprintf(fmisure,"%Lf %Lf \n",magn, ene);
    }
    
   
    
    printf("Il programma per beta %f e L %d ha runnato in %ld secondi\n",beta ,nlatt, time(NULL)-time_inizio);
    
    // chiudo i file aperti
    fclose(finput);
    fclose(flattice);
    fclose(fmisure);
  
    crea_lattice();
}



/*
 void inizializza_field() crea il lattice nxn in cui situiamo in ogni sito field[i][j] uno spin, 1 se è up o -1 se è down.
 L'obiettivo è inizializzarlo a uno stato iniziale.
 
 *iflag indica la tipica configurazione a determinate temperature.
 Il CASO [1] corrisponde alla temperatura nulla, in cui tutti gli spin o sono up o sono down. L'energia del mio sistema è invariante per simmetria di scambio, dunque è uguale porlo up o down, in questo caso ho scelto up.
 Il CASO [2] è tipico di una temperatura molto grande, maggiore della temperatura critica, in cui vi è disordine nel reticolo.
 Il CASO [3]
 
 
 */

void inizializza_field(){
    if(iflag==0)                            // CASO [1]
    {
        for(int i=0;i<nlatt;i++)
            for(int j=0;j<nlatt;j++)
                field[i][j]=1;
    }
    else if(iflag==1)                       //CASO [2]
    {
        for(int i=0;i<nlatt;i++)
            for(int j=0;j<nlatt;j++)
                field[i][j]=-1+2*(rand()%2);
    }
    else                                    //CASO [3]
    {
        for(int i=0;i<nlatt;i++)
           for(int j=0;j<nlatt;j++)
               fscanf(flattice,"%d\n" ,&field[i][j]);
    }
}

double magnetizzazion(){
    double sum=0;
    for(int i=0;i<nlatt;i++)
        for(int j=0;j<nlatt;j++)
            sum+=field[i][j];
    
    return sum/(nlatt*nlatt);
}

void geometry(){
    for(int i=0; i<nlatt;i++)
    {
        npp[i]=i+1;
        nmm[i]=i-1;
    }
    npp[nlatt-1]= 0;
    nmm[0]= nlatt-1;
}

double energy(){
    int jp, jm,ip, im;
    double xene = 0, force, nvol;
    for(int i=0;i<nlatt;i++)
        for(int j=0;j<nlatt;j++)
        {
            ip = npp[i];
            im = nmm[i];
            jp = npp[j];
            jm = nmm[j];
            force= field[i][jm] + field[i][jp] + field[im][j] + field[ip][j];
            xene = xene - 0.5*force*field[i][j];
            xene = xene - extfield*field[i][j];
            
        }
    nvol = nlatt*nlatt;
    xene = (double)xene/nvol;
    return xene;
}

void crea_lattice(){
    flattice = fopen("input_lattice.txt", "w");
    if (finput==NULL)
    {
        perror("Errore in apertura del file");
        exit(1);
    }
    for(int i=0;i<nlatt;i++)
       for(int j=0;j<nlatt;j++)
           fprintf(flattice,"%d\n",field[i][j]);
}


void update_metropolis(){
    int i,j;
    int jp, jm,ip, im;
    double force,p_rat, x;
    for(int ivol=0;ivol<nlatt*nlatt;ivol++)
    {
        i=rand()%nlatt;
        j=rand()%nlatt;
        
        ip = npp[i];
        im = nmm[i];
        jp = npp[j];
        jm = nmm[j];
        
        force= field[i][jm] + field[i][jp] + field[im][j] + field[ip][j];
        force=beta*(force+extfield);
        p_rat=exp(-2*field[i][j]*force);
        
        x=rand()/(RAND_MAX+1.0);
        if(x<p_rat)
            field[i][j]=-field[i][j];
    }
    
}

