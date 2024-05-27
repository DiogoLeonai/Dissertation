//Adaptive exponential integrate-and-fire
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
///////////parametro gerador aleatorio
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

///////Valores dos parâmetros/////////
#define pi 3.14
#define C  200.0
#define gl 12.0
#define El -70.0
#define Vt -50.0
#define deltat 2.0
#define tw 300.0
#define I 512.0
#define Vmax -40.0
#define taus 2.728
#define a 2.0

////////////////////////////////////////////////////////
double ran1(long *idum);
void nrerror(char error_text[]);
#define NR_END 1
#define FREE_ARG char*
/////// Ponteiros para alocar memoria ////////
int *vector(long nl,long nh);
double *dvector(long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_vector(int *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
//////////////////////////////////////////////////////

#define N 1 //número de neurônios

////parâmetros do runge-kutta
#define h 0.1 //passo
#define eq_num  2  //número de equações


double edo(double *df, int k, double variavel, double alpha, double beta);


FILE *yyy;

int main(){
  double x[eq_num+1], t, function[eq_num+1], coefrunge[eq_num+1][5], *IED, IEDmedia, desv, sum, CV, media, A, Vreset, b, *firing_rate, mediafiring, firing_media, t_ant, alpha, beta;
  int n, i, disp, escala_cv, sharp, broad, aux, aux2, *quantidade_sharp, k, teste, verificando, regular;

IED = dvector(0, 10010000+1);
firing_rate = dvector(0, 10010000+1);
quantidade_sharp = vector(0, 10010000+1);

  ////////////////////condições iniciais//////////////
 function[1]= -70.0;
 function[2]= 0.0; //corrente de adaptação
for(n=1; n<=10010000+1; n++){
            quantidade_sharp[n]=0;}
 
 disp=0;
 k=0;
 aux=0;
 regular = 0;
 verificando =0;
  yyy=fopen("ResetSpace_09.dat","wt");

  
  alpha = 0.9;
  beta = 0.9;
  
 for(Vreset = -62.05; Vreset<=-40.0; Vreset = Vreset + 0.15){
 for(b=0.0; b<=400; b=b+2.0){//espaco 200X200
   for(t=0.0001; t<=500000; t=t+h){
     
         
            
            ///////* Runge Kutta 4ª Ordem*////////   
            for(i=1; i<=eq_num; i++){
            coefrunge[i][1] = edo(function, i, t, alpha, beta);
            }
            
            
            for(i=1; i<=eq_num; i++){
            x[i] = function[i] + (h*coefrunge[i][1]*0.5);
            }
            for(i=1; i<=eq_num; i++){
            coefrunge[i][2] = edo(x, i, (t+h/2), alpha, beta);
            }
            
            
            for(i=1; i<=eq_num; i++){
            x[i] = function[i] + (h*coefrunge[i][2]*0.5);
            }
            for(i=1; i<=eq_num; i++){
            coefrunge[i][3] = edo(x, i, (t+h/2), alpha, beta);
            }
            
            
            for(i=1; i<=eq_num; i++){
            x[i] = function[i] + (h*coefrunge[i][3]);
            }
            for(i=1; i<=eq_num; i++){
            coefrunge[i][4] = edo(x, i, (t+h), alpha, beta);
            }
             
             for(i=1; i<=eq_num; i++){
             function[i] = function[i]+ (h/6)*(coefrunge[i][1]+2*coefrunge[i][2]+2*coefrunge[i][3]+coefrunge[i][4]);
             }
         
            /////condições de reinício///////
             if(function[1]>Vmax){
             
             disp++;
             function[1] = Vreset;
             function[2] = function[2]+ b;
             
             if(disp>1){
              IED[disp-1] = t - t_ant;
              firing_rate[disp-1] = 1/IED[disp-1];
              
              if(disp>8){
              media = media + IED[disp-1];
              mediafiring = mediafiring + firing_rate[disp-1];}}
             
             t_ant = t; // salvando o tempo do disparo anterior
             
            
             if(edo(function, 1, t, alpha, beta) < 0){
               broad++;
               
               if(aux2==1){
                k++;
                aux2=0;}
               
               aux = 1;}
             else{
               sharp++;
               if(aux==1){
                quantidade_sharp[broad] = quantidade_sharp[broad]+1;
                aux2=1;}
               
               }}
}
      
IEDmedia = (media/(disp-8-1));
firing_media = (1/IEDmedia);

for(n=8; n<disp; n++){
sum = sum + (IED[n]-IEDmedia)*(IED[n]-IEDmedia);
}

desv = sqrt( sum/(disp-10-1) );

CV = desv/IEDmedia;
sum = 0.0;



if(CV>0.5){

 if(quantidade_sharp[4]==quantidade_sharp[5] && quantidade_sharp[5]==quantidade_sharp[8] && quantidade_sharp[8]==quantidade_sharp[12]){
       regular = 1;}


if(regular==1){
escala_cv=3;}//regular bursting
else{
escala_cv=4;}//irregular bursting
}// disparo bursting
else{

for(n=4; n<disp; n++){
sum = sum + ((IED[n]-IED[n-1]) / (IED[n]+IED[n-1]));
}

A = (sum/(20-4-1));

   if(broad==0 || sharp==0){
      if(A>0.01){
        escala_cv=2; }//disparos adaptação
 
       if(A<0.01 && A>-0.01){
         escala_cv = 1;}//disparos tonicos
   }
   if(sharp!=0 && broad!=0){
   escala_cv=0;}//rajada inicial
  
 
 
}

fprintf(yyy," %lf \t %lf \t %i \t %lf \t %lf \t %lf \n", Vreset, b, escala_cv, A, CV, firing_media); 
fflush(yyy);
 for(n=1; n<disp; n++){
            IED[n] = 0;    
            }
            media = 0;
            IEDmedia = 0.0;
            sum = 0.0;
            desv = 0.0;
            CV = 0.0;
            A =0;
            disp = 0;
            broad =0;
            sharp=0;
            k=0;
            aux=0;
            aux2=0;
            verificando=0;
            regular = 0;
            mediafiring=0.0;
            firing_media=0.0;
    
            for(n=1; n<=10010000+1; n++){
            quantidade_sharp[n]=0;}
            
         //recolocando as condições iniciais   
         function[1] = -70.0;
         function[2] = 0.0;
}}




fclose(yyy);

//system("gnuplot espacoplot_irreg.plt");
free_dvector(IED, 0, 10010000+1);
free_dvector(firing_rate, 0, 10010000+1);
free_vector(quantidade_sharp, 0, 10010000+1);

return 0;}

double edo(double *df, int k, double variavel, double alpha, double beta){
  double eq[eq_num+1];
  
   eq[1] = alpha*pow(variavel, alpha-1)*(-(gl/C)*(df[1]-El)+(gl/C)*deltat*exp((df[1]-Vt)/deltat)- (df[2]/C) + (I/C)); //equação 1
   eq[2] = beta*pow(variavel, beta-1)*(((a/tw)*(df[1]-El))-(df[2]/tw)); //equação 2

   return eq[k];}
  
  
   ////////////////////////////////////////////////////////////////////////////////////////////
 
 double ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if(*idum<=0 || !iy)
    {
        if(-(*idum)<1) *idum=1;
        else *idum = -(*idum);
        for(j=NTAB+7;j>=0;j--)
        {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if(*idum<0) *idum +=IM;
            if(j<NTAB) iv[j]=*idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if(*idum<0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j]=*idum;
    if((temp=AM*iy)>RNMX) return RNMX;
    else return temp;
}
 
 
 void nrerror(char error_text[])
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

int *vector(long nl,long nh)
{
    int *v;

    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+NR_END;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int **m;

    m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
    if (!m) nrerror("allocation failure 1 in imatrix()");
    m += NR_END;
    m -= nrl;

    m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
    if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}


void free_vector(int *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}


void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}



  
  
