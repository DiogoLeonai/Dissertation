#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//****** Parametros do gerador aleatorio ********************//
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define PI acos(-1.0)

#define n_n 600000

//******   System parameters ********************//
#define NN 17324 // 86620      // Amount of neuros
#define N 2           // Amount of equations
#define h 0.05        // passo de integracao
#define transient 25000 // ms
#define NMD 1000      // max number of firing of a neuron = 10^4
#define conexMAX 350  // each neuron max conections number
#define P_exc_ini 1.0 // Percentage of excitatory connections

//*******Random generator******//
float ran1(long *idum);
#define NR_END 1
#define FREE_ARG char *

//*******Pointers*****///
void nrerror(char error_text[]);
int *vector(long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_vector(int *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void derivs(double *y, double *df, double *I, double *aa, double *Iext, double *Gexc, double *Gini);

FILE *q, *raster, *Zt, *LocalZ;


int main(void){
  double *df, *y, *I, *a, *b, *c, *x, *aa, *tpeak, *Iext, *cont_Iext, *Gexc, *Gini, tauex, tauin;

  int i, j, cont, t, n, contISI;
  int **listaexc, *conextotalex, **listaini, *conextotalin;
  double tempo, b_IF, g_exc, g_ini, V_reset;
  double kinicial, kfinal, R, real, compl, Rmedio, **nk, *phi;
  int contR, *k, *kmax;
  double ISI, ISI2, CV, desviopadrao, radius, F, g, Rmediolocalmedio, FracaodeR60p, FracaodeR70p, FracaodeR80p, Rmediolocalmedio_semborda;
  double Rlocal[625], Rmediolocal[625], locreal[625], loccompl[625]; 
  int bin, auxm, auxn, loccont[625];                    

  int aux, aux2, aux3, aux4;
  double Deltax, Deltay, Deltaz, xi, yi, zi, x1[NN + 1], y1[NN + 1], z1[NN + 1], ymax, Distance;

  long idum;

  I = dvector(1, NN + 1);
  Iext = dvector(1, NN + 1);
  cont_Iext = dvector(1, NN + 1);
  aa = dvector(1, NN + 1);

  y = dvector(1, N * NN + 1);
  df = dvector(1, N * NN + 1);
  x = dvector(1, N * NN + 1);
  a = dvector(1, N * NN + 1);
  b = dvector(1, N * NN + 1);
  c = dvector(1, N * NN + 1);
  k = vector(1, NN + 1);
  kmax = vector(1, NN + 1);
  phi = dvector(1, NN + 1);
  nk = dmatrix(1, NMD + 2, 1, NN + 2);
  tpeak = dvector(1, NN + 1);
  Gexc = dvector(1, NN + 1);
  Gini = dvector(1, NN + 1);

  conextotalex = vector(1, NN + 2);
  listaexc = imatrix(1, NN + 1, 1, conexMAX + 2);

  conextotalin = vector(1, NN + 2);
  listaini = imatrix(1, NN + 1, 1, conexMAX + 2);

  char output_filename[100], output_filename2[100];

  //run this simulation for the following combinations of (radius, g_exc); (20, 0.01) and (70, 0.05). These combinations reproduce the results from the dissertation, Global order parameter and the raster and also their local Kuramoto order parameter
  radius = 20.0;
  g_exc = 0.01;


  idum = -1234567899;



  sprintf(output_filename2, "POK_CV_F_N17324_R%.1f_gexc_g0_a2_I500_100pontos.dat", radius);
  q = fopen(output_filename2, "wt");

  sprintf(output_filename, "Raster_R%.1f_gexc_%.3f.dat", radius, g_exc);
  raster = fopen(output_filename, "w");

  sprintf(output_filename2, "Z(t)_R%.1f_gexc_%.3f.dat", radius, g_exc);
  Zt = fopen(output_filename2, "w");

  sprintf(output_filename, "LocalZ_R%.1f_gexc_%.3f.dat", radius, g_exc);
  LocalZ = fopen(output_filename, "w");
  
  ymax = 976.0;
  Deltax = 7.0;
  Deltay = 8.0;
  Deltaz = 10.0;
  xi = 3.5;
  yi = 4.0;
  zi = 105.0;
  x1[1] = xi;
  y1[1] = yi;
  z1[1] = zi;
  aux = 1; // layer
  i = 1;


  for (i = 2; i <= NN; i = i + 1){
    if (yi + Deltay < ymax) /// 2nd y-axis
      yi = yi + Deltay;
    else /// 3th x-axis
    {
      xi = xi + Deltax;
      yi = 4.0;
    }
    zi = 105.0;
    aux = 1;

    x1[i] = xi;
    y1[i] = yi;

    z1[i] = zi;

  }

  //****************** Creating connections *********************
  //-----------zerando tudo------------//
  for (i = 1; i <= NN; i++)
  {
    conextotalex[i] = 0.0;
    conextotalin[i] = 0.0;
    for (j = 1; j <= conexMAX; j++)
    {
      listaexc[i][j] = 0.0; //// listaexc[pre-syn][conection number] = post-syn
      listaini[i][j] = 0.0;
    }
  }
  //------- creating a list with the connections--------//
  aux = 0;
  for (i = 1; i <= NN; i++)
  {
    //~ aux=0;
    for (j = 1; j <= NN; j++)
      if (i != j)
      {
        Distance = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]) + (z1[i] - z1[j]) * (z1[i] - z1[j])); // micrometers
        if (Distance <= radius)
        {
          conextotalex[i] = conextotalex[i] + 1;
          listaexc[i][conextotalex[i]] = j;
          aux = aux + 1;
        }
      }

  }

  //*****************************************************************
  idum = -1234567990;
  n = n_n;     // total number of steps n*h
  tauex = 1.5; // ms synaptic conductances time constants
  tauin = tauex;
  g = 0.0;
  //****************  Condicoes iniciais  **************************
  for (i = 1; i <= NN; i++){
    x[1 + (i - 1) * N] = -75.0 + 30.0 * ran1(&idum); // V
    x[2 + (i - 1) * N] = 0.0 + 70.0 * ran1(&idum);   // w
    Gexc[i] = 20.0;
    Gini[i] = -80.0;
    tpeak[i] = -100.0; // time of the last action potential
    aa[i] = 2.0;
    I[i] = 500.0; // External current
  }

    g_ini = g * g_exc;

    V_reset = -58.0;
    b_IF = 70.0;

    ISI = 0;
    ISI2 = 0;
    contISI = 0;
    F = 0.0;

    //********************* LOOP DO TEMPO  ***************************
    for (i = 1; i <= NN; i++){
      tpeak[i] = tpeak[i] - n * h; // time of the last action potential
      k[i] = 0.0;                  // cont the number of action potential of a neuron
      kmax[i] = 0.0;               // count the total amount of action potential of a neuron
      nk[1][i] = 0.0;              // time which a action potential happen in a neuron
      Iext[i] = 0.0;
      cont_Iext[i] = 0.0;
    }
    tempo = 0.0;

    for (t = 1; t <= n; t++){
      tempo = tempo + h; // ms

      for (i = 1; i <= NN; i++){ 
      
        Gexc[i] = Gexc[i] * exp(-h / tauex); 
        Gini[i] = Gini[i] * exp(-h / tauin);

        if (x[1 + (i - 1) * N] > -40.0) // AEIF model
        {
          if(tempo > transient)
            fprintf(raster, "%.3f \t %i \n", tempo, i);

          if (tempo > transient && k[i] < NMD){ 
          
            k[i] = k[i] + 1;
            nk[k[i]][i] = tempo;

            if (tempo < transient + 500.0)
              if (k[i] == 1.0)
                kinicial = tempo;
          }

          if (tempo > transient)
          {
            ISI = ISI + tempo - tpeak[i];
            ISI2 = ISI2 + (tempo - tpeak[i]) * (tempo - tpeak[i]);
            contISI = contISI + 1;
          }

          tpeak[i] = tempo;

          if (k[i] >= NMD)
            t = n + 1;

          x[1 + (i - 1) * N] = V_reset;                   // V
          x[2 + (i - 1) * N] = x[2 + (i - 1) * N] + b_IF; // w
        }

        y[1 + (i - 1) * N] = x[1 + (i - 1) * N];
        y[2 + (i - 1) * N] = x[2 + (i - 1) * N];

        if (tempo > tpeak[i] && tempo <= tpeak[i] + h) // delay=h
        {
          if (i <= P_exc_ini * NN)
            for (j = 1; j <= conextotalex[i]; j++)
              Gexc[listaexc[i][j]] = Gexc[listaexc[i][j]] + g_exc; // exc conductance received
          else
            for (j = 1; j <= conextotalin[i]; j++)
              Gini[listaini[i][j]] = Gini[listaini[i][j]] + g_ini;
        }
      }

      // ------------Runge-Kutta 4ª order---------------
      derivs(y, df, I, aa, Iext, Gexc, Gini);
      for (i = 1; i <= N * NN; i++)
      {
        a[i] = h * df[i];
        y[i] = x[i] + a[i] / 2.0;
      }
      derivs(y, df, I, aa, Iext, Gexc, Gini);
      for (i = 1; i <= N * NN; i++)
      {
        b[i] = h * df[i];
        y[i] = x[i] + b[i] / 2.0;
      }
      derivs(y, df, I, aa, Iext, Gexc, Gini);
      for (i = 1; i <= N * NN; i++)
      {
        c[i] = h * df[i];
        y[i] = x[i] + c[i];
      }
      derivs(y, df, I, aa, Iext, Gexc, Gini);
      for (i = 1; i <= N * NN; i++)
        x[i] = x[i] + (a[i] + h * df[i]) / 6.0 + (b[i] + c[i]) / 3.0;
      //--------------------------------------------------------------------

    } // end time loop

    // Calculating Kuramoto order parameter //
    tempo = 0.0;
    Rmedio = 0.0;
    contR = 0.0;
    kfinal = n * h;

    for (auxm = 0; auxm < 25; auxm++)
      for (auxn = 0; auxn < 25; auxn++)
        Rmediolocal[25 * auxn + auxm] = 0.0;

    for (i = 1; i <= NN; i = i + 1)
    {
      phi[i] = 0.0;
      kmax[i] = k[i];
      k[i] = 1.0;

      if (kmax[i] > 1)
        if (kfinal > nk[kmax[i]][i] && kfinal > kinicial)
          kfinal = nk[kmax[i]][i];
    }

    for (tempo = transient; tempo <= n * h; tempo = tempo + 20.0 * h)
    {
      for (i = 1; i <= NN; i = i + 1)
        if (kmax[i] > 1)
          if (tempo > nk[k[i]][i] && tempo < nk[kmax[i]][i])
          {
            if (tempo <= nk[k[i] + 1][i])
              phi[i] = 2 * PI * (tempo - nk[k[i]][i]) / (nk[k[i] + 1][i] - nk[k[i]][i]); //+2*PI*k[i]

            if (tempo >= nk[k[i] + 1][i])
              k[i] = k[i] + 1;
          }

      if (tempo >= kinicial && tempo <= kfinal) // se o tempo for maior que kinicial e menor que kfinal
      {
        real = 0.0;
        compl = 0.0;
        cont = 0;

        for (bin = 0; bin < 25 * 25; bin++)
        {
          locreal[bin] = 0.0;
          loccompl[bin] = 0.0;
          Rlocal[bin] = 0.0;
          loccont[bin] = 0;
        }

        for (i = 1; i <= NN; i++) // network loop
          if (kmax[i] > 1)
          {
            real = real + cos(phi[i]); // i is the index of the neuron
            compl = compl +sin(phi[i]);
            cont = cont + 1;

            // Local order parameter
            auxn = ceil(y1[i] / 40.0);
            auxm = ceil(x1[i] / 40.0);
            bin = 25 * (auxn - 1) + (auxm - 1);
            locreal[bin] = locreal[bin] + cos(phi[i]);
            loccompl[bin] = loccompl[bin] + sin(phi[i]);
            loccont[bin] = loccont[bin] + 1;
          }

        real = real / cont;   // average of cos(phi[i]) over the network
        compl = compl / cont; // time average of sin(phi[i]) over the network

        R = sqrt(real * real + compl *compl ); // Global order parameter
        Rmedio = Rmedio + R;
        contR = contR + 1.0;

         fprintf(Zt,"%f \t %f \n",tempo,R);

        for (bin = 0; bin < 25 * 25; bin++){ // Local order parameter
        
          locreal[bin] = locreal[bin] / loccont[bin];   // average of cos(phi[i]) over the bin
          loccompl[bin] = loccompl[bin] / loccont[bin]; // average of sin(phi[i]) over the bin
          Rlocal[bin] = sqrt(locreal[bin] * locreal[bin] + loccompl[bin] * loccompl[bin]);
          // fprintf(rout2, "%.3f ", Rlocal[bin]);
          Rmediolocal[bin] = Rmediolocal[bin] + Rlocal[bin];
        }
      // fprintf(rout2, "\n");
      }
    }

    // Outputs local order parameter
    for (bin = 0; bin < 25 * 25; bin++)
      Rmediolocal[bin] = Rmediolocal[bin] / contR;

       for(i=1; i <= NN; i++){
          auxn = ceil(y1[i] / 40.0);
          auxm = ceil(x1[i] / 40.0);
         fprintf(LocalZ, "%f %f %.4f \n", x1[i], y1[i], Rmediolocal[25 * (auxn - 1) + (auxm - 1)]);
         }


    desviopadrao = sqrt(ISI2 / contISI - (ISI / contISI) * (ISI / contISI));

    if (contISI > 0)
    {
      CV = desviopadrao / (ISI / contISI);
      F = ((1000.0) / (ISI / contISI));
    }
    else
    {
      CV = 0.0;
      F = 0;
    }

    F = ((1000.0) / (ISI / contISI));
    if (contISI == 0)
      F = 0;

    ///////////////////////////////////////////////////////////
    // CALCULO DO Rmediolocal medio
    real = 0.0;
    cont = 0;
    aux2 = 0;
    aux3 = 0;
    aux4 = 0;
    for (auxm = 0; auxm < 25; auxm++)
      for (auxn = 0; auxn < 25; auxn++)
      {
        real = real + Rmediolocal[25 * auxn  + auxm ];
        cont = cont + 1;
        if (Rmediolocal[25 * auxn  + auxm ] < 0.7)
          aux3 = aux3 + 1;

      }
    ///////////////////////////////////////////////////////////
    Rmediolocalmedio = 1.0*real / cont;
    FracaodeR70p = 1.0*aux3;


    ///////////////////////////////////////////////////////////
    Rmediolocalmedio_semborda = 1.0*real / cont;


    printf("%.2f \t %.4f \t %.3f \t %.3f \t %.3f \t %.3f\n", radius, g_exc, g, Rmedio / contR, CV, Rmediolocalmedio); // aux : num total conex
    fprintf(q, "%.2f \t %.4f \t %.3f \t %.3f \t %.3f \t %.3f\n", radius, g_exc, g, Rmedio / contR, CV, Rmediolocalmedio);



    fclose(Zt);
    fclose(LocalZ);
    fclose(raster);
    fclose(q);


  
  free_dvector(y, 1, N * NN + 1);
  free_dvector(df, 1, N * NN + 1);
  free_dvector(x, 1, N * NN + 1);
  free_dvector(a, 1, N * NN + 1);
  free_dvector(b, 1, N * NN + 1);
  free_dvector(c, 1, N * NN + 1);
  free_dvector(I, 1, NN + 1);
  free_dvector(Iext, 1, NN + 1);
  free_dvector(cont_Iext, 1, NN + 1);
  free_dvector(aa, 1, NN + 1);
  free_vector(k, 1, NN + 1);
  free_vector(kmax, 1, NN + 1);
  free_dvector(phi, 1, NN + 1);
  free_dmatrix(nk, 1, NMD + 2, 1, NN + 2);
  free_dvector(tpeak, 1, NN + 1);
  free_dvector(Gexc, 1, NN + 1);
  free_dvector(Gini, 1, NN + 1);

  free_vector(conextotalex, 1, NN + 2);
  free_imatrix(listaexc, 1, NN + 1, 1, conexMAX + 2);
  free_vector(conextotalin, 1, NN + 2);
  free_imatrix(listaini, 1, NN + 1, 1, conexMAX + 2);


  return 0;
}

void derivs(double *y, double *df, double *I, double *aa, double *Iext, double *Gexc, double *Gini) // Equacoes diferenciais acopladas
{
  int i;
  double C, gL, VT, DeltaT, tauw, Vr_exc, Vr_ini, EL;

  C = 200.0;    // pF
  gL = 12.0;    // nS
  VT = -50.0;   // mV
  DeltaT = 2.0; // mV
  tauw = 300.0; // ms
  EL = -70.0;
  Vr_exc = 0.0;   // mV Potencial reverso exc
  Vr_ini = -80.0; // mV Potencial reverso inh

  //////////////////////////EQUAÇOES ACOPLADAS////////////////////*(Vr_exc-y[1+(i-1)*N])

  for (i = 1; i <= NN; i++)
  {
    df[1 + (i - 1) * N] = (1.0 / C) * (-gL * (y[1 + (i - 1) * N] - EL) + gL * DeltaT * exp((y[1 + (i - 1) * N] - VT) / DeltaT) - y[2 + (i - 1) * N] + Iext[i] + Gexc[i] * (Vr_exc - y[1 + (i - 1) * N]) + Gini[i] * (Vr_ini - y[1 + (i - 1) * N]) + I[i]);

    df[2 + (i - 1) * N] = (1.0 / tauw) * (aa[i] * (y[1 + (i - 1) * N] - EL) - y[2 + (i - 1) * N]);
  }
}

float ran1(long *idum)
{
  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy)
  {
    if (-(*idum) < 1)
      *idum = 1;
    else
      *idum = -(*idum);
    for (j = NTAB + 7; j >= 0; j--)
    {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0)
        *idum += IM;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0)
    *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM * iy) > RNMX)
    return RNMX;
  else
    return temp;
}

double *dvector(long nl, long nh)
{
  double *v;

  v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v - nl + NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

int *vector(long nl, long nh)
{
  int *v;

  v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
  if (!v)
    nrerror("allocation failure in dvector()");
  return v - nl + NR_END;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  m = (int **)malloc((size_t)((nrow + NR_END) * sizeof(int *)));
  if (!m)
    nrerror("allocation failure 1 in imatrix()");
  m += NR_END;
  m -= nrl;

  m[nrl] = (int *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
  if (!m[nrl])
    nrerror("allocation failure 2 in imatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  return m;
}

void free_vector(int *v, long nl, long nh)
{
  free((FREE_ARG)(v + nl - NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
  if (!m)
    nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl])
    nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl] + ncl - NR_END));
  free((FREE_ARG)(m + nrl - NR_END));
}
