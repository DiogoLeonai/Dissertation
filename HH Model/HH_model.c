//Hodgkin-Huxley model// HODGKIN, A. L.; HUXLEY, A. F. A quantitative description of membrane current and its application to conduction and excitation in nerve. Journal of Physiology, v. 117, p. 500–544, 1952.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define n_n 10000 // iteration 

//******   System Parameters ********************//
#define N 4         // Amount of equations
#define step 0.05        // integration step

void derivs(double *y, double *df, double Iext);

FILE *HH;
int main(void){
  double y[N+1], df[N+1], x[N+1], a[N+1], b[N+1], c[N+1], time, Iext;
  int i, t;

  HH = fopen("HH_model.dat", "wt");

  ////Initial conditions/////
  x[1] = 0.0; //V 
  x[2] = 0.0; // n
  x[3] = 0.0; // m
  x[4] = 0.0; // h
  //////////////////////////

  Iext = 0.0;
  time = 0.0;

  for (t = 1; t <= n_n; t++){ //iteration looping

    y[1] = x[1];
    Iext = 0.0;

    if(time > 140.0 && time < 142.0)
      Iext = 3.0; // mu A/cm^2

    if(time > 170.0 && time < 172.0)
      Iext = 8.0; // mu A/cm^2

    if(time > 182.0 && time < 184.0)
      Iext = 8.0; // mu A/cm^2

    if(time > 230.0 && time < 232.0)
      Iext = 8.0; // mu A/cm^2 
    
      

    fprintf(HH, "%.3f \t %.1f \t", time, Iext);
    for(i = 1; i <= N; i++){
      fprintf(HH, "%.3f \t", x[i]);
    }
    fprintf(HH, "\n");

    time = time + step; 

    // ------------ 4th order Runge-Kutta---------------
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
    {
      a[i] = step * df[i];
      y[i] = x[i] + a[i] / 2.0;
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
    {
      b[i] = step * df[i];
      y[i] = x[i] + b[i] / 2.0;
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
    {
      c[i] = step * df[i];
      y[i] = x[i] + c[i];
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
      x[i] = x[i] + (a[i] + step * df[i]) / 6.0 + (b[i] + c[i]) / 3.0;
    //--------------------------------------------------------------------
  
  }

  fclose(HH);

  return 0;
}

void derivs(double *y, double *df, double Iext){ //Differential equations
  int i;
  double V, n, m, h;
  double C, gl, El, gk, Ek, gna, Ena, alpha_n, beta_n, alpha_m, beta_m, alpha_h, beta_h;

  ////////Parameters/////
  C = 1.0;    // pF
  gl = 12.0;    // nS
  El = 10.6; //mV
  Ena = 120.0; //mV
  Ek = -12.0;
  gna = 120.0;
  gk = 36.0;
  gl = 0.3;


  V = y[1];
  n = y[2];
  m = y[3];
  h = y[4];

  //Alpha and beta parameters of the squid
  alpha_n = 0.01 * (10.0 - V) / (exp( 0.1*(10.0 - V) ) - 1.0);
  beta_n = 0.125 * exp(- V / 80.0);
  alpha_m = 0.1 * (25.0 - V) / (exp( 0.1*(25.0 - V) ) -1);
  beta_m = 4.0 * exp(-V/18.0);
  alpha_h = 0.07 * exp(-V/20.0);
  beta_h = 1 / (exp( 0.1*(30 - V)) + 1);

  //////////////////////////Diferential equations////////////////////
  df[1] = (1/C) * ( -gl * (V-El) - gk * pow(n, 4) * (V - Ek) - gna * pow(m, 3) * h * (V -Ena) + Iext); //dv/dt
  df[2] = alpha_n * (1 - n) - beta_n * n; //dn/dt
  df[3] = alpha_m * (1 - m) - beta_m * m; // dm/dt
  df[4] = alpha_h * (1 - h) - beta_h * h; // dh/dt
  
}
