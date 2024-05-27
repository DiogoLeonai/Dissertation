//Leaky integrate-and-fire model
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define n_n 10000 // iteration 

//******   System Parameters ********************//
#define N 1         // Amount of equations
#define h 0.05        // integration step

void derivs(double *y, double *df, double Iext);

FILE *lif;
int main(void){
  double y[N+1], df[N+1], x[N+1], a[N+1], b[N+1], c[N+1], time, V_reset, Iext;
  int i, t;

  lif = fopen("LIF.dat", "wt");

  ////Initial conditions/////
  x[1] = -70.0; 

  V_reset = -70.0;
  Iext = 0.0;
  time = 0.0;

  for (t = 1; t <= n_n; t++){ //iteration looping

    y[1] = x[1];

    if(time > 40.0)
      Iext = 509.0; //pA
    
    if(time > 350.0)
      Iext = 0.0;

    fprintf(lif, "%lf \t %lf \n", time, x[1]);

    time = time + h; 

    // ------------ 4th order Runge-Kutta---------------
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
    {
      a[i] = h * df[i];
      y[i] = x[i] + a[i] / 2.0;
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
    {
      b[i] = h * df[i];
      y[i] = x[i] + b[i] / 2.0;
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
    {
      c[i] = h * df[i];
      y[i] = x[i] + c[i];
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
      x[i] = x[i] + (a[i] + h * df[i]) / 6.0 + (b[i] + c[i]) / 3.0;
    //--------------------------------------------------------------------

    if(x[1]>-40.0){
      x[1] = V_reset;
      fprintf(lif, "%lf \t %lf \n", time, -20.0);
    }
    
  
  }

  fclose(lif);

  return 0;
}

void derivs(double *y, double *df, double Iext){ //Differential equations
  int i;
  double C, gL, VT, EL;

  C = 280.0;    // pF
  gL = 12.0;    // nS
  VT = -50.0;   // mV
  EL = -70.0; //mV


  //////////////////////////Equations////////////////////
  df[1] = (1.0 / C) * (-gL * (y[1] - EL) + Iext); //   dV/dt
  
}
