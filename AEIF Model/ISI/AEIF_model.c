//Leaky integrate-and-fire model
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define n_n 10000 // iteration 

//******   System Parameters ********************//
#define N 2         // Amount of equations
#define h 0.05        // integration step
#define C 200.0 //pF
#define gL 12.0 //nS
#define VT -50.0 // mV
#define DeltaT 2.0 //mV
#define tauw 300.0 // ms
#define EL -70.0 // mV
#define aa 2.0 //

void derivs(double *y, double *df, double Iext);

FILE *AEIF;
int main(void){
  double y[N+1], df[N+1], x[N+1], a[N+1], b[N+1], c[N+1], time, V_reset, bb, Iext;
  int i, t;

  AEIF = fopen("AEIF_model.dat", "wt");

  ////Initial conditions/////
  x[1] = -70.0; /// V
  x[2] = 0.0; /// w

  V_reset = -68.0;
  bb = 60.0;
  Iext = 0.0;
  time = 0.0;

  for (t = 1; t <= n_n; t++){ //iteration looping

    y[1] = x[1];
    y[2] = x[2];

    //if(time > 40.0)
      Iext = 509.0; //pA-------- Injecting the current
    
    //if(time > 350.0)
      //Iext = 0.0;


    fprintf(AEIF, "%.3f \t %.3f \t %.3f \n", time, x[1], x[2]);
    time = time + h; 

    // ------------ 4th order Runge-Kutta---------------
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++){
      a[i] = h * df[i];
      y[i] = x[i] + a[i] / 2.0;
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++){
      b[i] = h * df[i];
      y[i] = x[i] + b[i] / 2.0;
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++){
      c[i] = h * df[i];
      y[i] = x[i] + c[i];
    }
    derivs(y, df, Iext);
    for (i = 1; i <= N; i++)
      x[i] = x[i] + (a[i] + h * df[i]) / 6.0 + (b[i] + c[i]) / 3.0;
    //--------------------------------------------------------------------

    if(x[1]>-40.0){
      x[1] = V_reset;
      x[2] = x[2] + bb;
      fprintf(AEIF, "%.3f \t %.3f \t %.3f \n", time, 0.0, x[2]);
    }
    
  
  }

  fclose(AEIF);

  return 0;
}

void derivs(double *y, double *df, double Iext){ //Differential equations
  int i;

  //////////////////////////Coupled Equations////////////////////

    df[1] = (1.0 / C) * (-gL * (y[1] - EL) + gL * DeltaT * exp((y[1] - VT) / DeltaT) - y[2] + Iext);
    df[2] = (1.0 / tauw) * (aa * (y[1] - EL) - y[2]);
  
}
