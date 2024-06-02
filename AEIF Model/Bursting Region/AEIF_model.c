//Leaky integrate-and-fire model
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define n_n 100000 // iteration 
#define transient 1000.0 //ms

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

FILE *SpikeBurstRegion;
int main(void){
  double y[N+1], df[N+1], x[N+1], a[N+1], b[N+1], c[N+1], time, V_reset, bb, Iext, tpeak, meanISI, meanISI2, contISI, F, standard_deviation, CV;
  int i, t, Region;

  SpikeBurstRegion = fopen("Bursting_Region.dat", "wt");
  Iext = 500.00;

  for(V_reset = -60.0; V_reset <= -45; V_reset = V_reset + 0.015){
    for(bb = 0.0; bb <= 250.0; bb= bb + 0.25){
  ////Initial conditions/////
    x[1] = -70.0; /// V
    x[2] = 0.0; /// w
    
    contISI = 0.0;
    meanISI = 0.0;
    meanISI2 = 0.0;
    tpeak = 0.0; // time instant of the firing k-th
    time = 0.0;
    CV = 0.0;

    for (t = 1; t <= n_n; t++){ //iteration looping

      y[1] = x[1];
      y[2] = x[2];
      
      time = time + h; 

      if(x[1]>-40.0){
        x[1] = V_reset;
        x[2] = x[2] + bb;
        if(time > transient){
          meanISI = meanISI + (time - tpeak);
          meanISI2 = meanISI2 + (time - tpeak) * (time - tpeak);
          contISI = contISI + 1.0;
        }
        tpeak = time;
      }

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
    } //end temporal loop
    
    CV = 0.0;
    standard_deviation = 0.0;
    standard_deviation = sqrt(meanISI2 / contISI - (meanISI / contISI) * (meanISI / contISI));

    if (contISI > 0){
        CV = standard_deviation / (meanISI / contISI);
      }
    
    if(CV < 0.5){
      Region = 0; // Spike region
      }
    else{
      Region = 1; // Burst region
    }

    fprintf(SpikeBurstRegion, "%.3f \t %.3f \t %.3f \t %i\n", V_reset, bb, CV, Region);
    
  }//end loop bb

  }//end loop V_reset

  fclose(SpikeBurstRegion);

  return 0;
}

void derivs(double *y, double *df, double Iext){ //Differential equations
  int i;

  //////////////////////////Coupled Equations////////////////////

    df[1] = (1.0 / C) * (-gL * (y[1] - EL) + gL * DeltaT * exp((y[1] - VT) / DeltaT) - y[2] + Iext);
    df[2] = (1.0 / tauw) * (aa * (y[1] - EL) - y[2]);
  
}
