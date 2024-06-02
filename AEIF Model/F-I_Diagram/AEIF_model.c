//Leaky integrate-and-fire model
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define n_n 30000 // iteration 
#define transient 500.0 //ms

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

FILE *FxI;
int main(void){
  double y[N+1], df[N+1], x[N+1], a[N+1], b[N+1], c[N+1], time, V_reset, bb, Iext, tpeak, meanISI, contISI, F;
  int i, t, colorScale;

  FxI = fopen("FxI_Diagram.dat", "wt");

  for(Iext = 0.0; Iext <= 1500.0; Iext = Iext + 75){
  ////Initial conditions/////
    x[1] = -70.0; /// V
    x[2] = 0.0; /// w
    
    contISI = 0.0;
    meanISI = 0.0;
    tpeak = 0.0; // time instant of the firing k-th
    V_reset = -68.0;
    bb = 60.0;
    time = 0.0;

    for (t = 1; t <= n_n; t++){ //iteration looping

      y[1] = x[1];
      y[2] = x[2];
      
      time = time + h; 

      if(x[1]>-40.0){
        x[1] = V_reset;
        x[2] = x[2] + bb;
        if(time > transient){
          meanISI = meanISI + (time - tpeak);
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
    
    F = ((1000.0) / (meanISI / contISI)); // Frequency

    if(contISI == 0){ //if the neuron does not fire the frequency is zero
      F = 0.0;
      colorScale = 0;
      }

    if(F>0){

      if( F<4)
        colorScale = 1; //delta frequency
      else{
        if(F<8)
          colorScale = 2; //theta frequency
        else{
          if(F<12)
            colorScale = 3; //alpha frequency
          else{
            if(F<30)
              colorScale = 4; //beta frequency
            else{
              colorScale = 5; //gamma frequency
            }
          }
        }
      }
    }
    fprintf(FxI, "%.3f \t %.3f \t %i \n", Iext, F, colorScale);
    
  }//end loop Iext

  fclose(FxI);

  return 0;
}

void derivs(double *y, double *df, double Iext){ //Differential equations
  int i;

  //////////////////////////Coupled Equations////////////////////

    df[1] = (1.0 / C) * (-gL * (y[1] - EL) + gL * DeltaT * exp((y[1] - VT) / DeltaT) - y[2] + Iext);
    df[2] = (1.0 / tauw) * (aa * (y[1] - EL) - y[2]);
  
}
