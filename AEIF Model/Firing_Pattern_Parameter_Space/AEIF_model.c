//Leaky integrate-and-fire model
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define n_n 100000 // iteration 
#define transient 1000.0 //ms
#define NMF 1000 // number maximum of firing

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
  double y[N+1], df[N+1], x[N+1], a[N+1], b[N+1], c[N+1], time, V_reset, bb, Iext, tpeak, meanISI, meanISI2, contISI, F, standard_deviation, CV, ISI[NMF + 1], A, At;
  double dv, dw;
  int i, t, Region, m, d, f, broad, aux, aux2, k, sharp, amount_sharp[NMF + 1];

  SpikeBurstRegion = fopen("Bursting_Region.dat", "wt");
  Iext = 512.00;
  d = 5; // interspike interval initial which we calculate A
  f = 40; // interspike interval final which we calculate A
  At = 0.01; // threshold of the adaptive index

  for(V_reset = -60.0; V_reset <= -45.0; V_reset = V_reset + 0.15){
    for(bb = 0.0; bb <= 350.0; bb= bb + 3.5){
    ////Initial conditions/////
    x[1] = -70.0; /// V
    x[2] = 0.0; /// w
    
    contISI = 0.0;
    meanISI = 0.0;
    meanISI2 = 0.0;
    tpeak = 0.0; // time instant of the firing k-th
    time = 0.0;
    CV = 0.0;
    A = 0.0;
    sharp = 0;
    broad = 0;
    aux = 0;
    aux2 = 0;
    k=0;

    for(i = 0; i <= NMF; i++)
      amount_sharp[i] = 0;

    for(m = d; m <=f; m++)
      ISI[m] = 0.0;

    m = 0;
    for (t = 1; t <= n_n; t++){ //iteration looping


      
      
      time = time + h; 

      if(x[1]>-40.0){
        x[1] = V_reset;
        x[2] = x[2] + bb;
        m++;
        if(m > NMF)
          t = n_n + 1;

        if(m>=4)
          ISI[m] = time - tpeak; //calculating each ISI

        if(time > transient){
          meanISI = meanISI + (time - tpeak);
          meanISI2 = meanISI2 + (time - tpeak) * (time - tpeak);
          contISI = contISI + 1.0;
        }
        tpeak = time;

        dv = (1.0 / C) * (-gL * (x[1] - EL) + gL * DeltaT * exp((y[1] - VT) / DeltaT) - x[2] + Iext); //dv/dt

        if(dv < 0){
          broad++;

          if(aux2==1){
            k++;
           aux2=0;
          }

          aux = 1;
        }             
        else{
          sharp++;

          if(aux==1){
            amount_sharp[broad] = amount_sharp[broad] + 1;
            aux2=1;
          }
        
        }
      }
      y[1] = x[1];
      y[2] = x[2];
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

    for(m = d; m <= f; m++)
      A = A + (ISI[m] - ISI[m-1]) / (ISI[m] + ISI[m-1]);

    A = A / (f-d-1);

    if (contISI > 0){
        CV = standard_deviation / (meanISI / contISI);
      }
    
    if(CV < 0.5){

      if(A > At){
        Region = 1; //adaptation spike
      }
      else{
        Region = 2; // tonic spike
        }
      
      if(sharp!=0 && broad!=0){
        Region=3; // initial bursting
      }
      }
    else{

      for(i = 1; i <= NMF; i++){
        Region = 4; // Regular bursting

        if(amount_sharp[i] != amount_sharp[0]){
          i = NMF + 1;
          Region = 5; // irregular bursting
        }

      }

    }

    fprintf(SpikeBurstRegion, "%.3f \t %.3f \t %i\n", V_reset, bb, Region);
    //printf("%.3f \t %.3f \t %i \t %.4f\n", V_reset, bb, Region, A);
    
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
